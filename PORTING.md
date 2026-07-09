# Porting the prototype into EternaJS

The Python prototype (`python/crossed_junctions.py`) and the TypeScript draft
(`typescript/CrossedJunctionConstraint.ts`) implement the same algorithm. This maps
each piece to its EternaJS counterpart so wiring it into the game is mechanical.

## Function map

| Python (`crossed_junctions.py`) | EternaJS |
| --- | --- |
| `parse_dotbracket(dbn)` | Not needed in-game — pairs come from `undoBlock.getPairs(temp, pseudoknots)` (a `SecStruct`). For tests, `SecStruct.fromParens(str, true)` / `EPars.parenthesisToPairs`. |
| `crossed_pairs` / `crossed_residues` | `SecStruct.getCrossedPairs()` already exists (same interleave test). Iterate its `pairingPartner(i)` to get the residue set. |
| `_cross(bp1, bp2)` | Inlined `crosses()` in the draft; identical to the test inside `SecStruct.getCrossedPairs()`. |
| `stems(pairs)` | `SecStruct.stems()` already exists. |
| `separate_layers(pairs)` | `separateLayers()` in the draft, which calls `SecStruct.getParenthesis({pseudoknots: true})` and splits by bracket char (same routine `filterForPseudoknots`/`onlyPseudoknots` use). This guarantees the layering matches the game. |
| `_immediate_children` / `find_junctions` | `immediateChildren()` / `findJunctions()` in the draft. **The only genuinely new logic.** |
| `check(dbn)` | `CrossedJunctionConstraint.evaluate(context)`. |

## Constraint scaffolding (copy from `MinimumCrossedPercentConstraint.ts`)

- `NAME` string + `serialize()` returning `[NAME, param]`.
- `evaluate(context)` reads `context.undoBlocks[0]`, detects `pseudoknot` target type,
  calls `undoBlock.getPairs(EPars.DEFAULT_TEMPERATURE, pseudoknots)`.
- `getConstraintBoxConfig(status)` builds the green/red stat text and returns a
  `ConstraintBoxConfig` (icon, tooltip). Uses `BitmapManager.getBitmap(...)`.

## Soft-constraint design

The draft reports `crossedJunctions / totalJunctions` and is satisfied when the
fraction ≥ `minFraction` (default 1 = all junctions crossed). Because the puzzle
wants a **soft** target, expose `minFraction` as the serialized parameter and let the
puzzle author tune it (e.g. 0.5). Structures with no ≥3-way junction are trivially
satisfied (fraction defined as 1).

## Remaining wiring (not done here)

1. Confirm the `evaluate` signature on the target branch (`dev` uses
   `ConstraintContext`; older branches pass `(undoBlocks, targetConditions, puzzle)`).
2. Register `CrossedJunctionConstraint` in the constraint factory that deserializes
   `[NAME, param]` (same place `MIN_CROSSED_PERCENT` is registered).
3. Add a bespoke constraint-box icon (currently reuses `PseudoknotReqIcon`).
4. Optional `getHighlight()` to highlight bare junctions in the pose.
5. Tests: port `python/crossed_junctions.py`'s `__main__` cases and the 10
   `examples/*.dbn` into the EternaJS test suite.

## Behavioral notes worth preserving

- Layering **must match Eterna's `getParenthesis`**: walk pairs left-to-right by opening
  index and place each in the lowest layer where it doesn't cross a pair already there
  (0 -> `()`, 1 -> `[]`, ...). It is per-pair and left-to-right, **not** stem-based or
  longest-first. A short pseudoknot that opens before a longer stem it crosses belongs
  in layer 0 — that is the game's behaviour. (An earlier longest-first heuristic diverged
  from Eterna and was the bug this fix addresses.)
- Because `getCrossedPairs()` marks **both** members of a crossed pair, a pseudoknot
  that crosses a junction's own stem also marks that stem's boundary residues — which
  are junction members. So crossing a junction's branch/closing helix satisfies it,
  not only landing a strand on an unpaired linker. This is intended.
