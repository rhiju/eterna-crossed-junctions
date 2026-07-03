# CLAUDE.md — project context for Claude Code

You are working in **eterna-crossed-junctions**, a prototype for a new Eterna game
constraint. Read this first, then `README.md`. If you're a colleague or Eterna player
picking this up: just `git clone` the repo, open it with Claude Code (or `claude` in
the repo dir), and this file gives Claude the context to help you.

## What this is

A **"crossed-junction" constraint** for RNA secondary structures: every multiway
junction (interior multiloop, or the exterior loop when it joins ≥3 stems) should be
**crossed** — touched by a pseudoknot. Motivation: a bare 3-/4-way junction with no
tertiary contact tends to be floppy and hard to resolve by cryo-EM. Idea from the
Das/Janelia cryo-EM group. Full rule and motivation are in `README.md`.

Two phases:
- **Phase 1 (done):** a self-contained Python prototype + verified examples + a review
  over real Eterna targets + draw_rna figures.
- **Phase 2 (draft):** a TypeScript port to ship as an EternaJS `Constraint`
  (`typescript/CrossedJunctionConstraint.ts`, mapping in `PORTING.md`).

## Layout

- `python/crossed_junctions.py` — core library (stdlib only). Start here. Key funcs:
  `parse_dotbracket`, `crossed_residues` (port of `SecStruct.getCrossedPairs`),
  `separate_layers` (stem-based, longest-first), `find_junctions`, `check`.
- `python/generate_examples.py` — writes/verifies `examples/*.dbn` (len 100, stems ≥4bp).
- `python/review_targets.py` — runs the check over the OpenKnotAIDesignData targets.
- `python/draw_crossed.py` — draw_rna figures (junctions green=crossed / red=bare, PKs as arcs).
- `typescript/`, `PORTING.md` — the EternaJS port.

## How to run (verify everything works)

```bash
cd python
python3 crossed_junctions.py     # unit sanity checks
python3 generate_examples.py     # regenerate + verify the 10 examples -> "ALL 10 EXAMPLES OK"
python3 review_targets.py        # pass/fail table over real targets (downloads CSVs on first run)
python3 draw_crossed.py          # figures for default real targets -> figures/
```

Figures need `matplotlib`, `numpy`, and `draw_rna`
(https://github.com/eternagame/draw_rna). The core check/generator/review are
stdlib-only. `python3` (not `python`).

## Conventions / gotchas (please preserve)

- **Layering is by whole stems, longest first.** A naive left-to-right greedy can pull
  a short pseudoknot into layer 0 and push a real backbone helix out, distorting the
  junctions. Don't "simplify" this back to per-pair greedy.
- Because `getCrossedPairs` marks **both** members of a crossed pair, a pseudoknot
  crossing a junction's own stem also credits that junction (its stem-boundary residues
  are junction members). This is intended.
- The constraint is **soft**: `check()` returns per-junction detail; the TS draft reports
  `crossedJunctions / totalJunctions`. The pass threshold (`minFraction`) is a tuning
  decision, not fixed.
- Keep `crossed_junctions.py` dependency-free so it ports cleanly to TypeScript.

## Giving feedback / iterating

If you have feedback from colleagues or Eterna players, describe the desired behavior
change to Claude (e.g. "a 3-way exterior loop shouldn't count", "require 50% of
junctions not all", "handle chain breaks `&`"). The natural places to change:
- the **rule** → `find_junctions` / `check` in `crossed_junctions.py`
  (keep the TS `findJunctions`/`evaluate` in sync);
- **new test structures** → add to `generate_examples.py` (they self-verify);
- **figures** → `draw_crossed.py`.
Re-run the four commands above to confirm nothing regressed before committing.

## Status

Repo: https://github.com/rhiju/eterna-crossed-junctions (public). Phase 1 complete and
self-verifying. Phase 2 (EternaJS constraint) is a reviewed draft awaiting wiring into
the game — see `PORTING.md` for the remaining steps.
