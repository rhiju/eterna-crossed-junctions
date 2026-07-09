"""
Crossed-junction constraint for RNA secondary structures.

Checks whether every multiway junction (interior multiloop, or the exterior loop
when it joins >=3 stems) is "crossed" -- i.e. touched by a pseudoknot base pair.

The idea (from the #janelia_rnacryoem Slack thread): a bare 3-/4-way junction with
no designed tertiary contact tends to be conformationally floppy and hard to
resolve by cryo-EM; a pseudoknot crossing the junction pins it in place.

Crossing rule (locked): a junction PASSES if *any* residue bordering it -- a base
in one of its closing/branching stem pairs, OR an unpaired linker base directly in
the junction -- participates in a crossed pair of the full structure.

Self-contained: stdlib only. Ports cleanly to the EternaJS TypeScript constraint.

Reference implementations mirrored here:
  - crossed-pair test:   EternaJS SecStruct.getCrossedPairs()
                         (== OpenKnotScore figure_out_which_bps_are_crossed.m)
  - multi-layer parsing: EternaJS SecStruct.setPairs(parens, pseudoknots=True)
  - layer separation:    arnie utils _group_into_non_conflicting_bp()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Set, Tuple

Pairs = List[int]   # pairs[i] = partner index of base i, or -1 if unpaired
BP = Tuple[int, int]


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

# Bracket layers. Uppercase letters open, lowercase close (Aa, Bb, ...), matching
# arnie's convert_bp_list_to_dotbracket convention. EternaJS auto-detects letter
# polarity from the first letter seen; real OpenKnot targets only use ()[]{}<>, so
# the two conventions agree on all inputs we care about.
_OPEN = "([{<"
_CLOSE = ")]}>"
_MATCH = {")": "(", "]": "[", "}": "{", ">": "<"}


def parse_dotbracket(dbn: str) -> Pairs:
    """Convert (possibly multi-layer) dot-bracket notation to a pairs array."""
    n = len(dbn)
    pairs: Pairs = [-1] * n
    stacks = {b: [] for b in _OPEN}          # bracket stacks
    letter_stacks: dict = {}                  # 'A' -> [positions]
    # Letter-pair polarity (which case opens) is auto-detected from the first
    # letter seen, matching EternaJS SecStruct.setPairs: some structures use
    # Aa (upper opens), others aA (lower opens).
    opener_upper = None

    for i, c in enumerate(dbn):
        if c in (".", "-", "&", ":", ","):
            continue
        if c in _OPEN:
            stacks[c].append(i)
        elif c in _CLOSE:
            st = stacks[_MATCH[c]]
            if not st:
                raise ValueError(f"Unbalanced notation: extra '{c}' at position {i}")
            j = st.pop()
            pairs[i] = j
            pairs[j] = i
        elif c.isalpha():
            if opener_upper is None:
                opener_upper = c.isupper()
            if c.isupper() == opener_upper:   # opener
                letter_stacks.setdefault(c.upper(), []).append(i)
            else:                             # closer
                st = letter_stacks.get(c.upper())
                if not st:
                    raise ValueError(f"Unbalanced notation: extra '{c}' at position {i}")
                j = st.pop()
                pairs[i] = j
                pairs[j] = i
        else:
            raise ValueError(f"Unknown character {c!r} at position {i}")

    for b, st in stacks.items():
        if st:
            raise ValueError(f"Unbalanced notation: unclosed '{b}'")
    for L, st in letter_stacks.items():
        if st:
            raise ValueError(f"Unbalanced notation: unclosed '{L}'")
    return pairs


def pairs_to_bp_list(pairs: Pairs) -> List[BP]:
    """Sorted list of (i, j) base pairs with i < j."""
    return sorted((i, pairs[i]) for i in range(len(pairs)) if pairs[i] > i)


# --------------------------------------------------------------------------- #
# Crossed-pair detection  (port of SecStruct.getCrossedPairs)
# --------------------------------------------------------------------------- #

def _cross(bp1: BP, bp2: BP) -> bool:
    """True if two base pairs are topologically interleaved (crossed)."""
    a, b = bp1
    c, d = bp2
    return (a < c < b < d) or (c < a < d < b)


def crossed_pairs(pairs: Pairs) -> Set[BP]:
    """Set of base pairs that cross at least one other pair."""
    bp_list = pairs_to_bp_list(pairs)
    crossed: Set[BP] = set()
    for i in range(len(bp_list)):
        for j in range(i + 1, len(bp_list)):
            if _cross(bp_list[i], bp_list[j]):
                crossed.add(bp_list[i])
                crossed.add(bp_list[j])
    return crossed


def crossed_residues(pairs: Pairs) -> Set[int]:
    """Set of residues that participate in any crossed base pair."""
    res: Set[int] = set()
    for a, b in crossed_pairs(pairs):
        res.add(a)
        res.add(b)
    return res


# --------------------------------------------------------------------------- #
# Layer separation  (greedy, mirrors arnie _group_into_non_conflicting_bp)
# --------------------------------------------------------------------------- #

def separate_layers(pairs: Pairs) -> List[Pairs]:
    """
    Split base pairs into non-crossing layers, matching EternaJS exactly.

    This is a faithful port of the layer/degree assignment in
    `SecStruct.getParenthesis({pseudoknots: true})`: walk pairs left-to-right by
    OPENING index and place each into the LOWEST-numbered layer where it does not
    cross a pair already placed in that layer (layer 0 -> '()', 1 -> '[]', 2 ->
    '{}', 3 -> '<>', 4+ -> letter pairs). Layer 0 is Eterna's "backbone"
    (== filterForPseudoknots); layers 1+ are the pseudoknots (== onlyPseudoknots).

    NOTE: this is per-pair and left-to-right, NOT stem-based or longest-first. A
    short pseudoknot stem that opens before a longer stem it crosses lands in
    layer 0 (Eterna's behaviour); do not "optimize" this to longest-first, or the
    layer assignment will diverge from the game.
    """
    n = len(pairs)
    # For each layer, the closing bases of placed pairs that could still be
    # crossed, kept in ascending order (mirrors Eterna's closingBasesPerDegree).
    closing_per_layer: List[List[int]] = []
    pairs_per_layer: List[List[BP]] = []

    for a in range(n):
        b = pairs[a]
        if b <= a:                      # skip unpaired (b == -1) and closing halves
            continue
        placed = False
        for d, closing in enumerate(closing_per_layer):
            # Drop closing bases we've already walked past.
            while closing and closing[0] < a:
                closing.pop(0)
            # Usable layer if nothing open here, or the nearest open pair closes
            # after this one (so no cross).
            if not closing or closing[0] > b:
                pairs_per_layer[d].append((a, b))
                if not closing or b != closing[0] - 1:
                    closing.insert(0, b)
                placed = True
                break
        if not placed:
            closing_per_layer.append([b])
            pairs_per_layer.append([(a, b)])

    layer_arrays: List[Pairs] = []
    for layer in pairs_per_layer:
        arr = [-1] * n
        for a, b in layer:
            arr[a] = b
            arr[b] = a
        layer_arrays.append(arr)
    return layer_arrays


# --------------------------------------------------------------------------- #
# Junction (multiloop / exterior-loop) detection within a nested layer
# --------------------------------------------------------------------------- #

@dataclass
class Junction:
    kind: str                 # "multiloop" or "exterior"
    layer: int                # which layer it was found in
    nway: int                 # number of stems meeting at the junction
    closing: BP | None        # closing base pair (None for exterior loop)
    members: Set[int]         # bordering residues (stem-boundary + unpaired linker)
    passed: bool = False
    crossed_by: List[int] = field(default_factory=list)


def _immediate_children(pairs: Pairs, lo: int, hi: int) -> Tuple[List[BP], List[int]]:
    """
    Immediate contents of the loop spanning open interval (lo, hi), exclusive.
    Returns (child_stem_outer_pairs, unpaired_positions) at this nesting level.
    Assumes `pairs` is nested (no crossings), so a forward scan is unambiguous.
    """
    children: List[BP] = []
    unpaired: List[int] = []
    k = lo + 1
    while k < hi:
        p = pairs[k]
        if p == -1:
            unpaired.append(k)
            k += 1
        elif p > k:                 # start of a child stem; jump past it
            children.append((k, p))
            k = p + 1
        else:                       # p < k: closing bracket of enclosing loop
            k += 1
    return children, unpaired


def find_junctions(layer_pairs: Pairs, layer_index: int) -> List[Junction]:
    """Find qualifying junctions in one nested layer."""
    n = len(layer_pairs)
    junctions: List[Junction] = []

    # Interior multiloops: any base pair whose loop contains >= 2 child stems.
    for i in range(n):
        j = layer_pairs[i]
        if j <= i:
            continue
        children, unpaired = _immediate_children(layer_pairs, i, j)
        if len(children) >= 2:                       # closing + >=2 branches => >=3-way
            members = {i, j}
            for a, b in children:
                members.add(a)
                members.add(b)
            members.update(unpaired)
            junctions.append(Junction(
                kind="multiloop", layer=layer_index,
                nway=len(children) + 1, closing=(i, j), members=members,
            ))

    # Exterior loop: stems that open at the top level (not enclosed by any pair).
    ext_children, ext_unpaired = _immediate_children(layer_pairs, -1, n)
    if len(ext_children) >= 3:                       # ">two stems"
        members = set(ext_unpaired)
        for a, b in ext_children:
            members.add(a)
            members.add(b)
        junctions.append(Junction(
            kind="exterior", layer=layer_index,
            nway=len(ext_children), closing=None, members=members,
        ))

    return junctions


# --------------------------------------------------------------------------- #
# Top-level check
# --------------------------------------------------------------------------- #

@dataclass
class Report:
    dbn: str
    pairs: Pairs
    layers: List[Pairs]
    junctions: List[Junction]
    crossed_res: Set[int]
    satisfied: bool

    def summary(self) -> str:
        lines = []
        n_pass = sum(j.passed for j in self.junctions)
        verdict = "PASS" if self.satisfied else "FAIL"
        lines.append(
            f"[{verdict}] {len(self.junctions)} qualifying junction(s), "
            f"{n_pass} crossed, {len(self.junctions) - n_pass} bare "
            f"(len {len(self.dbn)}, {len(self.crossed_res)} crossed residues)"
        )
        for k, jn in enumerate(self.junctions):
            tag = "crossed" if jn.passed else "BARE"
            close = f"pair {jn.closing}" if jn.closing else "exterior"
            why = (f"via residue(s) {jn.crossed_by[:6]}"
                   if jn.passed else "no bordering residue is in a crossed pair")
            lines.append(
                f"    J{k}: {jn.nway}-way {jn.kind} (layer {jn.layer}, {close}) "
                f"-> {tag}; {why}"
            )
        return "\n".join(lines)


def check(dbn: str) -> Report:
    """Parse a structure and evaluate the crossed-junction constraint."""
    pairs = parse_dotbracket(dbn)
    cross_res = crossed_residues(pairs)
    layers = separate_layers(pairs)

    junctions: List[Junction] = []
    for li, layer in enumerate(layers):
        junctions.extend(find_junctions(layer, li))

    for jn in junctions:
        hit = sorted(jn.members & cross_res)
        jn.crossed_by = hit
        jn.passed = len(hit) > 0

    satisfied = all(jn.passed for jn in junctions)
    return Report(dbn, pairs, layers, junctions, cross_res, satisfied)


# --------------------------------------------------------------------------- #
# Helpers for the example generator / validation
# --------------------------------------------------------------------------- #

def stems(pairs: Pairs) -> List[List[BP]]:
    """Group base pairs into contiguous helices (stems). Port of SecStruct.stems()."""
    result: List[List[BP]] = []
    for i in range(len(pairs)):
        j = pairs[i]
        if j <= i:
            continue
        attached = False
        for stem in result:
            for (a, b) in stem:
                if (a == i - 1 and b == j + 1) or (a == i + 1 and b == j - 1) \
                        or (b == i - 1 and a == j + 1) or (b == i + 1 and a == j - 1):
                    stem.append((i, j))
                    attached = True
                    break
            if attached:
                break
        if not attached:
            result.append([(i, j)])
    return result


def min_stem_length(dbn: str) -> int:
    """Length (in bp) of the shortest helix; 0 if unpaired."""
    st = stems(parse_dotbracket(dbn))
    return min((len(s) for s in st), default=0)


# --------------------------------------------------------------------------- #
# Sanity checks
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    cases = [
        # A 3-way junction, no pseudoknot -> should FAIL (bare junction).
        ("bare 3-way",
         "....((((....))))((((....))))((((....))))....",
         False),
        # Same, but a pseudoknot ([ at pos 1) reaches into the first hairpin loop,
        # crossing helix 1; its exterior end (pos 1) borders the junction -> PASS.
        ("crossed 3-way",
         ".[..((((.]..))))((((....))))((((....))))....",
         True),
        # Simple hairpin, no junction at all -> vacuously PASS.
        ("hairpin only", "((((....))))", True),
    ]
    for name, dbn, expected in cases:
        rep = check(dbn)
        ok = "OK" if rep.satisfied == expected else "!!! MISMATCH"
        print(f"[{ok}] {name}: expected {'PASS' if expected else 'FAIL'}")
        print(rep.summary())
        print()

    # Layering parity with Eterna's getParenthesis: a SHORT stem that opens before
    # a LONGER stem it crosses must land in layer 0 (left-to-right, not longest-first).
    dbn = "[[[[......(((((((((....]]]]......)))))))))"
    layers = separate_layers(parse_dotbracket(dbn))
    layer0_opener0 = layers[0][0] > 0    # is the pair opening at index 0 in layer 0?
    print(f"[{'OK' if layer0_opener0 else '!!! MISMATCH'}] "
          f"layer parity: short early stem is in layer 0 (matches Eterna getParenthesis)")
