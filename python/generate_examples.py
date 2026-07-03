"""
Generate example structures (length 100, every stem >= 4 bp) for the
crossed-junction constraint, and self-verify each with crossed_junctions.check().

  positive_1..5 : two 3-way multiloops, EACH crossed -- via a pseudoknot that
                  bridges the two junctions, or a separate pseudoknot per junction.
  negative_1..5 : two 3-way multiloops, but NOT all crossed -- at least one
                  qualifying junction is left bare.

Every structure is built from the same nested two-multiloop scaffold:

    tail - P1( lA1 P2(hp) lA2 P3( lB1 P4(hp) lB2 P5(hp) lB3 ) lA3 ) tail

  Junction A: closed by P1, branches P2 & P3.
  Junction B: closed by P3, branches P4 & P5.

Pseudoknot strands ([], {}, <>) are dropped onto exact residues. The constraint
only credits a junction when a residue *bordering* it (a stem-boundary base or an
unpaired linker base) is in a crossed pair -- so to satisfy a junction we land a
PK strand on one of its unpaired linkers (lA2 / lB2), never merely inside a
branch hairpin loop.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, List

import crossed_junctions as cj

LENGTH = 100
HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.normpath(os.path.join(HERE, "..", "examples"))


class Builder:
    """Assembles a dot-bracket string, remembering named regions for PK overlay."""

    def __init__(self):
        self.s: List[str] = []
        self.marks: dict = {}

    def add(self, text: str, mark: str | None = None) -> int:
        start = len(self.s)
        self.s.extend(text)
        if mark:
            self.marks[mark] = (start, len(self.s))
        return start

    def dots(self, n: int, mark: str | None = None):
        self.add("." * n, mark)

    def opens(self, n: int):
        self.add("(" * n)

    def closes(self, n: int):
        self.add(")" * n)

    def first_n(self, mark: str, n: int) -> str:
        """Create and return a sub-mark covering the first n residues of `mark`."""
        a, b = self.marks[mark]
        assert b - a >= n
        name = f"{mark}[:{n}]"
        self.marks[name] = (a, a + n)
        return name

    def place_pk(self, open_mark: str, close_mark: str, opench: str, closech: str):
        o0, o1 = self.marks[open_mark]
        c0, c1 = self.marks[close_mark]
        assert (o1 - o0) == (c1 - c0), "PK strands must match in length"
        for p in range(o0, o1):
            assert self.s[p] == ".", f"PK opener site {p} not unpaired"
            self.s[p] = opench
        for p in range(c0, c1):
            assert self.s[p] == ".", f"PK closer site {p} not unpaired"
            self.s[p] = closech

    def build(self) -> str:
        return "".join(self.s)


def scaffold(s1=8, s3=8, sp=7, loop=4, la2=4, lb2=4, tail5=2) -> Builder:
    """Two-multiloop nested scaffold; 3' tail auto-sized so len == LENGTH."""
    base = tail5 + 2 * s1 + 2 * s3 + 6 * sp + 3 * loop + la2 + lb2 + 4
    tail3 = LENGTH - base
    assert tail3 >= 0, f"parameters overflow LENGTH by {-tail3}"

    b = Builder()
    b.dots(tail5, "tail5")
    b.opens(s1)                                    # P1 open
    b.dots(1, "lA1")
    b.opens(sp); b.dots(loop, "P2loop"); b.closes(sp)      # P2 hairpin
    b.dots(la2, "lA2")                             # junction-A linker (PK site)
    b.opens(s3)                                    # P3 open
    b.dots(1, "lB1")
    b.opens(sp); b.dots(loop, "P4loop"); b.closes(sp)      # P4 hairpin
    b.dots(lb2, "lB2")                             # junction-B linker (PK site)
    b.opens(sp); b.dots(loop, "P5loop"); b.closes(sp)      # P5 hairpin
    b.dots(1, "lB3")
    b.closes(s3)                                   # P3 close
    b.dots(1, "lA3")
    b.closes(s1)                                   # P1 close
    b.dots(tail3)
    return b


# --------------------------------------------------------------------------- #
# The ten examples
# --------------------------------------------------------------------------- #

def _pos1():
    # One PK bridging the junctions: opener on A-linker, closer on B-linker.
    # The pair crosses P3, so BOTH junctions gain a crossed linker.
    b = scaffold()
    b.place_pk("lA2", "lB2", "[", "]")
    return b.build()


def _pos2():
    # A separate PK per junction; each lands its closer on that junction's linker
    # while its opener sits in a branch hairpin loop (crossing that branch).
    b = scaffold(s1=8, s3=6, sp=5, loop=8)
    b.place_pk(b.first_n("P2loop", 4), "lA2", "[", "]")   # crosses P2 -> junction A
    b.place_pk(b.first_n("P4loop", 4), "lB2", "{", "}")   # crosses P4 -> junction B
    return b.build()


def _pos3():
    # Bridge PK, wider loops / shorter branch stems.
    b = scaffold(sp=6, loop=6)
    b.place_pk("lA2", "lB2", "[", "]")
    return b.build()


def _pos4():
    # Bridge PK with a long P1 neck and short branch stems.
    b = scaffold(s1=12, s3=6, sp=6, loop=4)
    b.place_pk("lA2", "lB2", "[", "]")
    return b.build()


def _pos5():
    # Bridge PK using <> and a 5-bp helix on 5-nt linkers.
    b = scaffold(s1=6, s3=6, sp=8, loop=4, la2=5, lb2=5, tail5=1)
    b.place_pk("lA2", "lB2", "<", ">")
    return b.build()


def _neg1():
    # PK crosses junction B only (opener in P4 loop, closer on B-linker).
    # Junction A has no crossed border -> bare.
    b = scaffold(s1=8, s3=6, sp=5, loop=8)
    b.place_pk(b.first_n("P4loop", 4), "lB2", "[", "]")
    return b.build()


def _neg2():
    # No pseudoknots at all -> both junctions bare.
    return scaffold().build()


def _neg3():
    # PK crosses junction A only -> junction B bare.
    b = scaffold(s1=8, s3=6, sp=5, loop=8)
    b.place_pk(b.first_n("P2loop", 4), "lA2", "[", "]")
    return b.build()


def _neg4():
    # Kissing loops between B's own branches (P4 loop <-> P5 loop). This crosses
    # P4 and P5 -- both junction-B stems -- so B is satisfied, but the crossing is
    # confined inside P3 and never touches a junction-A stem or linker -> A bare.
    b = scaffold(s1=8, s3=6, sp=5, loop=8)
    b.place_pk(b.first_n("P4loop", 4), b.first_n("P5loop", 4), "[", "]")
    return b.build()


def _neg5():
    # A PK against the P1 neck (5' tail <-> P2 loop) crosses P1 and P2 -- both
    # junction-A stems -- satisfying A, while junction B (P3/P4/P5) is untouched
    # -> B bare.
    b = scaffold(s1=8, s3=6, sp=5, loop=8, tail5=4)
    b.place_pk(b.first_n("tail5", 4), b.first_n("P2loop", 4), "[", "]")
    return b.build()


@dataclass
class Example:
    name: str
    fn: Callable[[], str]
    expect_pass: bool


EXAMPLES = [
    Example("positive_1_bridge", _pos1, True),
    Example("positive_2_inside_each", _pos2, True),
    Example("positive_3_wide_loops", _pos3, True),
    Example("positive_4_long_neck", _pos4, True),
    Example("positive_5_five_bp_pk", _pos5, True),
    Example("negative_1_only_B_crossed", _neg1, False),
    Example("negative_2_no_pk", _neg2, False),
    Example("negative_3_only_A_crossed", _neg3, False),
    Example("negative_4_B_kissing_A_bare", _neg4, False),
    Example("negative_5_A_neck_B_bare", _neg5, False),
]


def main():
    os.makedirs(EXAMPLES_DIR, exist_ok=True)
    print(f"{'name':28s} {'len':>4} {'minstem':>7} {'#junc':>5} "
          f"{'verdict':>7} {'expect':>6}  ok")
    all_ok = True
    for e in EXAMPLES:
        dbn = e.fn()
        rep = cj.check(dbn)
        ms = cj.min_stem_length(dbn)
        ok = (rep.satisfied == e.expect_pass and len(dbn) == LENGTH and ms >= 4)
        all_ok = all_ok and ok
        print(f"{e.name:28s} {len(dbn):>4} {ms:>7} {len(rep.junctions):>5} "
              f"{'PASS' if rep.satisfied else 'FAIL':>7} "
              f"{'PASS' if e.expect_pass else 'FAIL':>6}  {'OK' if ok else 'XX'}")
        with open(os.path.join(EXAMPLES_DIR, e.name + ".dbn"), "w") as fh:
            fh.write(dbn + "\n")

    print("\n" + ("ALL 10 EXAMPLES OK" if all_ok else "*** SOME EXAMPLES FAILED ***"))
    if not all_ok:
        raise SystemExit(1)

    # Detailed view, for the record.
    print("\n--- detail ---")
    for e in EXAMPLES:
        dbn = e.fn()
        print(f"\n{e.name}:")
        print(dbn)
        print(cj.check(dbn).summary())


if __name__ == "__main__":
    main()
