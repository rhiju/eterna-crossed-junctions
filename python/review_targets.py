"""
Run the crossed-junction constraint over the real Eterna OpenKnotAIDesignData
targets and print a pass/fail table with reasons, as a real-world sanity review.

Data source (Git-LFS): eternagame/OpenKnotAIDesignData/Targets/*.csv
CSVs are cached under ../data/. Use --detail to print per-junction breakdowns.
"""

from __future__ import annotations

import argparse
import csv
import os
import urllib.request

import crossed_junctions as cj

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
BASE_URL = ("https://media.githubusercontent.com/media/eternagame/"
            "OpenKnotAIDesignData/main/Targets/")
FILES = ["Rounds1and2_targets.csv", "Round3_targets.csv", "Round4_targets.csv"]

# Column names differ across files.
ID_COLS = ["ID", "Title"]
STRUCT_COLS = ["Secstruct", "Dot-bracket", "Structure"]


def ensure_data(fname: str) -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, fname)
    if not os.path.exists(path) or os.path.getsize(path) < 200:
        url = BASE_URL + fname
        print(f"  downloading {fname} ...")
        urllib.request.urlretrieve(url, path)
    return path


def pick(row: dict, names) -> str | None:
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n].strip()
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--detail", action="store_true", help="print per-junction detail")
    args = ap.parse_args()

    print(f"{'target':40s} {'len':>4} {'#junc':>5} {'pass':>4} {'verdict':>7}")
    print("-" * 66)

    totals = {"n": 0, "pass": 0, "with_junc": 0, "err": 0}
    detail_lines = []

    for fname in FILES:
        path = ensure_data(fname)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                struct = pick(row, STRUCT_COLS)
                ident = pick(row, ID_COLS) or "?"
                if not struct:
                    continue
                totals["n"] += 1
                try:
                    rep = cj.check(struct)
                except Exception as e:  # malformed notation -> report, keep going
                    totals["err"] += 1
                    print(f"{ident[:40]:40s} {'--':>4} {'--':>5} {'--':>4} "
                          f"  ERROR: {e}")
                    continue

                njunc = len(rep.junctions)
                npass = sum(j.passed for j in rep.junctions)
                if njunc:
                    totals["with_junc"] += 1
                verdict = "PASS" if rep.satisfied else "FAIL"
                if rep.satisfied:
                    totals["pass"] += 1
                print(f"{ident[:40]:40s} {len(struct):>4} {njunc:>5} "
                      f"{npass:>4} {verdict:>7}")
                if args.detail and njunc:
                    detail_lines.append(f"\n{ident}:\n{struct}\n{rep.summary()}")

    print("-" * 66)
    print(f"{totals['n']} targets: {totals['with_junc']} have qualifying junction(s); "
          f"{totals['pass']} satisfy the constraint; {totals['err']} parse errors.")
    print("(Targets with no >=3-way junction trivially satisfy the constraint.)")

    if args.detail:
        print("\n=== per-junction detail ===")
        print("\n".join(detail_lines))


if __name__ == "__main__":
    main()
