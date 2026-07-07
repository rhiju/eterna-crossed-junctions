"""
Render a CSV of Eterna submissions to a multipage PDF: a 4x5 grid of secondary
structures per page, with multiloop/exterior junctions colored GREEN (crossed) or
RED (bare), pseudoknot pairs drawn as arcs, and each panel titled by its 'title'.

Usage:
    python3 draw_submissions_pdf.py [input.csv] [output.pdf]

Defaults to the CSV in ../submissions/ (gitignored player data). The structure
column is 'Structure' and the title column is 'title'.
"""

from __future__ import annotations

import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Helvetica"
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import draw_crossed as dc

HERE = os.path.dirname(os.path.abspath(__file__))
SUB_DIR = os.path.normpath(os.path.join(HERE, "..", "submissions"))

COLS, ROWS = 4, 5              # panels per page
PER_PAGE = COLS * ROWS
STRUCT_COL = "Structure"
TITLE_COL = "title"


def load_rows(path: str):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            struct = (row.get(STRUCT_COL) or "").strip()
            title = (row.get(TITLE_COL) or row.get("nid") or "").strip()
            if struct:
                rows.append((title, struct))
    return rows


def main():
    args = sys.argv[1:]
    if args and args[0].endswith(".csv"):
        in_csv = args[0]
    else:
        found = sorted(glob.glob(os.path.join(SUB_DIR, "*.csv")))
        if not found:
            sys.exit(f"No CSV given and none found in {SUB_DIR}")
        in_csv = found[0]
    out_pdf = (args[1] if len(args) > 1
               else os.path.join(SUB_DIR, os.path.splitext(os.path.basename(in_csv))[0] + ".pdf"))

    rows = load_rows(in_csv)
    npages = (len(rows) + PER_PAGE - 1) // PER_PAGE
    print(f"{len(rows)} structures -> {npages} pages -> {out_pdf}")

    n_pass = n_fail = n_err = 0
    with PdfPages(out_pdf) as pdf:
        for page in range(npages):
            fig, axes = plt.subplots(ROWS, COLS, figsize=(16, 20))
            axes = axes.ravel()
            batch = rows[page * PER_PAGE:(page + 1) * PER_PAGE]
            for k, ax in enumerate(axes):
                if k >= len(batch):
                    ax.axis("off")
                    continue
                title, struct = batch[k]
                label = title if len(title) <= 34 else title[:31] + "..."
                try:
                    rep = dc.render_on_ax(ax, label, struct, title_fontsize=10, ends=False)
                    if rep.satisfied:
                        n_pass += 1
                    else:
                        n_fail += 1
                except Exception as e:  # malformed notation etc.
                    n_err += 1
                    ax.axis("off")
                    ax.text(0.5, 0.5, f"{label}\n(could not render: {e})",
                            ha="center", va="center", fontsize=8, color="#aa0000",
                            transform=ax.transAxes, wrap=True)

            fig.legend(handles=dc.legend_handles(), fontsize=12, ncol=3,
                       loc="lower center", frameon=False,
                       bbox_to_anchor=(0.5, 0.005))
            fig.suptitle(
                f"Eterna submissions — crossed-junction check   "
                f"(page {page + 1}/{npages})",
                fontsize=16, weight="bold", y=0.997)
            fig.tight_layout(rect=(0, 0.02, 1, 0.985))
            pdf.savefig(fig, dpi=110)
            plt.close(fig)
            print(f"  page {page + 1}/{npages} done")

    print(f"\n{n_pass} satisfy, {n_fail} fail, {n_err} render errors.")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()
