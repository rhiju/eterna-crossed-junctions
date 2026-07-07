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

import crossed_junctions as cj
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


def summary_page(pdf, rows, src_name):
    """Cover page: overall counts, a color key, and the failing submissions."""
    results = []   # (title, satisfied, n_pass, n_total, error)
    for title, struct in rows:
        try:
            rep = cj.check(struct)
            results.append((title, rep.satisfied,
                            sum(j.passed for j in rep.junctions),
                            len(rep.junctions), None))
        except Exception as e:
            results.append((title, False, 0, 0, str(e)))

    n = len(results)
    n_pass = sum(1 for r in results if r[4] is None and r[1])
    n_fail = sum(1 for r in results if r[4] is None and not r[1])
    n_err = sum(1 for r in results if r[4] is not None)
    fails = [r for r in results if r[4] is None and not r[1]]

    fig = plt.figure(figsize=(16, 20))
    ax = fig.add_axes([0.06, 0.04, 0.88, 0.92])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def T(x, y, s, size=13, weight="normal", color="#222222", family=None):
        ax.text(x, y, s, transform=ax.transAxes, fontsize=size, weight=weight,
                color=color, va="top", ha="left", family=family)

    y = 0.985
    T(0.0, y, "Eterna submissions — crossed-junction check", size=26, weight="bold")
    y -= 0.035
    T(0.0, y, f"source: {src_name}", size=11, color="#666666")
    y -= 0.045
    T(0.0, y, f"{n} submissions:  ", size=16, weight="bold")
    ax.text(0.24, y, f"{n_pass} satisfy", transform=ax.transAxes, fontsize=16,
            weight="bold", color=dc.GREEN, va="top")
    ax.text(0.42, y, f"{n_fail} fail", transform=ax.transAxes, fontsize=16,
            weight="bold", color=dc.RED, va="top")
    ax.text(0.55, y, f"{n_err} unparsed", transform=ax.transAxes, fontsize=16,
            color="#888888", va="top")

    # ---- color key -------------------------------------------------------
    y -= 0.06
    T(0.0, y, "What the drawing shows", size=18, weight="bold")
    y -= 0.03
    T(0.0, y, "Each panel lays out the nested “backbone” of one submission; "
             "the non-nested pseudoknot pairs are drawn as arcs over it.",
      size=12, color="#444444")

    def marker(yy, kind, color, label, sub):
        if kind == "circle":
            ax.scatter([0.03], [yy], s=170, c=[color], edgecolors="#444444",
                       linewidths=0.8, transform=ax.transAxes, zorder=3)
        else:  # line
            ax.plot([0.015, 0.05], [yy, yy], color=color, lw=3,
                    transform=ax.transAxes)
        ax.text(0.08, yy, label, transform=ax.transAxes, fontsize=13,
                weight="bold", va="center", ha="left")
        ax.text(0.08, yy - 0.016, sub, transform=ax.transAxes, fontsize=11,
                color="#555555", va="center", ha="left")

    y -= 0.05
    marker(y, "circle", dc.GREEN, "green nucleotide  =  crossed junction",
           "a base bordering a ≥3-way junction (multiloop or exterior loop) that IS touched by a pseudoknot")
    y -= 0.055
    marker(y, "circle", dc.RED, "red nucleotide  =  bare junction",
           "a base bordering a ≥3-way junction that is NOT crossed by any pseudoknot — the case we flag")
    y -= 0.055
    marker(y, "circle", dc.GRAY, "gray nucleotide  =  everything else",
           "a base not bordering a ≥3-way junction (helices, hairpin/internal loops, tails)")
    y -= 0.055
    marker(y, "line", dc.BB_PAIR, "gray lines  =  backbone",
           "the chain plus the nested base pairs (layer 0) — the secondary structure that is drawn")
    y -= 0.055
    marker(y, "line", dc.PK_COLORS[0], "colored arcs  =  pseudoknot pairs",
           "the non-nested (crossing) pairs shown as strings between partners; one color per extra layer "
           "(purple, blue, orange, teal)")

    y -= 0.06
    T(0.0, y, "Rule:", size=14, weight="bold")
    ax.text(0.06, y, "a submission PASSES when every ≥3-way junction has at least one green "
                     "(crossed) base; it FAILS if any such junction is entirely bare (red).",
            transform=ax.transAxes, fontsize=12, color="#444444", va="top")

    # ---- failing list ----------------------------------------------------
    y -= 0.055
    T(0.0, y, f"Failing submissions ({n_fail})", size=18, weight="bold", color=dc.RED)
    y -= 0.03
    if n_err:
        T(0.0, y, f"(plus {n_err} that could not be parsed)", size=11, color="#888888")
        y -= 0.022
    col_x = [0.0, 0.5]
    start_y = y
    per_col = 14
    for i, (title, _sat, npass, ntot, _e) in enumerate(fails):
        cx = col_x[i // per_col] if i // per_col < len(col_x) else col_x[-1]
        ry = start_y - 0.026 * (i % per_col)
        label = title if len(title) <= 40 else title[:37] + "..."
        ax.text(cx, ry, f"•  {label}", transform=ax.transAxes, fontsize=12,
                va="top", ha="left")
        ax.text(cx + 0.34, ry, f"{npass}/{ntot} crossed", transform=ax.transAxes,
                fontsize=11, color=dc.RED, va="top", ha="left")

    pdf.savefig(fig, dpi=110)
    plt.close(fig)
    return n_pass, n_fail, n_err


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
    print(f"{len(rows)} structures -> {npages} grid pages + summary -> {out_pdf}")

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
                    dc.render_on_ax(ax, label, struct, title_fontsize=10, ends=False)
                except Exception as e:  # malformed notation etc.
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

        n_pass, n_fail, n_err = summary_page(pdf, rows, os.path.basename(in_csv))
        print("  summary page done (at end)")

    print(f"\n{n_pass} satisfy, {n_fail} fail, {n_err} render errors.")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()
