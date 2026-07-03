"""
Draw structures with draw_rna, highlighting junctions by crossed/bare status.

The nested "backbone" (layer 0 of the crossed-junction layering) is laid out with
draw_rna's renderer. Junction residues are colored GREEN if their junction is
crossed and RED if it is bare. The remaining (non-nested) pseudoknot pairs are
drawn as arcs ("strings") connecting the paired residues.

Usage:
    python3 draw_crossed.py                 # a default set of real targets
    python3 draw_crossed.py examples        # the 10 generated examples
    python3 draw_crossed.py "Kissing multiloops" SV_c   # named targets
"""

from __future__ import annotations

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Helvetica"
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path

# draw_rna renderer + layout constants. The repo nests a 'draw_rna' package inside
# the 'draw_rna' repo dir and its modules use absolute `import draw_rna.svg`. So we
# point the path at the repo dir and clear the cached (outer) package, making
# `draw_rna` resolve to the inner package that actually holds render_rna/svg.
import draw_rna as _pkg  # noqa: E402
sys.path.insert(0, os.path.dirname(_pkg.__file__))
for _m in [m for m in sys.modules if m == "draw_rna" or m.startswith("draw_rna.")]:
    del sys.modules[_m]
from draw_rna import render_rna  # noqa: E402
NODE_R, PRIMARY_SPACE, PAIR_SPACE = 10, 20, 20

import crossed_junctions as cj
import review_targets as rt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
EX_DIR = os.path.normpath(os.path.join(HERE, "..", "examples"))

GREEN = "#3aa14f"   # crossed junction
RED = "#d1483a"     # bare junction
GRAY = "#cccccc"    # everything else
CHAIN = "#999999"
BB_PAIR = "#777777"
PK_COLORS = ["#8a2be2", "#1e90ff", "#e77300", "#00a2a2"]  # per PK layer


def load_targets() -> dict:
    """Map target id/title -> dot-bracket, from the cached OpenKnot CSVs."""
    out = {}
    for fname in rt.FILES:
        path = rt.ensure_data(fname)
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                struct = rt.pick(row, rt.STRUCT_COLS)
                ident = rt.pick(row, rt.ID_COLS)
                if struct and ident:
                    out[ident] = struct
    return out


def backbone_dbn(layer0: list) -> str:
    """Pure '()' dot-bracket for a nested layer."""
    chars = []
    for i, p in enumerate(layer0):
        chars.append("(" if p > i else (")" if 0 <= p < i else "."))
    return "".join(chars)


def coords_for(dbn_backbone: str):
    r = render_rna.RNARenderer()
    r.setup_tree(dbn_backbone, NODE_R, PRIMARY_SPACE, PAIR_SPACE, 1, 0)
    # draw_rna's y grows downward; negate so figures read upright.
    xs = list(r.xarray_)
    ys = [-y for y in r.yarray_]
    return xs, ys


def arc(ax, x0, y0, x1, y1, color, lw=1.6):
    """Quadratic-bezier arc bulging perpendicular to the chord."""
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    dx, dy = x1 - x0, y1 - y0
    dist = (dx * dx + dy * dy) ** 0.5 or 1.0
    # perpendicular offset, scaled by chord length
    ox, oy = -dy / dist, dx / dist
    bulge = 0.28 * dist
    cx, cy = mx + ox * bulge, my + oy * bulge
    path = Path([(x0, y0), (cx, cy), (x1, y1)],
                [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(path, fill=False, edgecolor=color, lw=lw, alpha=0.85))


def draw_one(name: str, struct: str, outpath: str):
    rep = cj.check(struct)
    layers = rep.layers
    xs, ys = coords_for(backbone_dbn(layers[0]))
    n = len(struct)

    # residue -> color by junction status (bare wins over crossed for shared bases)
    crossed_res, bare_res = set(), set()
    for jn in rep.junctions:
        (crossed_res if jn.passed else bare_res).update(jn.members)
    node_color = {}
    for i in range(n):
        if i in bare_res:
            node_color[i] = RED
        elif i in crossed_res:
            node_color[i] = GREEN

    fig, ax = plt.subplots(figsize=(9, 9))

    # backbone chain
    ax.plot(xs, ys, "-", color=CHAIN, lw=1.0, zorder=1)
    # backbone pairs (layer 0)
    for i, p in enumerate(layers[0]):
        if p > i:
            ax.plot([xs[i], xs[p]], [ys[i], ys[p]], "-", color=BB_PAIR, lw=1.0, zorder=1)
    # pseudoknot pairs (layers 1+) as arcs
    for li in range(1, len(layers)):
        color = PK_COLORS[(li - 1) % len(PK_COLORS)]
        for i, p in enumerate(layers[li]):
            if p > i:
                arc(ax, xs[i], ys[i], xs[p], ys[p], color)

    # nodes
    for i in range(n):
        c = node_color.get(i, GRAY)
        rad = NODE_R * 0.62 if i in node_color else NODE_R * 0.42
        ax.add_patch(Circle((xs[i], ys[i]), rad, facecolor=c,
                            edgecolor="#444444" if i in node_color else "none",
                            lw=0.6, zorder=3))
    # 5'/3' labels
    ax.annotate("5'", (xs[0], ys[0]), fontsize=13, ha="right", va="top", weight="bold")
    ax.annotate("3'", (xs[-1], ys[-1]), fontsize=13, ha="left", va="bottom", weight="bold")

    n_pass = sum(j.passed for j in rep.junctions)
    verdict = "PASS" if rep.satisfied else "FAIL"
    ax.set_title(f"{name}   [{verdict}]   {n_pass}/{len(rep.junctions)} junctions crossed",
                 fontsize=20, weight="bold", color=GREEN if rep.satisfied else RED)

    # legend
    handles = [
        plt.Line2D([], [], marker="o", ls="", mfc=GREEN, mec="#444", label="crossed junction"),
        plt.Line2D([], [], marker="o", ls="", mfc=RED, mec="#444", label="bare junction"),
        plt.Line2D([], [], color=PK_COLORS[0], lw=1.8, label="pseudoknot pair"),
    ]
    ax.legend(handles=handles, fontsize=12, loc="upper left", frameon=True)

    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  {verdict}  {n_pass}/{len(rep.junctions)}  ->  {os.path.relpath(outpath, HERE)}")


DEFAULT_TARGETS = [
    "Kissing multiloops", "W04", "AK_PK100-3",   # pass
    "SV_c", "AK_PK240-3", "Guide RNA",           # fail
]


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    args = sys.argv[1:]

    if args == ["examples"]:
        items = []
        for f in sorted(os.listdir(EX_DIR)):
            if f.endswith(".dbn"):
                with open(os.path.join(EX_DIR, f)) as fh:
                    items.append((f[:-4], fh.read().strip()))
    else:
        targets = load_targets()
        names = args or DEFAULT_TARGETS
        items = []
        for nm in names:
            if nm not in targets:
                print(f"  (skip) target not found: {nm}")
                continue
            items.append((nm, targets[nm]))

    print(f"Drawing {len(items)} structure(s) -> {os.path.relpath(FIG_DIR, HERE)}/")
    for name, struct in items:
        safe = name.replace(" ", "_").replace("/", "_")
        draw_one(name, struct, os.path.join(FIG_DIR, safe + ".png"))


if __name__ == "__main__":
    main()
