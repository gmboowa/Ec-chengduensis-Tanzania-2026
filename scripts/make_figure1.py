#!/usr/bin/env python3
"""Figure 1 - recombination-filtered core-genome phylogeny of 43 E. chengduensis genomes.

Reproduces Figure 1 of the manuscript. Needs Python 3 with matplotlib and openpyxl.

    python scripts/make_figure1.py --out Figure1 --no-note --contrast high --dpi 1200

Inputs default to trees/Ec_43_phylo_iqtree.nwk, tables/Additional_file_1_Table_S1.xlsx
and tables/virulence_per_genome_43.tsv; override with --tree, --table, --virulence.

Branch lengths and topology are taken unchanged from the IQ-TREE output. The only
display choice is where the unrooted tree is drawn from (see --root).
"""
import argparse, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
import openpyxl
from newick import parse, tips, ladderize, layout, Node

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
TABLE = os.path.join(_ROOT, "tables", "Additional_file_1_Table_S1.xlsx")
VIR   = os.path.join(_ROOT, "tables", "virulence_per_genome_43.tsv")
NWK   = os.path.join(_ROOT, "trees", "Ec_43_phylo_iqtree.nwk")

# palette sampled from the figure currently in the manuscript
ORANGE, BLUE = "#A84400", "#004C80"
AMR_LOW, AMR_MED, AMR_HIGH = "#CCE2ED", "#77A9C6", "#123E60"
ROW_ALT, ROW_STUDY = "#F6F7F8", "#FBEEED"
INK, MUTED, RULE = "#1A1A1A", "#5A5A5A", "#B9BDC1"
LEADER, BRANCH_LW, CIRC_LW, FS_BASE = "#C9CDD1", 0.85, 1.15, 10.4

def high_contrast():
    """Deeper fills, pure-black ink and heavier strokes, for print and small sizes."""
    global AMR_LOW, AMR_MED, AMR_HIGH, ROW_ALT, ROW_STUDY
    global INK, RULE, LEADER, BRANCH_LW, CIRC_LW, FS_BASE
    AMR_LOW, AMR_MED, AMR_HIGH = "#BBD8EA", "#4F8CB5", "#07263F"
    ROW_ALT, ROW_STUDY = "#ECEFF2", "#FBE3E0"
    INK, RULE, LEADER = "#000000", "#7A8085", "#9AA0A6"
    BRANCH_LW, CIRC_LW, FS_BASE = 1.25, 1.55, 11.2
STUDY = ("A55822", "A55955")

# ----------------------------------------------------------------- data
def load(table=None, vir=None):
    wb = openpyxl.load_workbook(table or TABLE)
    ws = wb["Table S1"]
    h = [c.value for c in ws[1]]
    ix = {k: h.index(k) for k in
          ("Genome", "Country", "Isolation source", "AMR genes (n)",
           "Plasmid replicons", "ANI to type strain (%)", "Retained in final panel")}
    d = {}
    for r in ws.iter_rows(min_row=2, values_only=True):
        if r[ix["Retained in final panel"]] != "Yes":
            continue
        src = (r[ix["Isolation source"]] or "Not reported").strip()
        d[r[ix["Genome"]]] = dict(
            country=r[ix["Country"]],
            source=src[0].upper() + src[1:],
            amr=int(r[ix["AMR genes (n)"]]),
            plasmid=not str(r[ix["Plasmid replicons"]]).lower().startswith("no plasmid"),
            ani=float(r[ix["ANI to type strain (%)"]]),
        )
    for line in open(vir or VIR):
        p = line.rstrip("\n").split("\t")
        if p[0] in d:
            d[p[0]]["ybt"] = p[0] != "genome" and int(p[1]) > 0
    return d

# ----------------------------------------------------------------- rooting
def longest_branch(root):
    best = None
    def walk(n):
        nonlocal best
        if n.parent is not None and (best is None or n.length > best.length):
            best = n
        for c in n.children:
            walk(c)
    walk(root)
    return best

def reroot_midpoint_of(target):
    """Root the tree at the midpoint of the branch leading to `target`."""
    L = target.length
    p = target.parent
    R = Node()
    p.children.remove(target)
    target.parent = R
    target.length = L / 2.0
    R.children.append(target)
    cur, prev, carry = p, R, L / 2.0
    while cur is not None:
        nxt = cur.parent
        nxt_len = cur.length
        if nxt is not None:
            nxt.children.remove(cur)
        cur.parent = prev
        cur.length = carry
        prev.children.append(cur)
        prev, carry, cur = cur, nxt_len, nxt
    # drop any internal node left with a single child
    def collapse(n):
        for c in list(n.children):
            collapse(c)
        if n.parent is not None and len(n.children) == 1:
            c = n.children[0]
            c.length += n.length
            c.parent = n.parent
            n.parent.children[n.parent.children.index(n)] = c
    collapse(R)
    return R

def total_length(n):
    return n.length + sum(total_length(c) for c in n.children)

# ----------------------------------------------------------------- drawing
def draw(tree, data, out, scalebar, title_root, dpi=600):
    ladderize(tree, ascending=True)
    # put the clade holding the study isolates on top
    if len(tree.children) == 2:
        a = {t.name for t in tips(tree.children[0])}
        if not (set(STUDY) & a):
            tree.children.reverse()
    n = layout(tree)
    order = list(tips(tree))
    assert n == 43, n
    xmax = max(t.x for t in order)

    # geometry, in axes coordinates 0..1
    TREE0, TREE1 = 0.005, 0.245
    ANI_X, ANI_W = 0.262, 0.016
    GEN_X, CTY_X, SRC_X = 0.292, 0.452, 0.590
    AMR_X, PLA_X, YBT_X = 0.846, 0.912, 0.968
    TOP, BOT = 0.962, 0.168                      # rows live between these
    rowh = (TOP - BOT) / n
    sx = lambda v: TREE0 + (TREE1 - TREE0) * (v / xmax)
    sy = lambda i: TOP - (i + 0.5) * rowh

    fig = plt.figure(figsize=(14.0, 16.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    FS = FS_BASE

    # row bands
    for i, t in enumerate(order):
        y = TOP - (i + 1) * rowh
        if t.name in STUDY:
            ax.add_patch(Rectangle((ANI_X - 0.004, y), 1 - ANI_X - 0.004, rowh,
                                   fc=ROW_STUDY, ec="none", zorder=0))
        elif i % 2 == 0:
            ax.add_patch(Rectangle((ANI_X - 0.004, y), 1 - ANI_X - 0.004, rowh,
                                   fc=ROW_ALT, ec="none", zorder=0))

    # tree
    def branches(nd):
        y = sy(nd.y)
        if nd.parent is not None:
            ax.plot([sx(nd.parent.x), sx(nd.x)], [y, y], color=INK, lw=BRANCH_LW,
                    solid_capstyle="butt", zorder=3)
        if nd.children:
            ys = [sy(c.y) for c in nd.children]
            ax.plot([sx(nd.x), sx(nd.x)], [min(ys), max(ys)], color=INK, lw=BRANCH_LW,
                    solid_capstyle="butt", zorder=3)
            for c in nd.children:
                branches(c)
    branches(tree)

    for i, t in enumerate(order):
        y = sy(i)
        ax.plot([sx(t.x) + 0.004, ANI_X - 0.006], [y, y], color=LEADER,
                lw=0.75, ls=(0, (1, 2.0)), zorder=1)
        d = data[t.name]
        # ANI band
        ax.add_patch(Rectangle((ANI_X, y - rowh * 0.34), ANI_W, rowh * 0.68,
                               fc=ORANGE if d["ani"] < 97 else BLUE, ec="none", zorder=4))
        study = t.name in STUDY
        ax.text(GEN_X, y, t.name, va="center", ha="left", fontsize=FS,
                color="#B2182B" if study else INK,
                fontweight="bold" if study else "normal", zorder=5)
        ax.text(CTY_X, y, d["country"], va="center", ha="left", fontsize=FS,
                color=INK, zorder=5)
        ax.text(SRC_X, y, d["source"], va="center", ha="left", fontsize=FS,
                color=INK, zorder=5)
        a = d["amr"]
        fc = AMR_LOW if a <= 5 else AMR_MED if a <= 8 else AMR_HIGH
        ax.add_patch(FancyBboxPatch((AMR_X - 0.019, y - rowh * 0.30), 0.038, rowh * 0.60,
                                    boxstyle="round,pad=0,rounding_size=0.004",
                                    fc=fc, ec="none", zorder=4))
        ax.text(AMR_X, y, str(a), va="center", ha="center", fontsize=FS,
                color="white" if a > 8 else INK, zorder=5)
        for xc, on, col in ((PLA_X, d["plasmid"], BLUE), (YBT_X, d.get("ybt", False), ORANGE)):
            ax.add_patch(Circle((xc, y), rowh * 0.26, fc=col if on else "white",
                                ec=col, lw=CIRC_LW, zorder=5,
                                transform=ax.transData))

    # headers
    hy = TOP + 0.011
    for x, lab, ha in ((GEN_X, "Genome", "left"), (CTY_X, "Country", "left"),
                       (SRC_X, "Isolation source", "left")):
        ax.text(x, hy, lab, va="bottom", ha=ha, fontsize=FS + 0.8, fontweight="bold", color=INK)
    ax.text(ANI_X + ANI_W / 2, hy, "ANI", va="bottom", ha="center",
            fontsize=FS + 0.8, fontweight="bold", color=INK)
    ax.text(AMR_X, hy, "AMR\ngenes", va="bottom", ha="center", fontsize=FS + 0.8,
            fontweight="bold", color=INK, linespacing=1.15)
    ax.text(PLA_X, hy, "Plasmid\nreplicon", va="bottom", ha="center", fontsize=FS + 0.8,
            fontweight="bold", color=INK, linespacing=1.15)
    ax.text(YBT_X, hy, "ybtP/\nybtQ", va="bottom", ha="center", fontsize=FS + 0.8,
            fontweight="bold", color=INK, style="italic", linespacing=1.15)
    ax.plot([TREE0, 0.995], [TOP + 0.006, TOP + 0.006], color=RULE, lw=1.3)
    ax.text(TREE0, hy, "Phylogeny", va="bottom", ha="left",
            fontsize=FS + 0.8, fontweight="bold", color=INK)

    # ---- legend block
    ax.plot([TREE0, 0.995], [BOT - 0.018, BOT - 0.018], color=RULE, lw=1.3)
    L1, L2 = BOT - 0.052, BOT - 0.094
    bar = scalebar / xmax * (TREE1 - TREE0)
    ax.plot([TREE0, TREE0 + bar], [L1, L1], color=INK, lw=1.3)
    for xe in (TREE0, TREE0 + bar):
        ax.plot([xe, xe], [L1 - 0.005, L1 + 0.005], color=INK, lw=1.3)
    ax.text(TREE0, L1 - 0.015, f"{scalebar:g} substitutions per site",
            va="top", ha="left", fontsize=FS, color=INK)
    ax.text(TREE0, L2 - 0.008, "Red labels: this study", va="top", ha="left",
            fontsize=FS, color="#B2182B")

    LX = 0.300
    ax.text(LX, L1 + 0.015, "ANI to type strain GCF_001984825.2", va="bottom",
            ha="left", fontsize=FS, fontweight="bold", color=INK)
    for dx, col, lab in ((0.000, ORANGE, "95.67\u201395.80%  (n = 22)"),
                         (0.250, BLUE, "99.03\u201399.99%  (n = 21)")):
        ax.add_patch(Rectangle((LX + dx, L1 - 0.007), 0.016, 0.014, fc=col, ec="none"))
        ax.text(LX + dx + 0.023, L1, lab, va="center", ha="left", fontsize=FS, color=INK)

    ax.text(LX, L2 + 0.015, "Detected resistance genes", va="bottom", ha="left",
            fontsize=FS, fontweight="bold", color=INK)
    for dx, col, lab in ((0.000, AMR_LOW, "\u22645"), (0.062, AMR_MED, "6\u20138"),
                         (0.124, AMR_HIGH, "\u22659")):
        ax.add_patch(FancyBboxPatch((LX + dx, L2 - 0.009), 0.034, 0.018,
                                    boxstyle="round,pad=0,rounding_size=0.004",
                                    fc=col, ec="none"))
        ax.text(LX + dx + 0.017, L2, lab, va="center", ha="center", fontsize=FS - 0.4,
                color="white" if col == AMR_HIGH else INK)
    ax.add_patch(Circle((LX + 0.250, L2), 0.0060, fc=INK, ec=INK, lw=CIRC_LW))
    ax.text(LX + 0.268, L2, "Detected", va="center", ha="left", fontsize=FS, color=INK)
    ax.add_patch(Circle((LX + 0.380, L2), 0.0060, fc="white", ec=INK, lw=CIRC_LW))
    ax.text(LX + 0.398, L2, "Not detected", va="center", ha="left", fontsize=FS, color=INK)
    if title_root:
        ax.text(0.995, 0.016, title_root, va="bottom", ha="right",
                fontsize=FS - 1.0, color=MUTED)

    for ext, kw in (("png", dict(dpi=dpi)), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{out}.{ext}", facecolor="white", **kw)
    plt.close(fig)
    return order

# ----------------------------------------------------------------- main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", choices=("longest", "file"), default="longest")
    ap.add_argument("--scalebar", type=float, default=0.05)
    ap.add_argument("--out", default="Figure1")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--no-note", action="store_true")
    ap.add_argument("--contrast", choices=("normal", "high"), default="normal")
    ap.add_argument("--table", default=TABLE, help="Additional file 1: Table S1 (.xlsx)")
    ap.add_argument("--virulence", default=VIR, help="per-genome virulence table (.tsv)")
    ap.add_argument("--tree", default=NWK, help="Newick tree from IQ-TREE")
    a = ap.parse_args()
    if a.contrast == "high":
        high_contrast()
    data = load(a.table, a.virulence)
    tree = parse(open(a.tree).read())
    before = total_length(tree)
    if a.root == "longest":
        tgt = longest_branch(tree)
        tree = reroot_midpoint_of(tgt)
        note = ("Unrooted maximum-likelihood tree, drawn from the midpoint of its "
                "longest branch; branch lengths and topology unchanged.")
    else:
        note = ("Unrooted maximum-likelihood tree, drawn as written by IQ-TREE "
                "(basal trifurcation).")
    after = total_length(tree)
    assert abs(before - after) < 1e-12, (before, after)
    order = draw(tree, data, a.out, a.scalebar, "" if a.no_note else note, a.dpi)
    print("tips drawn     :", len(order))
    print("tree length    : %.6f before / %.6f after rooting" % (before, after))
    print("tip order (top): ", ", ".join(t.name for t in order[:4]))
    print("wrote          :", a.out + ".png/.pdf/.svg")
