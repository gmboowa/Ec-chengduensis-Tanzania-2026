#!/usr/bin/env python3
"""Annotated phylogram of the two ANI clusters, in the style of the iTOL figure.

The branch separating the two clusters is about a thousand times longer than
any branch within either of them. Drawn on one axis, every within-cluster
relationship collapses to a vertical line and the figure shows nothing. Each
cluster is therefore drawn in its own panel at its own scale, with its own
scale bar, and the separation between them is stated rather than drawn.

Pure matplotlib with a hand-rolled Newick parser, so the figure regenerates
from the deposited files with no phylogenetics package installed.

  python3 draw_tree.py <tree.nwk> <Table_S1...xlsx> <virulence_per_genome.tsv> <out.png>

The two ANI bands are separated by a branch roughly a thousand times longer
than any within-band branch. Drawn to scale, every within-band relationship
collapses to a point, so that one branch is compressed and marked with a break;
the compression is stated in the caption and applies to no other branch.
"""
import csv
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


# ----------------------------------------------------------------- newick
class Node:
    __slots__ = ("name", "length", "support", "children", "x", "y")

    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.support = None
        self.children = []
        self.x = self.y = 0.0

    @property
    def is_leaf(self):
        return not self.children


def parse_newick(text):
    text = text.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        m = re.match(r"[^,:()]*", text[pos:])
        label = m.group(0)
        pos += len(label)
        if node.is_leaf:
            node.name = label
        elif label:
            try:
                node.support = float(label)
            except ValueError:
                node.name = label
        if pos < len(text) and text[pos] == ":":
            pos += 1
            m = re.match(r"-?[\d.eE+-]+", text[pos:])
            node.length = float(m.group(0))
            pos += len(m.group(0))
        return node

    return parse_node()


def leaves(node):
    if node.is_leaf:
        yield node
    for c in node.children:
        yield from leaves(c)


def nodes(node):
    yield node
    for c in node.children:
        yield from nodes(c)


# ----------------------------------------------------------------- metadata
def load_metadata(xlsx, virulence_tsv):
    from openpyxl import load_workbook
    ws = load_workbook(xlsx)["Table S1"]
    hdr = [c.value for c in ws[1]]
    meta = {}
    for i in range(2, ws.max_row + 1):
        d = dict(zip(hdr, [ws.cell(i, j).value for j in range(1, ws.max_column + 1)]))
        meta[d["Genome"]] = d
    vir = {}
    with open(virulence_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vir[r["genome"]] = r["virulence_genes"]
    return meta, vir


# ----------------------------------------------------------------- palettes
# Okabe-Ito, extended; chosen for colour-vision deficiency and for printing.
OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
         "#56B4E9", "#F0E442", "#999999", "#882255", "#44AA99",
         "#117733", "#DDCC77", "#AA4499", "#6699CC", "#661100"]
GREYS = {"present": "#333333", "absent": "#FFFFFF"}


def amr_band(n):
    if n is None:
        return "unknown"
    return "low (<=5)" if n <= 5 else ("medium (6-8)" if n <= 8 else "high (>=9)")



def midpoint_root(root):
    """Re-root at the midpoint of the longest tip-to-tip path.

    IQ-TREE writes an unrooted tree as a trifurcation. Drawn as-is, whichever
    taxa happen to sit at that trifurcation appear to hang off the left edge
    unattached. Midpoint rooting is the conventional display choice and is
    what iTOL applies by default; it changes no branch length and asserts no
    biological root.
    """
    adj = {}
    def walk(n):
        for c in n.children:
            adj.setdefault(n, []).append((c, c.length))
            adj.setdefault(c, []).append((n, c.length))
            walk(c)
    walk(root)
    tips = list(leaves(root))

    def far(src):
        seen, stack, best = {src: 0.0}, [src], (src, 0.0)
        while stack:
            cur = stack.pop()
            for nb, w in adj.get(cur, []):
                if nb not in seen:
                    seen[nb] = seen[cur] + w
                    stack.append(nb)
                    if nb.is_leaf and seen[nb] > best[1]:
                        best = (nb, seen[nb])
        return best, seen

    (a, _), _ = far(tips[0])
    (b, total), dist_a = far(a)

    # walk back from b towards a until half the path is covered
    prev, path = {}, []
    seen, stack = {a}, [a]
    while stack:
        cur = stack.pop()
        for nb, w in adj.get(cur, []):
            if nb not in seen:
                seen.add(nb); prev[nb] = (cur, w); stack.append(nb)
    node = b
    while node is not a:
        par, w = prev[node]
        path.append((node, par, w))
        node = par

    acc = 0.0
    for child, parent, w in path:
        if acc + w >= total / 2:
            offset = total / 2 - acc
            break
        acc += w
    else:
        return root

    new = Node()
    def rebuild(n, came_from, stem):
        out = Node()
        out.name, out.support, out.length = n.name, n.support, stem
        for nb, w in adj.get(n, []):
            if nb is not came_from:
                out.children.append(rebuild(nb, n, w))
        return out

    new.children = [rebuild(child, parent, offset),
                    rebuild(parent, child, w - offset)]
    return new


def prune(node, keep):
    """Copy the tree, dropping tips not in `keep` and collapsing dead ends."""
    if node.is_leaf:
        return node if node.name in keep else None
    kids = [prune(c, keep) for c in node.children]
    kids = [k for k in kids if k is not None]
    if not kids:
        return None
    new = Node()
    new.name, new.length, new.support = node.name, node.length, node.support
    if len(kids) == 1:
        only = kids[0]
        only.length += node.length
        return only
    new.children = kids
    return new


def layout(root):
    tips = list(leaves(root))
    for i, t in enumerate(tips):
        t.y = i
    def place(n, x0):
        n.x = x0 + n.length
        for c in n.children:
            place(c, n.x)
    place(root, 0.0)
    for n in sorted(nodes(root), key=lambda z: -z.x):
        if not n.is_leaf:
            n.y = sum(c.y for c in n.children) / len(n.children)
    return tips


def nice_bar(xmax):
    import math
    raw = xmax / 4.0
    if raw <= 0:
        return xmax
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def draw_panel(ax, root, meta, vir, cmap, smap, amap, title, n_total):
    tips = layout(root)
    xmax = max(n.x for n in nodes(root)) or 1.0

    def draw(node, parent_x):
        ax.plot([parent_x, node.x], [node.y, node.y], color="#2b2b2b", lw=1.0,
                solid_capstyle="butt", zorder=2)
        if node.children:
            ys = [c.y for c in node.children]
            ax.plot([node.x, node.x], [min(ys), max(ys)], color="#2b2b2b",
                    lw=1.0, solid_capstyle="butt", zorder=2)
            for c in node.children:
                draw(c, node.x)
    draw(root, 0.0)

    for n in nodes(root):
        if n.children and n.support is not None and n.support >= 90:
            ax.plot(n.x, n.y, "o", ms=3.2, mfc="#2b2b2b", mec="none", zorder=4)

    label_x = xmax * 1.04
    col_w, gap = xmax * 0.042, xmax * 0.016
    track_x = label_x + xmax * 0.42

    for t in tips:
        m = meta.get(t.name, {})
        study = m.get("Origin") == "This study"
        ax.text(label_x, t.y, t.name, va="center", ha="left", fontsize=8.0,
                fontweight="bold" if study else "normal",
                color="#B2182B" if study else "#222222")
        x = track_x
        vals = [cmap.get(m.get("Country") or "n/a", "#E8E8E8"),
                smap.get((m.get("Isolation source") or "not provided").lower(), "#E8E8E8"),
                amap[amr_band(m.get("AMR genes (n)"))]]
        for colour in vals:
            ax.add_patch(Rectangle((x, t.y - 0.40), col_w, 0.80, color=colour,
                                   lw=0.3, ec="white", zorder=3))
            x += col_w + gap
        rep = (m.get("Plasmid replicons") or "")
        for present in (bool(rep) and "No plasmid" not in rep,
                        "ybt" in (vir.get(t.name, "") or "")):
            ax.add_patch(Rectangle((x, t.y - 0.40), col_w, 0.80,
                                   color=GREYS["present"] if present else GREYS["absent"],
                                   lw=0.5, ec="#888888", zorder=3))
            x += col_w + gap

    for i, h in enumerate(["Country", "Source", "AMR load", "Plasmid", "ybtPQ"]):
        ax.text(track_x + i * (col_w + gap) + col_w / 2, len(tips) - 0.35, h,
                rotation=45, ha="left", va="bottom", fontsize=7.6, color="#555555")

    bar = nice_bar(xmax)
    y0 = -1.9
    ax.plot([0, bar], [y0, y0], color="#2b2b2b", lw=1.4)
    for xx in (0, bar):
        ax.plot([xx, xx], [y0 - 0.22, y0 + 0.22], color="#2b2b2b", lw=1.4)
    ax.text(bar / 2, y0 - 0.65, f"{bar:g} substitutions per site", ha="center",
            va="top", fontsize=8, color="#333333")

    ax.set_title(title, loc="left", fontsize=10.5, color="#222222", pad=10)
    ax.set_xlim(-xmax * 0.04, track_x + 5 * (col_w + gap) + xmax * 0.03)
    ax.set_ylim(y0 - 1.6, len(tips) + 1.2)
    ax.axis("off")
    return tips


def main(argv):
    nwk, xlsx, vir_tsv, out = argv[1:5]
    root = parse_newick(open(nwk).read())
    meta, vir = load_metadata(xlsx, vir_tsv)

    # the longest branch separates the two clusters
    deep = max((n for n in nodes(root) if n is not root), key=lambda n: n.length)
    upper_names = {t.name for t in leaves(deep)}
    all_names = {t.name for t in leaves(root)}
    lower_names = all_names - upper_names

    # The reference is the same genome as SRR5247447 (0 SNPs between them);
    # showing both would duplicate a tip.
    DROP = {"Reference"}
    upper_names -= DROP
    lower_names -= DROP

    upper = midpoint_root(prune(root, upper_names))
    lower = midpoint_root(prune(root, lower_names))
    upper.length = lower.length = 0.0

    tips_all = [t for t in leaves(root) if t.name != "Reference"]
    countries = sorted({meta.get(t.name, {}).get("Country") or "n/a" for t in tips_all})
    cmap = {c: OKABE[i % len(OKABE)] for i, c in enumerate(countries)}
    cmap["n/a"] = "#E8E8E8"
    sources = sorted({(meta.get(t.name, {}).get("Isolation source") or "not provided").lower()
                      for t in tips_all})
    smap = {s: OKABE[(i + 4) % len(OKABE)] for i, s in enumerate(sources)}
    smap["not provided"] = "#E8E8E8"
    amrs = ["low (<=5)", "medium (6-8)", "high (>=9)", "unknown"]
    amap = dict(zip(amrs, ["#BBDDEE", "#6699CC", "#114477", "#E8E8E8"]))

    n_lo, n_hi = len(lower_names), len(upper_names)
    h = max(6.5, 0.30 * max(n_lo, n_hi) + 2.6)
    fig, axes = plt.subplots(1, 2, figsize=(17.5, h), dpi=200)

    draw_panel(axes[0], lower, meta, vir, cmap, smap, amap,
               f"A   Lower ANI cluster  (95.67-95.80% to the type strain, n={n_lo})",
               len(tips_all))
    draw_panel(axes[1], upper, meta, vir, cmap, smap, amap,
               f"B   Upper ANI cluster  (99.03-99.99% to the type strain, n={n_hi})",
               len(tips_all))

    handles = [Line2D([], [], marker="s", ls="", ms=7, mfc=cmap[c], mec="none", label=c)
               for c in countries]
    l1 = fig.legend(handles=handles, title="Country", loc="upper left",
                    bbox_to_anchor=(0.995, 0.97), frameon=False, fontsize=7.6,
                    title_fontsize=8.8, labelspacing=0.32)
    fig.add_artist(l1)
    h2 = [Line2D([], [], marker="s", ls="", ms=7, mfc=amap[a], mec="none", label=a)
          for a in amrs[:3]]
    h2 += [Line2D([], [], marker="o", ls="", ms=5, mfc="#2b2b2b", mec="none",
                  label="bootstrap >= 90"),
           Line2D([], [], marker="s", ls="", ms=7, mfc="#333333", mec="none",
                  label="feature present"),
           Line2D([], [], marker="s", ls="", ms=7, mfc="white", mec="#888888",
                  label="feature absent")]
    fig.legend(handles=h2, title="Annotation", loc="upper left",
               bbox_to_anchor=(0.995, 0.46), frameon=False, fontsize=7.6,
               title_fontsize=8.8, labelspacing=0.32)

    fig.suptitle(
        f"Recombination-filtered core-genome phylogeny of {len(tips_all)} "
        f"Enterobacter chengduensis genomes.  The two clusters are separated by "
        f"122,329-123,080 core-genome SNPs and are drawn in separate panels, each "
        f"to its own scale -- the scale bars differ -- and each is midpoint-rooted. "
        f"Isolates from this study are in red.",
        fontsize=9.6, color="#444444", y=0.995, x=0.012, ha="left")

    fig.tight_layout(rect=(0, 0, 0.88, 0.96))
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    pdf = out.rsplit(".", 1)[0] + ".pdf"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    print(f"wrote {out} and {pdf}  (lower n={n_lo}, upper n={n_hi})")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
