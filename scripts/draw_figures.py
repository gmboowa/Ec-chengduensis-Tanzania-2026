#!/usr/bin/env python3
"""Figure 1 (single tree with a scale break) and the pairwise-distance figure.

  python3 draw_figures.py <tree.nwk> <matrix.tsv> <Table_S1.xlsx> <virulence.tsv> <outdir>

The branch separating the two clusters is ~1,000x longer than any branch within
either. It is drawn with an axis break so both the separation and the internal
structure are visible; the break is marked and the compression stated. No other
branch is altered.
"""
import csv
import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from draw_tree import (Node, parse_newick, leaves, nodes, load_metadata,
                       midpoint_root, amr_band, OKABE, GREYS)

REF = "Reference"


def strip_ref(node):
    if node.is_leaf:
        return None if node.name == REF else node
    kids = [k for k in (strip_ref(c) for c in node.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].length += node.length
        return kids[0]
    new = Node()
    new.name, new.length, new.support = node.name, node.length, node.support
    new.children = kids
    return new


def figure1(nwk, xlsx, vir_tsv, out):
    root = midpoint_root(strip_ref(parse_newick(open(nwk).read())))
    meta, vir = load_metadata(xlsx, vir_tsv)

    tips = list(leaves(root))
    # order tips so the two clusters stay contiguous
    for i, t in enumerate(tips):
        t.y = i

    # Midpoint rooting splits the long inter-cluster branch into two halves, so
    # every branch over the cutoff is compressed, not just the single longest.
    lens = sorted((n.length for n in nodes(root) if n is not root), reverse=True)
    rest = next((v for v in lens if v < lens[0] * 0.25), lens[-1])
    CUT = rest * 2.5
    KEEP = rest * 1.4
    broken = set()

    def place(n, x0):
        eff = n.length
        if n.length > CUT:
            eff = KEEP
            broken.add(n)
        n.x = x0 + eff
        for c in n.children:
            place(c, n.x)
    place(root, 0.0)
    for n in sorted(nodes(root), key=lambda z: -z.x):
        if not n.is_leaf:
            n.y = sum(c.y for c in n.children) / len(n.children)
    xmax = max(n.x for n in nodes(root))
    deep = max(broken, key=lambda n: n.length) if broken else None
    parents = {c: p for p in nodes(root) for c in p.children}

    fig, ax = plt.subplots(figsize=(14.5, 0.30 * len(tips) + 3.0), dpi=200)

    def draw(node, px):
        ax.plot([px, node.x], [node.y, node.y], color="#2b2b2b", lw=1.0,
                solid_capstyle="butt", zorder=2)
        if node.children:
            ys = [c.y for c in node.children]
            ax.plot([node.x, node.x], [min(ys), max(ys)], color="#2b2b2b", lw=1.0,
                    solid_capstyle="butt", zorder=2)
            for c in node.children:
                draw(c, node.x)
    draw(root, 0.0)

    # break marks across the compressed branch
    for n in broken:
        px = parents[n].x if n in parents else 0.0
        bx = (px + n.x) / 2
        for off in (-xmax * 0.005, xmax * 0.005):
            ax.plot([bx + off - xmax * 0.0035, bx + off + xmax * 0.0035],
                    [n.y - 0.55, n.y + 0.55], color="#2b2b2b", lw=1.2, zorder=5)
    if deep is not None:
        px = parents[deep].x if deep in parents else 0.0
        ax.annotate("branches marked // are compressed:\nthe two clusters differ by "
                    "121,375-122,118 SNPs",
                    xy=((px + deep.x) / 2, deep.y),
                    xytext=((px + deep.x) / 2, len(tips) + max(1.6, len(tips) * 0.05)),
                    ha="center", fontsize=8.6, color="#555555",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.8))

    for n in nodes(root):
        if n.children and n.support is not None and n.support >= 90:
            ax.plot(n.x, n.y, "o", ms=3.0, mfc="#2b2b2b", mec="none", zorder=4)

    label_x = xmax * 1.02
    col_w, gap = xmax * 0.021, xmax * 0.008
    track_x = label_x + xmax * 0.21

    countries = sorted({meta.get(t.name, {}).get("Country") or "n/a" for t in tips})
    cmap = {c: OKABE[i % len(OKABE)] for i, c in enumerate(countries)}
    sources = sorted({(meta.get(t.name, {}).get("Isolation source") or "not provided").lower()
                      for t in tips})
    smap = {s: OKABE[(i + 4) % len(OKABE)] for i, s in enumerate(sources)}
    smap["not provided"] = "#E8E8E8"
    amrs = ["low (<=5)", "medium (6-8)", "high (>=9)", "unknown"]
    amap = dict(zip(amrs, ["#BBDDEE", "#6699CC", "#114477", "#E8E8E8"]))

    for t in tips:
        m = meta.get(t.name, {})
        study = m.get("Origin") == "This study"
        ax.text(label_x, t.y, t.name, va="center", ha="left", fontsize=8.0,
                fontweight="bold" if study else "normal",
                color="#B2182B" if study else "#222222")
        x = track_x
        rep = (m.get("Plasmid replicons") or "")
        cells = [cmap.get(m.get("Country") or "n/a", "#E8E8E8"),
                 smap.get((m.get("Isolation source") or "not provided").lower(), "#E8E8E8"),
                 amap[amr_band(m.get("AMR genes (n)"))]]
        for colour in cells:
            ax.add_patch(Rectangle((x, t.y - 0.40), col_w, 0.80, color=colour,
                                   lw=0.3, ec="white", zorder=3))
            x += col_w + gap
        for present in (bool(rep) and "No plasmid" not in rep,
                        "ybt" in (vir.get(t.name, "") or "")):
            ax.add_patch(Rectangle((x, t.y - 0.40), col_w, 0.80,
                                   color=GREYS["present"] if present else GREYS["absent"],
                                   lw=0.5, ec="#888888", zorder=3))
            x += col_w + gap

    for i, h in enumerate(["Country", "Source", "AMR load", "Plasmid", "ybtPQ"]):
        ax.text(track_x + i * (col_w + gap) + col_w / 2, len(tips) - 0.3, h,
                rotation=45, ha="left", va="bottom", fontsize=7.6, color="#555555")

    import math
    raw = (xmax - KEEP * len(broken)) / 4 if broken else xmax / 4
    raw = max(raw, xmax * 0.06)
    mag = 10 ** math.floor(math.log10(raw))
    bar = next((m * mag for m in (1, 2, 5, 10) if raw <= m * mag), 10 * mag)
    y0 = -2.2
    ax.plot([0, bar], [y0, y0], color="#2b2b2b", lw=1.5)
    for xx in (0, bar):
        ax.plot([xx, xx], [y0 - 0.3, y0 + 0.3], color="#2b2b2b", lw=1.5)
    ax.text(bar / 2, y0 - 0.9, f"{bar:g} substitutions per site  (applies either side "
            f"of the break)", ha="left", va="top", fontsize=8.2, color="#333333")

    handles = [Line2D([], [], marker="s", ls="", ms=7, mfc=cmap[c], mec="none", label=c)
               for c in countries]
    l1 = ax.legend(handles=handles, title="Country", loc="upper left",
                   bbox_to_anchor=(1.26, 1.0), frameon=False, fontsize=7.5,
                   title_fontsize=8.6, labelspacing=0.3)
    ax.add_artist(l1)
    h2 = [Line2D([], [], marker="s", ls="", ms=7, mfc=amap[x], mec="none", label=x)
          for x in amrs[:3]]
    h2 += [Line2D([], [], marker="o", ls="", ms=5, mfc="#2b2b2b", mec="none",
                  label="bootstrap >= 90"),
           Line2D([], [], marker="s", ls="", ms=7, mfc="#333333", mec="none",
                  label="feature present"),
           Line2D([], [], marker="s", ls="", ms=7, mfc="white", mec="#888888",
                  label="feature absent")]
    ax.legend(handles=h2, title="Annotation", loc="upper left",
              bbox_to_anchor=(1.26, 0.45), frameon=False, fontsize=7.5,
              title_fontsize=8.6, labelspacing=0.3)

    ax.set_xlim(-xmax * 0.02, track_x + 5 * (col_w + gap) + xmax * 0.02)
    ax.set_ylim(y0 - 2.0, len(tips) + max(4.0, len(tips) * 0.11))
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight", facecolor="white")
    print(f"wrote {out}.png / .pdf  ({len(tips)} tips)")


def figure_distances(matrix_tsv, xlsx, out):
    rows = list(csv.reader(open(matrix_tsv, newline=""), delimiter="\t"))
    hdr = rows[0][1:]
    D = {r[0]: {h: int(v) for h, v in zip(hdr, r[1:])} for r in rows[1:] if r}
    names = [n for n in hdr if n != REF]

    from openpyxl import load_workbook
    ws = load_workbook(xlsx)["Table S1"]
    th = [c.value for c in ws[1]]
    ani = {}
    for i in range(2, ws.max_row + 1):
        d = dict(zip(th, [ws.cell(i, j).value for j in range(1, ws.max_column + 1)]))
        if d["Retained in final panel"] == "Yes":
            ani[d["Genome"]] = d["ANI to type strain (%)"]

    hi = [n for n in names if ani.get(n, 0) >= 99]
    lo = [n for n in names if n in ani and ani[n] < 99]
    within = [D[a][b] for a, b in itertools.combinations(hi, 2)] + \
             [D[a][b] for a, b in itertools.combinations(lo, 2)]
    between = [D[a][b] for a in hi for b in lo]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.2, 6.4), dpi=200,
                                   gridspec_kw={"height_ratios": [2, 1.15]})

    bins = range(0, 130001, 2000)
    ax1.hist(within, bins=bins, color="#0072B2", alpha=0.85, label="within a cluster")
    ax1.hist(between, bins=bins, color="#D55E00", alpha=0.85, label="between clusters")
    ax1.set_yscale("log")
    ax1.set_xlabel("pairwise core-genome SNP distance")
    ax1.set_ylabel("number of genome pairs (log)")
    ax1.axvspan(max(within), min(between), color="#999999", alpha=0.13, lw=0)
    ax1.text((max(within) + min(between)) / 2, ax1.get_ylim()[1] * 0.45,
             f"no pair falls in this range\n{max(within):,} to {min(between):,} SNPs",
             ha="center", va="center", fontsize=9, color="#444444")
    ax1.legend(frameon=False, fontsize=9)
    ax1.set_title("A   All pairwise core-genome SNP distances", loc="left",
                  fontsize=10.5, pad=8)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    vals = sorted(ani[n] for n in names if n in ani)
    ax2.hist(vals, bins=[95.5 + 0.05 * i for i in range(95)], color="#009E73", alpha=0.9)
    ax2.set_xlabel("ANI to the type-strain assembly GCF_001984825.2 (%)")
    ax2.set_ylabel("genomes")
    ax2.axvspan(max(v for v in vals if v < 99), min(v for v in vals if v >= 99),
                color="#999999", alpha=0.13, lw=0)
    ax2.set_title("B   ANI of the same genomes to the type strain", loc="left",
                  fontsize=10.5, pad=8)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    fig.suptitle("Both measures are taken against the same reference genome and are "
                 "therefore concordant rather than independent.",
                 fontsize=8.8, color="#555555", y=0.015, x=0.012, ha="left")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight", facecolor="white")
    print(f"wrote {out}.png / .pdf  (within n={len(within)}, between n={len(between)})")


if __name__ == "__main__":
    nwk, mat, xlsx, vir, outdir = sys.argv[1:6]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    figure1(nwk, xlsx, vir, f"{outdir}/Figure1_phylogeny")
    figure_distances(mat, xlsx, f"{outdir}/FigureS_distance_distribution")
