#!/usr/bin/env python3
"""Generate iTOL annotation files from the panel metadata.

  python3 make_itol.py <matrix.tsv> <Table_S1.xlsx> <virulence.tsv> <outdir>

Writes one file per track. Drag them onto the tree in iTOL; no manual entry.
"""
import csv
import itertools
import sys
from pathlib import Path

from openpyxl import load_workbook

OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
         "#F0E442", "#999999", "#882255", "#44AA99", "#117733", "#DDCC77",
         "#AA4499", "#6699CC", "#661100"]
REF = "Reference"


def main(argv):
    matrix, xlsx, vir_tsv, outdir = argv[1:5]
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)

    rows = list(csv.reader(open(matrix, newline=""), delimiter="\t"))
    hdr = rows[0][1:]
    D = {r[0]: {h: int(v) for h, v in zip(hdr, r[1:])} for r in rows[1:] if r}
    names = [n for n in hdr if n != REF]

    ws = load_workbook(xlsx)["Table S1"]
    th = [c.value for c in ws[1]]
    meta = {}
    for i in range(2, ws.max_row + 1):
        d = dict(zip(th, [ws.cell(i, j).value for j in range(1, ws.max_column + 1)]))
        meta[d["Genome"]] = d
    vir = {r["genome"]: r["virulence_genes"]
           for r in csv.DictReader(open(vir_tsv, newline=""), delimiter="\t")}

    def colorstrip(fname, label, value_of, palette_offset=0, width=28):
        vals = sorted({value_of(n) for n in names})
        cmap = {v: OKABE[(i + palette_offset) % len(OKABE)] for i, v in enumerate(vals)}
        for k in ("not provided", "n/a", "unknown"):
            if k in cmap:
                cmap[k] = "#E8E8E8"
        with open(out / fname, "w") as fh:
            fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
            fh.write(f"DATASET_LABEL\t{label}\nCOLOR\t#333333\n")
            fh.write(f"STRIP_WIDTH\t{width}\nMARGIN\t4\nBORDER_WIDTH\t1\n")
            fh.write("BORDER_COLOR\t#ffffff\nSHOW_STRIP_LABELS\t0\n")
            fh.write("LEGEND_TITLE\t" + label + "\n")
            fh.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in vals) + "\n")
            fh.write("LEGEND_COLORS\t" + "\t".join(cmap[v] for v in vals) + "\n")
            fh.write("LEGEND_LABELS\t" + "\t".join(str(v) for v in vals) + "\n")
            fh.write("DATA\n")
            for n in names:
                fh.write(f"{n}\t{cmap[value_of(n)]}\t{value_of(n)}\n")
        print(f"  {fname}  ({len(vals)} categories)")

    colorstrip("itol_country.txt", "Country",
               lambda n: meta.get(n, {}).get("Country") or "n/a")
    colorstrip("itol_source.txt", "Isolation source",
               lambda n: (meta.get(n, {}).get("Isolation source") or "not provided").lower(),
               palette_offset=4)

    def amr(n):
        v = meta.get(n, {}).get("AMR genes (n)")
        if v is None:
            return "unknown"
        return "low (<=5)" if v <= 5 else ("medium (6-8)" if v <= 8 else "high (>=9)")

    order = ["low (<=5)", "medium (6-8)", "high (>=9)", "unknown"]
    shades = {"low (<=5)": "#BBDDEE", "medium (6-8)": "#6699CC",
              "high (>=9)": "#114477", "unknown": "#E8E8E8"}
    present = [v for v in order if any(amr(n) == v for n in names)]
    with open(out / "itol_amr_load.txt", "w") as fh:
        fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tAMR load\n")
        fh.write("COLOR\t#114477\nSTRIP_WIDTH\t28\nMARGIN\t4\nBORDER_WIDTH\t1\n")
        fh.write("BORDER_COLOR\t#ffffff\nSHOW_STRIP_LABELS\t0\n")
        fh.write("LEGEND_TITLE\tAcquired resistance genes\n")
        fh.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in present) + "\n")
        fh.write("LEGEND_COLORS\t" + "\t".join(shades[v] for v in present) + "\n")
        fh.write("LEGEND_LABELS\t" + "\t".join(present) + "\n")
        fh.write("DATA\n")
        for n in names:
            fh.write(f"{n}\t{shades[amr(n)]}\t{amr(n)}\n")
    print("  itol_amr_load.txt")

    with open(out / "itol_binary.txt", "w") as fh:
        fh.write("DATASET_BINARY\nSEPARATOR TAB\n")
        fh.write("DATASET_LABEL\tPlasmid replicon / ybtPQ\nCOLOR\t#333333\n")
        fh.write("FIELD_SHAPES\t2\t2\nFIELD_LABELS\tPlasmid replicon\tybtP/ybtQ\n")
        fh.write("FIELD_COLORS\t#0072B2\t#D55E00\n")
        fh.write("SHOW_LABELS\t1\nHEIGHT_FACTOR\t1\nMARGIN\t6\nDATA\n")
        for n in names:
            rep = meta.get(n, {}).get("Plasmid replicons") or ""
            has_rep = 1 if (rep and "No plasmid" not in rep) else 0
            has_ybt = 1 if "ybt" in (vir.get(n, "") or "") else 0
            fh.write(f"{n}\t{has_rep}\t{has_ybt}\n")
    print("  itol_binary.txt")

    # label styling for the isolates from this study, and clade ranges
    study = [n for n in names if meta.get(n, {}).get("Origin") == "This study"]
    hi = [n for n in names if (meta.get(n, {}).get("ANI to type strain (%)") or 0) >= 99]
    lo = [n for n in names
          if n in meta and (meta[n].get("ANI to type strain (%)") or 0) < 99]

    def span(group):
        """Two tips whose common ancestor is the whole cluster."""
        a, b = max(itertools.combinations(group, 2), key=lambda p: D[p[0]][p[1]])
        return a, b

    with open(out / "itol_styles.txt", "w") as fh:
        fh.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        lo_a, lo_b = span(lo)
        hi_a, hi_b = span(hi)
        fh.write(f"{lo_a}|{lo_b}\trange\t#EAF2F8\tLower ANI cluster "
                 f"(95.67-95.80%), n={len(lo)}\n")
        fh.write(f"{hi_a}|{hi_b}\trange\t#FDF2EC\tUpper ANI cluster "
                 f"(99.03-99.99%), n={len(hi)}\n")
        for n in study:
            fh.write(f"{n}\tlabel\t#B2182B\tbold\n")
            fh.write(f"{n}\tbranch\t#B2182B\tnormal\t2\n")
    print(f"  itol_styles.txt  (study isolates: {', '.join(study)})")

    with open(out / "README_itol.txt", "w") as fh:
        fh.write(
            "Load Ec_43_phylo_iqtree.nwk in iTOL, then drag each file onto the tree.\n\n"
            "  itol_country.txt    country of origin\n"
            "  itol_source.txt     isolation source as recorded in NCBI BioSample\n"
            "  itol_amr_load.txt   acquired resistance genes, banded\n"
            "  itol_binary.txt     plasmid replicon and ybtP/ybtQ presence\n"
            "  itol_styles.txt     cluster shading and red labels for this study\n\n"
            "Branch lengths differ by about three orders of magnitude between and\n"
            "within clusters, so to scale the within-cluster structure is not legible.\n"
            "Either switch to Normal (cladogram) display, or collapse the upper\n"
            "cluster, and give the divergence in the caption; the distance\n"
            "distribution figure carries the quantitative comparison.\n")
    print("  README_itol.txt")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
