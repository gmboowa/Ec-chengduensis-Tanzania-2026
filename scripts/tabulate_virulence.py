#!/usr/bin/env python3
"""Tabulate the AMRFinderPlus --plus virulence calls already in the Terra table.

AMRFinderPlus ran with --plus across all 52 genomes in the TheiaProk run, so the
virulence genes are already recorded; they were never pulled into a table. This
reads the data-table export and produces the per-genome and summary tables the
Results need, restricted to the 43 retained genomes.

  python3 tabulate_virulence.py <Ec_all52 export.tsv> <Table_S1_merged...xlsx>
"""
import csv
import sys
from collections import Counter

from openpyxl import load_workbook

export, panel_xlsx = sys.argv[1], sys.argv[2]

ws = load_workbook(panel_xlsx)["Table S1"]
hdr = [c.value for c in ws[1]]
keep = set()
for i in range(2, ws.max_row + 1):
    row = dict(zip(hdr, [ws.cell(i, j).value for j in range(1, ws.max_column + 1)]))
    if row["Retained in final panel"] == "Yes":
        keep.add(row["Genome"])

with open(export, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

cols = rows[0].keys()
id_col = next(c for c in cols if c.endswith("_id") or c == "entity:Ec_all52_id")
vir_cols = [c for c in cols if "virulence" in c.lower()]
if not vir_cols:
    sys.exit(f"no virulence column in the export; columns were:\n  " +
             "\n  ".join(sorted(cols)))
vir = vir_cols[0]
print(f"using column: {vir}")

NULLS = {"", "none", "na", "n/a"}
SENTINEL = "no virulence genes detected"
per_genome, counts = [], Counter()
for r in rows:
    g = r[id_col]
    if g not in keep:
        continue
    raw = (r.get(vir) or "").strip()
    if raw.lower() in NULLS or SENTINEL in raw.lower():
        genes = []
    else:
        genes = [x.strip() for x in raw.replace(";", ",").split(",") if x.strip()]
    per_genome.append((g, len(genes), ",".join(sorted(genes)) or "-"))
    for x in genes:
        counts[x] += 1

print(f"retained genomes found in the export: {len(per_genome)}/{len(keep)}")
missing = keep - {g for g, _, _ in per_genome}
if missing:
    print("NOT IN EXPORT:", sorted(missing))

with open("out/virulence_per_genome_43.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome", "virulence_genes_n", "virulence_genes"])
    w.writerows(sorted(per_genome))

with open("out/virulence_gene_counts_43.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene", "genomes", "percent_of_43"])
    for gene, n in counts.most_common():
        w.writerow([gene, n, f"{100 * n / len(per_genome):.1f}"])

print("\ngene frequencies across the retained panel:")
for gene, n in counts.most_common(25):
    print(f"  {gene:<16}{n:>3}  ({100 * n / len(per_genome):.0f}%)")
print(f"\ngenomes with no virulence gene detected: "
      f"{sum(1 for _, n, _ in per_genome if n == 0)}")
print("wrote out/virulence_per_genome_43.tsv and out/virulence_gene_counts_43.tsv")
