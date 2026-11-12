#!/usr/bin/env python3
"""Validate the iTOL bundle against the tree it annotates."""
import re
import sys
from pathlib import Path

from draw_tree import parse_newick, leaves

d = Path(sys.argv[1])
problems = []

nwk = d / "Ec_43_phylo_iqtree.nwk"
root = parse_newick(nwk.read_text())
tips = [t.name for t in leaves(root)]
print(f"tree: {len(tips)} tips, {len(set(tips))} unique")
if len(tips) != len(set(tips)):
    problems.append("duplicate tip names")
if "Reference" in tips:
    problems.append("Reference still present")
if not nwk.read_text().strip().endswith(";"):
    problems.append("newick does not end with ';'")
blank = [t for t in tips if not t or " " in t]
if blank:
    problems.append(f"tip names with spaces or blank: {blank}")

HEX = re.compile(r"^#[0-9A-Fa-f]{6}$")

for f in sorted(d.glob("itol_*.txt")):
    text = f.read_text()
    lines = [l.rstrip("\n") for l in text.splitlines()]
    kind = lines[0].strip()
    if "SEPARATOR TAB" not in text:
        problems.append(f"{f.name}: no SEPARATOR TAB")
    try:
        di = lines.index("DATA")
    except ValueError:
        problems.append(f"{f.name}: no DATA line")
        continue
    data = [l for l in lines[di + 1:] if l.strip()]
    ids = [l.split("\t")[0] for l in data]

    if kind == "TREE_COLORS":
        for row in data:
            parts = row.split("\t")
            target = parts[0]
            for one in target.split("|"):
                if one not in tips:
                    problems.append(f"{f.name}: '{one}' is not a tip")
            if parts[1] not in ("range", "label", "branch", "clade"):
                problems.append(f"{f.name}: unknown type {parts[1]}")
            if not HEX.match(parts[2]):
                problems.append(f"{f.name}: bad colour {parts[2]}")
        print(f"{f.name}: {len(data)} directives, targets resolve")
        continue

    missing = sorted(set(tips) - set(ids))
    extra = sorted(set(ids) - set(tips))
    if missing:
        problems.append(f"{f.name}: {len(missing)} tips unannotated: {missing[:4]}")
    if extra:
        problems.append(f"{f.name}: {len(extra)} ids not in tree: {extra[:4]}")
    if len(ids) != len(set(ids)):
        problems.append(f"{f.name}: duplicate ids")

    for row in data:
        for tok in row.split("\t")[1:]:
            if tok.startswith("#") and not HEX.match(tok):
                problems.append(f"{f.name}: bad colour {tok}")

    if kind == "DATASET_COLORSTRIP":
        leg = {k: lines[[i for i, l in enumerate(lines)
                         if l.startswith(k)][0]].split("\t")[1:]
               for k in ("LEGEND_SHAPES", "LEGEND_COLORS", "LEGEND_LABELS")
               if any(l.startswith(k) for l in lines)}
        if leg and len({len(v) for v in leg.values()}) != 1:
            problems.append(f"{f.name}: legend arrays differ in length "
                            f"{ {k: len(v) for k, v in leg.items()} }")
        used = {r.split("\t")[1] for r in data}
        declared = set(leg.get("LEGEND_COLORS", []))
        if declared and not used <= declared:
            problems.append(f"{f.name}: colours used but not in legend: "
                            f"{sorted(used - declared)}")
        print(f"{f.name}: {len(data)} rows, "
              f"{len(leg.get('LEGEND_LABELS', []))} legend entries")

    if kind == "DATASET_BINARY":
        n_sh = len([l for l in lines if l.startswith("FIELD_SHAPES")][0].split("\t")) - 1
        n_lb = len([l for l in lines if l.startswith("FIELD_LABELS")][0].split("\t")) - 1
        n_co = len([l for l in lines if l.startswith("FIELD_COLORS")][0].split("\t")) - 1
        if not (n_sh == n_lb == n_co):
            problems.append(f"{f.name}: FIELD_* counts differ {n_sh}/{n_lb}/{n_co}")
        for row in data:
            vals = row.split("\t")[1:]
            if len(vals) != n_sh:
                problems.append(f"{f.name}: '{row.split(chr(9))[0]}' has {len(vals)} "
                                f"values, expected {n_sh}")
            if any(v not in ("0", "1", "-1") for v in vals):
                problems.append(f"{f.name}: non-binary value in {row[:40]}")
        print(f"{f.name}: {len(data)} rows, {n_sh} fields")

print()
if problems:
    print(f"{len(problems)} PROBLEM(S):")
    for p in problems:
        print("  -", p)
    raise SystemExit(1)
print("all checks pass")
