#!/usr/bin/env python3
"""Add GTDB-Tk ani_rep columns to Supplementary Table S1."""
import csv, openpyxl
from copy import copy
from openpyxl.styles import Font

SRC = 'tables/Table_S1_merged_genome_panel.xlsx'
ANI = 'results/gtdbtk/gtdbtk_ani_closest.tsv'
OUT = 'tables/Table_S1_merged_genome_panel.xlsx'

gt = {}
for r in csv.DictReader(open(ANI), delimiter='\t'):
    if not r.get('user_genome'):
        continue
    gt[r['user_genome'].replace('_filtered_contigs', '')] = r

wb = openpyxl.load_workbook(SRC)
ws = wb['Table S1']
hdr = [c.value for c in ws[1]]

NEW = ['GTDB closest representative', 'GTDB species',
       'GTDB ANI (%)', 'GTDB AF', 'GTDB species circumscription met']
if NEW[0] in hdr:
    raise SystemExit('GTDB columns already present')

# insert before the last two columns so QC flags / retention stay rightmost
anchor = hdr.index('QC flags')
ws.insert_cols(anchor + 1, len(NEW))

style_src = ws.cell(row=1, column=1)
for n, name in enumerate(NEW):
    c = ws.cell(row=1, column=anchor + 1 + n, value=name)
    c.font = copy(style_src.font)
    c.fill = copy(style_src.fill)
    c.border = copy(style_src.border)
    c.alignment = copy(style_src.alignment)

body = ws.cell(row=2, column=1)
missing = []
for row in range(2, ws.max_row + 1):
    genome = ws.cell(row=row, column=1).value
    g = gt.get(genome)
    if g is None:
        missing.append(genome)
        vals = ['', '', '', '', '']
    else:
        vals = [g['reference_genome'],
                g['reference_taxonomy'].split('s__')[1],
                float(g['skani_ani']),
                float(g['skani_af']),
                'Yes' if g['satisfies_gtdb_circumscription_criteria'] == 'True' else 'No']
    for n, v in enumerate(vals):
        c = ws.cell(row=row, column=anchor + 1 + n, value=v)
        c.font = copy(body.font)
        c.border = copy(body.border)
        c.alignment = copy(body.alignment)
        if n == 1 and v:
            c.font = Font(name=body.font.name, sz=body.font.sz, italic=True)
        if n == 2 and v != '':
            c.number_format = '0.00'
        if n == 3 and v != '':
            c.number_format = '0.000'

if missing:
    raise SystemExit('no GTDB result for: %s' % missing)

for n, w in enumerate([26, 24, 13, 11, 16]):
    ws.column_dimensions[ws.cell(row=1, column=anchor + 1 + n).column_letter].width = w

# ---- legend ----
lg = wb['Legend']
ref_cell = lg.cell(row=12, column=1)
notes = [
 'GTDB closest representative — accession of the closest GTDB species representative '
 'identified by GTDB-Tk v2.4.0 ani_rep (GTDB release R220) using skani.',
 'GTDB species — species of that representative under the GTDB taxonomy.',
 'GTDB ANI (%) and GTDB AF — skani average nucleotide identity and alignment fraction '
 'to the closest representative. These are not directly comparable to the MUMmer-based '
 '"ANI to type strain (%)" column, which uses a different aligner.',
 'GTDB species circumscription met — whether the genome falls within the ANI and alignment-'
 'fraction circumscription radius of that species cluster, as reported by GTDB-Tk.',
]
start = lg.max_row + 1
for n, t in enumerate(notes):
    c = lg.cell(row=start + n, column=1, value=t)
    c.font = copy(ref_cell.font)
    c.alignment = copy(ref_cell.alignment)

wb.save(OUT)
print('wrote', OUT)
print('columns now', ws.max_column, 'rows', ws.max_row)
