import sys, csv, re
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

meta_p, qc_p, out_p = sys.argv[1:4]

meta = {}
with open(meta_p, newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        meta[row['Genome']] = row

wb_in = load_workbook(qc_p)
ws_in = wb_in.active
qc_hdr = [c.value for c in ws_in[1]]
qc = {}
for r in ws_in.iter_rows(min_row=2, values_only=True):
    qc[r[0]] = dict(zip(qc_hdr, r))

# run accession per genome, from the search CSV via BioSample
runs = {}
with open(sys.argv[4], newline='') as fh:
    for row in csv.DictReader(fh):
        runs[row['BioSample']] = row['Run']

COLS = [
    ('Genome', 12), ('Origin', 13), ('Run accession', 14), ('BioSample', 15),
    ('Strain', 16), ('Country', 14), ('Geographic location', 22),
    ('Isolation source', 22), ('Host', 16), ('Collection date', 14),
    ('GAMBIT species', 26), ('ANI to type strain (%)', 13),
    ('ANI group', 18), ('Contigs (n)', 10), ('BUSCO complete (%)', 12),
    ('Coverage (x)', 11), ('AMR genes (n)', 11), ('AMR genes', 54),
    ('Plasmid replicons', 26), ('QC flags', 30), ('Retained in final panel', 13),
]

ITAL_GENE_COL = 'AMR genes'
ARIAL = 'Arial'

wb = Workbook()

# ---------- legend sheet ----------
lg = wb.active
lg.title = 'Legend'
lines = [
    ('Supplementary Table S1.', True),
    ('Genome panel, isolate metadata, assembly quality metrics, taxonomic assignment '
     'and inclusion status for the 52 Enterobacter chengduensis-candidate genomes '
     'analysed in this study.', False),
    ('', False),
    ('This table consolidates the genome panel, isolate metadata, assembly '
     'quality metrics, taxonomic assignment and inclusion status for all 52 '
     'genomes retrieved, including those excluded from the final analysis '
     'panel.', False),
    ('', False),
    ('Column definitions', True),
    ('Genome - isolate identifier; for public genomes, the SRA run accession.', False),
    ('Origin - "This study" for the two Tanzanian isolates deposited under BioProject '
     'PRJNA1447346; "Public genome" for genomes retrieved from the SRA.', False),
    ('Run accession / BioSample - NCBI accessions, as retrieved on 11 November 2026.', False),
    ('Country, Geographic location, Isolation source, Host, Collection date - as recorded '
     'in the submitter-provided NCBI BioSample attributes. Blank cells indicate the '
     'attribute was not populated, or was recorded as "missing" or "not applicable", '
     'by the original submitter; these were not inferred.', False),
    ('GAMBIT species - species assignment from GAMBIT (TheiaProk_Illumina_PE_PHB v4.3.0).', False),
    ('ANI to type strain (%) - average nucleotide identity computed with MUMmer within the '
     'TheiaProk workflow against the Enterobacter chengduensis type strain WCHECl-C4, '
     'assembly GCF_001984825.2 (ASM198482v2).', False),
    ('ANI group - ≥99%, "divergent (95–99%)", or "other species" (<95%).', False),
    ('Contigs (n), BUSCO complete (%), Coverage (x) - assembly metrics from QUAST, BUSCO '
     'and the read-mapping step of the same workflow.', False),
    ('AMR genes - acquired resistance determinants detected by AMRFinderPlus 4.2.7; gene '
     'names are italicised.', False),
    ('Plasmid replicons - PlasmidFinder 2.1.6.', False),
    ('QC flags - reasons for exclusion, where applicable.', False),
    ('Retained in final panel - "Yes" if the genome passed the ANI (≥95%) and GAMBIT '
     'species-assignment criteria; "No" otherwise.', False),
]
r = 1
for text, bold in lines:
    c = lg.cell(r, 1, text)
    c.font = Font(name=ARIAL, size=11, bold=bold)
    c.alignment = Alignment(wrap_text=True, vertical='top')
    lg.row_dimensions[r].height = None if not text else max(15, 15 * (1 + len(text) // 110))
    r += 1
lg.column_dimensions['A'].width = 118

# ---------- data sheet ----------
ws = wb.create_sheet('Table S1')
head_fill = PatternFill('solid', fgColor='DCE6F1')
thin = Side(style='thin', color='BFBFBF')
for j, (name, width) in enumerate(COLS, 1):
    c = ws.cell(1, j, name)
    c.font = Font(name=ARIAL, size=10, bold=True)
    c.fill = head_fill
    c.alignment = Alignment(wrap_text=True, vertical='bottom', horizontal='center')
    c.border = Border(bottom=thin)
    ws.column_dimensions[get_column_letter(j)].width = width
ws.row_dimensions[1].height = 32

def gene_cell(cell, value):
    """AMR gene lists are entirely gene names -> italicise the whole cell."""
    cell.value = value
    cell.font = Font(name=ARIAL, size=10, italic=True)

order = sorted(qc, key=lambda g: (qc[g]['Origin'] != 'This study', g))
for i, g in enumerate(order, 2):
    m, q = meta[g], qc[g]
    flags = q['QC flags'] or '-'
    retained = 'No' if flags not in ('-', '', None) else 'Yes'
    vals = [
        g, q['Origin'], runs.get(m['BioSample'], ''), m['BioSample'], m['Strain'],
        m['Country'], m['geo_loc_name'], m['Isolation source'], m['Host'],
        m['Collection date'], q['GAMBIT species'], q['ANI (%)'], q['ANI group'],
        q['Contigs'], q['BUSCO complete (%)'], q['Coverage (x)'],
        q['AMR genes (n)'], q['AMR genes'], q['Plasmid replicons'], flags, retained,
    ]
    for j, (v, (name, _)) in enumerate(zip(vals, COLS), 1):
        cell = ws.cell(i, j)
        if name in (ITAL_GENE_COL, 'GAMBIT species'):
            gene_cell(cell, v)
        else:
            cell.value = v
            cell.font = Font(name=ARIAL, size=10)
        cell.alignment = Alignment(vertical='top',
                                   wrap_text=name in ('AMR genes', 'QC flags',
                                                      'Plasmid replicons',
                                                      'Geographic location',
                                                      'Isolation source'))
        if name in ('ANI to type strain (%)',
                    'BUSCO complete (%)', 'Coverage (x)'):
            cell.number_format = '0.00'

ws.freeze_panes = 'C2'
ws.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{ws.max_row}"

wb.save(out_p)
print("rows:", ws.max_row - 1, "cols:", ws.max_column)
print("retained Yes:", sum(1 for i in range(2, ws.max_row+1)
                           if ws.cell(i, len(COLS)).value == 'Yes'))
print("saved", out_p)
