import sys, csv, re, xml.etree.ElementTree as ET
from openpyxl import load_workbook

xml_path, csv_path, xlsx_path = sys.argv[1:4]

NULLS = {'', 'missing', 'not applicable', 'na', 'n/a', 'not collected',
         'not provided', 'unknown', 'none', 'not determined', 'restricted access'}

def clean(v):
    v = (v or '').strip()
    return '' if v.lower() in NULLS else v

# --- BioSample attributes ---
bios = {}
for s in ET.parse(xml_path).getroot().findall('.//BioSample'):
    acc = s.get('accession')
    ids = {i.get('db'): (i.text or '').strip() for i in s.findall('.//Id')}
    d = {}
    for a in s.findall('.//Attribute'):
        name = (a.get('harmonized_name') or a.get('attribute_name') or '').strip()
        val = clean(a.text)
        if name and val and name not in d:
            d[name] = val
    bios[acc] = {
        'biosample': acc,
        'sra': ids.get('SRA', ''),
        'strain': d.get('strain') or d.get('isolate') or d.get('sample_name') or '',
        'geo': d.get('geo_loc_name', ''),
        'isolation_source': d.get('isolation_source') or d.get('env_medium') or d.get('host_tissue_sampled') or '',
        'host': d.get('host', ''),
        'collection_date': d.get('collection_date', ''),
        'raw': d,
    }

# --- Run -> BioSample ---
run2bs, run2name = {}, {}
with open(csv_path, newline='') as fh:
    for row in csv.DictReader(fh):
        run2bs[row['Run']] = row.get('BioSample', '')
        run2name[row['Run']] = row.get('SampleName', '')

# --- panel ---
wb = load_workbook(xlsx_path)
ws = wb.active
genomes = [ws.cell(r, 1).value for r in range(2, ws.max_row + 1)]

def country(geo):
    return geo.split(':')[0].strip() if geo else ''

rows, missing = [], []
for g in genomes:
    bs = run2bs.get(g, '')
    rec = bios.get(bs)
    if rec is None:
        # study isolates: match on strain name in BioSample attributes
        hits = [v for v in bios.values() if g in (v['strain'], v['raw'].get('isolate',''), v['raw'].get('sample_name',''))]
        rec = hits[0] if hits else None
    if rec is None:
        missing.append(g)
        rows.append((g, bs, '', '', '', '', '', ''))
    else:
        rows.append((g, rec['biosample'], rec['strain'], country(rec['geo']),
                     rec['geo'], rec['isolation_source'], rec['host'],
                     rec['collection_date']))

print("panel genomes:", len(genomes))
print("resolved to a BioSample record:", sum(1 for r in rows if r[1] and r[3]))
if missing:
    print("UNRESOLVED:", missing)

print("\ncountry coverage:")
from collections import Counter
print(" ", Counter(r[3] or '(blank)' for r in rows).most_common())
print("\nisolation_source blanks:", sum(1 for r in rows if not r[5]))
print("host blanks:", sum(1 for r in rows if not r[6]))
print("collection_date blanks:", sum(1 for r in rows if not r[7]))

with open('tables/panel_biosample_meta.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['Genome','BioSample','Strain','Country','geo_loc_name',
                'Isolation source','Host','Collection date'])
    w.writerows(rows)
print("\nwrote out/panel_biosample_meta.tsv")
