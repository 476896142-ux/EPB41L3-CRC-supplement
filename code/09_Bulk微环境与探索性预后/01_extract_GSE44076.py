"""从已有完整 GEO 表提取 MCP-counter 所需标记及 EPB41L3，不重新下载。"""
from pathlib import Path
import csv, gzip, statistics, math
root = Path('D:/Bioinfor/colorectal/EPB41L3_bulk_validation')
raw = Path('D:/Bioinfor/colorectal/GEO_ML_hub10/raw')
with (root/'reference/genes.txt').open() as f:
    genes = {r['HUGO symbols'] for r in csv.DictReader(f, delimiter='\t')} | {'EPB41L3'}
with (raw/'hgu219_unique_probe_map.csv').open() as f:
    mapping = {r['PROBEID']:r['SYMBOL'] for r in csv.DictReader(f) if r['SYMBOL'] in genes}
rows=[]
for path in sorted((raw/'GSE44076').glob('GSM*.soft.gz')):
    values={}; inside=False
    for line in gzip.open(path,'rt'):
        line=line.strip()
        if line=='!sample_table_begin': inside=True; continue
        if line=='!sample_table_end': inside=False
        if inside and '\t' in line:
            probe,value=line.split('\t',1)
            if probe in mapping:
                value=float(value)
                if math.isfinite(value): values.setdefault(mapping[probe],[]).append(value)
    rows.append({'Sample':path.name.split('.')[0], **{g:statistics.median(v) for g,v in values.items()}})
assert len(rows)==246
with (root/'GSE44076_marker_expression.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['Sample']+sorted(set(mapping.values())))
    writer.writeheader();writer.writerows(rows)
print('Extracted',len(rows),'samples;',len(mapping),'mapped probes;',len(set(mapping.values())),'genes')
