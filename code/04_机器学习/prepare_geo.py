import csv,gzip,re,statistics,math
from pathlib import Path
base=Path('D:/Bioinfor/colorectal/GEO_ML_hub10/raw')
out=Path('E:/codex/2026-10-31/dui/outputs/GEO_ML_hub10');out.mkdir(exist_ok=True)
genes=['NTRK3','EPHB1','KLB','GREM1','HHIP','ANGPTL1','LAMC2','WDR72','GPER1','EPB41L3']
for acc,platform,expected in [('GSE44076','hgu219',246)]:
 ann=list(csv.DictReader(open(base/(platform+'_unique_probe_map.csv'))))
 mapping={r['PROBEID']:r['SYMBOL'] for r in ann if r['SYMBOL'] in genes}
 files=sorted((base/acc).glob('GSM*.soft.gz'))
 assert len(files)==expected,(acc,len(files),expected)
 expressions=[];metadata=[];probe_records=[]
 for path in files:
  text=gzip.decompress(path.read_bytes()).decode(); lines=text.splitlines()
  d={'Sample':path.name.split('.')[0],'Cohort':acc};values={};in_table=False
  for line in lines:
   if line.startswith('!Sample_title = '):d['Title']=line.split(' = ',1)[1]
   if line.startswith('!Sample_characteristics_ch1 = '):
    s=line.split(' = ',1)[1]
    if ': ' in s:
     k,v=s.split(': ',1);d[k.lower()]=v
   if line=='!sample_table_begin':in_table=True;continue
   if line=='!sample_table_end':in_table=False
   if in_table and '\t' in line:
    p,v=line.split('\t',1)
    if p in mapping:
     try:values[p]=float(v)
     except ValueError:pass
  if acc=='GSE44076':
   kind=d['sample type'];d['Class']='Tumor' if kind.lower() in ['tumor','tumour'] else 'Normal'
   assert kind.lower() in ['tumor','tumour','normal','mucosa'],kind
   d['Patient']=acc+'_'+d['individual id']
  else:
   assert '(micro)' not in d['Title'] and d['tissue'] in ['adenocarcinoma','normal colon']
   d['Class']='Tumor' if d['tissue']=='adenocarcinoma' else 'Normal';d['Patient']=d['Sample']
  row={'Sample':d['Sample']}
  for gene in genes:
   used=[(p,v) for p,v in values.items() if mapping[p]==gene and math.isfinite(v)]
   assert used,(acc,d['Sample'],gene,'No finite probe')
   row[gene]=statistics.median(v for p,v in used)
   for p,v in used:probe_records.append({'Sample':d['Sample'],'Gene':gene,'Probe':p,'Value':v})
  metadata.append(d);expressions.append(row)
 for name,rows in [('expression',expressions),('metadata',metadata),('probe_values',probe_records)]:
  fields=list(dict.fromkeys(k for row in rows for k in row))
  with open(out/(acc+'_'+name+'.csv'),'w',newline='',encoding='utf8') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 from collections import Counter
 print(acc,Counter(r['Class'] for r in metadata),Counter(r.get('sample type','macro') for r in metadata),flush=True)
 print('Candidate expression range',min(r[g] for r in expressions for g in genes),max(r[g] for r in expressions for g in genes))

