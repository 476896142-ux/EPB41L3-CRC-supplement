import urllib.request,gzip,csv,time,re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
base=Path('D:/Bioinfor/colorectal/GEO_ML_hub10/raw')
series=(base/'GSE44076_series_metadata.txt').read_text(encoding='utf8')
gsms=re.findall(r'!Series_sample_id = (GSM\d+)',series)
val=list(csv.DictReader(open(base/'GSE20916_metadata.csv',encoding='utf8')))
val=[r['geo_accession'] for r in val if '(micro)' not in r['title'] and r.get('tissue') in ['adenocarcinoma','normal colon']]
jobs=[('GSE44076',g) for g in gsms]+[('GSE20916',g) for g in val]
print('Samples:',len(gsms),len(val),flush=True)
def fetch(job):
 acc,gsm=job; folder=base/acc;folder.mkdir(exist_ok=True);p=folder/(gsm+'.soft.gz')
 if p.exists():
  try:
   if '!sample_table_end' in gzip.decompress(p.read_bytes()).decode():return gsm
  except:pass
 for trial in range(4):
  try:
   url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
   data=urllib.request.urlopen(url,timeout=60).read()
   if b'!sample_table_end' not in data or len(data)<10000:raise ValueError('Incomplete sample')
   p.write_bytes(gzip.compress(data));return gsm
  except Exception:
   if trial==3:raise
   time.sleep(2+trial)
with ThreadPoolExecutor(max_workers=3) as pool:
 futures=[pool.submit(fetch,j) for j in jobs]
 for i,f in enumerate(as_completed(futures),1):
  f.result()
  if i%20==0 or i==len(jobs):print('Completed',i,'/',len(jobs),flush=True)
