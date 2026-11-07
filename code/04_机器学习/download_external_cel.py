import csv,gzip,re,urllib.request,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
base=Path('D:/Bioinfor/colorectal/GEO_ML_hub10/raw');folder=base/'GSE20916_CEL';folder.mkdir(exist_ok=True)
meta=list(csv.DictReader(open(base/'GSE20916_metadata.csv',encoding='utf8')))
rows=[r for r in meta if '(micro)' not in r['title'] and r['tissue'] in ['adenocarcinoma','normal colon']]
def fetch(r):
 gsm=r['geo_accession'];dest=folder/(gsm+'.CEL.gz')
 if dest.exists():
  try:gzip.decompress(dest.read_bytes());return gsm
  except:pass
 soft=base/'GSE20916'/(gsm+'.soft.gz');soft.parent.mkdir(exist_ok=True)
 if not soft.exists():
  b=urllib.request.urlopen(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full',timeout=60).read()
  assert b'!sample_table_end' in b;soft.write_bytes(gzip.compress(b))
 text=gzip.decompress(soft.read_bytes()).decode()
 ftp=re.search(r'!Sample_supplementary_file = (\S+\.CEL\.gz)',text).group(1)
 filename=ftp.rsplit('/',1)[1]
 url='https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/916/E-GEOD-20916/Files/'+filename.removesuffix('.gz')
 part=folder/(gsm+'.CEL.part')
 for attempt in range(30):
  offset=part.stat().st_size if part.exists() else 0
  try:
   req=urllib.request.Request(url,headers={'Range':f'bytes={offset}-'})
   with urllib.request.urlopen(req,timeout=60) as response:
    if offset and response.status!=206: offset=0
    total=offset+int(response.headers['Content-Length'])
    with open(part,'ab' if offset else 'wb') as output:
     while True:
      block=response.read(262144)
      if not block:break
      output.write(block)
   if part.stat().st_size!=total:continue
   b=part.read_bytes()
   assert len(b)>10000000 and (b.startswith(b'[CEL]') or b[:4]==b'\x40\x00\x00\x00')
   dest.write_bytes(gzip.compress(b));part.unlink();return gsm
  except Exception as e:
   print('Retry',gsm,attempt+1,str(e)[:80],flush=True)
   if attempt==29:raise
   time.sleep(3)

with ThreadPoolExecutor(max_workers=6) as pool:
 futures=[pool.submit(fetch,r) for r in rows]
 for i,f in enumerate(as_completed(futures),1):
  f.result();print('CEL completed',i,'/',len(rows),flush=True)
