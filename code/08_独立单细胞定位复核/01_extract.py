"""GSE146771 10x原作者注释定位复核：下载、条码匹配及固定标志提取。"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,gzip,json,hashlib,time
import pandas as pd
import numpy as np

out=Path(__file__).parent
src=out.parent/'cohort_pairing_audit/Miro_atlas_SPP1_audit'
name='GSE146771_CRC.Leukocyte.10x.TPM.txt.gz'
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/'+name
size=454578705
def part(i):
    a=i*16*1024**2;b=min(a+16*1024**2,size)-1;p=out/f'download.part{i}'
    if p.exists() and p.stat().st_size==b-a+1:return p
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'})
            with urllib.request.urlopen(req,timeout=90) as r:
                assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/')
                data=r.read()
            assert len(data)==b-a+1
            p.write_bytes(data);print('Downloaded part',i,flush=True);return p
        except Exception:
            if attempt==3:raise
            time.sleep(2)
f=out/name
if not f.exists():
    with ThreadPoolExecutor(max_workers=6) as pool:parts=list(pool.map(part,range((size+16*1024**2-1)//(16*1024**2))))
    with f.open('wb') as w:
        for p in parts:w.write(p.read_bytes())
    assert f.stat().st_size==size
    for p in parts:p.unlink()

# 不利用EPB41L3表达重新选择细胞；沿用原作者通过质控的细胞和亚群标签。
meta=pd.read_csv(src/'GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz',sep='\t')
genes=['EPB41L3','SPP1','LST1','TYROBP','FCER1G','CD68','C1QA','C1QB','C1QC','APOE','TREM2','LPL','NLRP3','IL1B','FCN1','S100A8','S100A9','EPCAM','COL1A1','CD3D','MS4A1']
values={};rows=0;maxval=0
with gzip.open(f,'rt') as r:
    cells=[v.strip('"') for v in r.readline().split()]
    assert len(cells)==len(set(cells))==len(meta)==43817
    assert set(cells)==set(meta.CellName)
    for line in r:
        gene,data=line.strip().split(None,1);gene=gene.strip('"');rows+=1
        if gene in genes:
            assert gene not in values
            v=np.fromstring(data,sep=' ');assert len(v)==len(cells) and np.isfinite(v).all() and (v>=0).all()
            values[gene]=v;print('Extracted',gene,flush=True)
assert set(values)==set(genes)
meta=meta.set_index('CellName').loc[cells].reset_index()
joined=pd.concat([meta,pd.DataFrame(values)],axis=1)
joined.to_csv(out/'cells_fixed_markers.csv.gz',index=False)
assert ((joined.Sub_Cluster=='hM13_TAM-SPP1')&(joined.Tissue=='T')).sum()==544
info={'source_url':url,'download_bytes':size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'expression_rows':rows,'cells':len(cells),'barcode_match':'exact_one_to_one_CellName','markers':genes,'normalization':'Published processed values used as supplied; no additional log or count normalization','raw_UMI_used':False}
(out/'extraction_verification.json').write_text(json.dumps(info,indent=2))
print(info,flush=True)
