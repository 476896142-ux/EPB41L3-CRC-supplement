"""从 GEO 补充文件下载压缩 CEL；完整解压校验后才替换目标文件。"""
import csv, gzip, re, urllib.request, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

base = Path('D:/Bioinfor/colorectal/GEO_ML_hub10/raw')
folder = base / 'GSE20916_CEL'
rows = list(csv.DictReader(open(base / 'GSE20916_metadata.csv', encoding='utf8')))
rows = [r for r in rows if '(micro)' not in r['title'] and r['tissue'] in ['adenocarcinoma', 'normal colon']]

def fetch(row):
    gsm = row['geo_accession']; dest = folder / (gsm + '.CEL.gz')
    if dest.exists():
        try:
            assert len(gzip.decompress(dest.read_bytes())) > 10000000
            return gsm
        except Exception:
            pass
    if (folder / (gsm + '.CEL.part')).exists():
        return gsm + ' (EBI in progress)'
    text = gzip.decompress((base / 'GSE20916' / (gsm + '.soft.gz')).read_bytes()).decode()
    url = re.search(r'!Sample_supplementary_file = (\S+\.CEL\.gz)', text).group(1).replace('ftp://', 'https://')
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            assert len(gzip.decompress(data)) > 10000000
            temporary = folder / (gsm + '.ncbi.part')
            temporary.write_bytes(data)
            temporary.replace(dest)
            return gsm
        except Exception as error:
            if attempt == 3:
                return gsm + ' failed: ' + str(error)
            time.sleep(3)

with ThreadPoolExecutor(max_workers=3) as pool:
    futures = [pool.submit(fetch, row) for row in reversed(rows)]
    for future in as_completed(futures):
        print(future.result(), flush=True)
