"""核对原始注释并提取三基因原始 UMI，供 R 与缓存逐值比对。"""
from pathlib import Path
import csv
import gzip
import hashlib

root = Path("D:/Bioinfor/colorectal/single_raw")
out = root / "three_gene_analysis"
genes = {"EPB41L3", "KLB", "ANGPTL1"}
records = []
for accession in ["GSE132465", "GSE132257"]:
    annotation = next((root / "raw" / accession).glob("*annotation*.gz"))
    with gzip.open(annotation, "rt") as stream:
        raw = list(csv.DictReader(stream, delimiter="\t"))
    with open(root / "results" / f"annotation_{accession}.csv", encoding="utf-8-sig") as stream:
        cached = list(csv.DictReader(stream))
    assert len(raw) == len(cached)
    for left, right in zip(raw, cached):
        assert left["Index"] == right[""]
        assert all(left[key] == right[key] for key in left if key != "Index")
    matrix = next((root / "raw" / accession).glob("*UMI*.gz"))
    selected = {}
    with gzip.open(matrix, "rt") as stream:
        header = next(stream).strip().replace('"', '').split("\t")
        for line in stream:
            name = line.split("\t", 1)[0].strip('"')
            if name in genes:
                selected[name] = line.rstrip().replace('"', '').split("\t")[1:]
                if len(selected) == 3:
                    break
    assert set(selected) == genes
    n = len(selected["KLB"])
    cells = header[-n:]
    with open(out / f"{accession}_raw_three_genes.csv", "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cell", "EPB41L3", "KLB", "ANGPTL1"])
        writer.writerows(zip(cells, selected["EPB41L3"], selected["KLB"], selected["ANGPTL1"]))
    records.append([accession, len(raw), "exact match", matrix.name, hashlib.sha256(matrix.read_bytes()).hexdigest()])
    print(accession, len(raw), "annotation verified; raw targets extracted", flush=True)
with open(out / "input_verification.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["dataset", "annotated_cells", "annotation_check", "raw_matrix", "sha256"])
    writer.writerows(records)
