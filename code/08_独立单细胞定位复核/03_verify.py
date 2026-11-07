"""独立重算关键分母与患者等权统计；不把静态检查当作分析完成。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
p=Path(__file__).parent
d=pd.read_csv(p/'cells_fixed_markers.csv.gz')
z=d[(d.Tissue=='T')&(d.Sub_Cluster=='hM13_TAM-SPP1')]
a=z.groupby('Sample').EPB41L3.agg(cells='size',detected=lambda s:(s>0).sum(),mean_expression='mean')
a['detection']=a.detected/a.cells
r=pd.read_csv(p/'SPP1_patient_EPB41L3.csv');r=r[r.Tissue=='T'].set_index('Sample').loc[a.index]
assert len(d)==43817 and len(z)==544 and len(a)==10
assert np.all(a.cells==r.cells) and np.all(a.detected==r.detected)
assert np.allclose(a.mean_expression,r.mean_expression) and np.allclose(a.detection,r.detection)
b=a[a.cells>=20]
s=pd.read_csv(p/'patient_equal_localization.csv')
s=s[(s.min_cells==20)&(s.Sub_Cluster=='hM13_TAM-SPP1')&(s.gene=='EPB41L3')].iloc[0]
assert len(b)==3 and np.isclose(s.detection,b.detection.mean())
assert np.isclose(s.mean_expression,b.mean_expression.mean())
assert all((p/f'Figure_GSE146771_localization.{ext}').stat().st_size>10000 for ext in ['pdf','png'])
info={'status':'PASS','cells':len(d),'SPP1_tumor_cells':len(z),'SPP1_tumor_patients':len(a),'SPP1_tumor_patients_ge20':len(b),'gene_detection_patient_summary':'independently_recomputed','patient_equal_summary':'independently_recomputed','all_SPP1_tumor_detected':int((z.EPB41L3>0).sum()),'normal_SPP1_cells':int(((d.Tissue=='N')&(d.Sub_Cluster=='hM13_TAM-SPP1')).sum())}
(p/'verification.json').write_text(json.dumps(info,indent=2));print(info)
