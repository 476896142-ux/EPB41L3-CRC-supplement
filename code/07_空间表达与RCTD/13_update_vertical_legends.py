from pathlib import Path
import shutil
root=Path('D:/Bioinfor/colorectal/GSE225857_spatial')
out=root/'combined_figures'
s='''**Figure S1. Spatial expression, detection and autocorrelation of candidate genes in colorectal tumors.** (A) Spatial expression of EPB41L3, KLB and ANGPTL1. Rows represent sections C1–C4; columns show the original H&E image and log-normalized gene expression. Gray spots indicate zero counts. Each gene uses a shared color scale across sections. Expression-panel H&E backgrounds were uniformly blended with white to 30% image intensity. Scale bars, 500 μm, calculated from deposited scaling factors and the standard 55-μm Visium spot diameter. Retained spot counts were 2,054, 3,077, 3,330 and 3,969, respectively. (B) Dot size and color indicate the percentage of spots with detectable expression. (C) Moran’s I for each gene and section; asterisks indicate Benjamini–Hochberg-adjusted permutation P < 0.05 across 12 tests, with 999 permutations per test.

**中文图注：** 图 S1｜候选基因的空间表达、检出率与空间自相关。A，C1–C4 切片的 H&E 及 EPB41L3、KLB、ANGPTL1 表达分布；灰点表示零计数，同一基因跨切片统一色标，比例尺为500 μm。B，点大小与颜色表示表达检出率。C，各基因的 Moran’s I；星号表示999次置换检验后、12项检验的 BH 校正 P < 0.05。
'''
r='''**Figure R1. RCTD-derived cellular distributions and candidate-gene associations in colorectal tumors.** (A) Spatial distributions of myeloid, stromal and epithelial relative weights across sections C1–C4. Full-mode RCTD weights were normalized within each spot; myeloid and stromal subtypes were summed according to the single-cell reference annotations. Maps use a fixed 0–1 scale and H&E backgrounds displayed at 30% image intensity. Scale bars, 500 μm. (B) Total-count-adjusted rank correlations between candidate-gene expression and relative cell-type weights, calculated separately for each section. (C) Spearman correlations between relative weights and corresponding marker-expression scores. Numbers in B and C indicate correlation coefficients; color scales range from −1 to 1 and 0 to 1, respectively. Marker scores are the mean log-normalized expression of EPCAM/KRT8/KRT18 (epithelial), COL1A1/COL1A2/DCN/LUM (stromal), LST1/TYROBP/FCER1G/CTSS (myeloid), and ACTA2/TAGLN/CNN1 (smooth muscle).

**中文图注：** 图 R1｜RCTD 细胞组成空间分布及候选基因关联。A，C1–C4 切片中髓系、基质及上皮细胞的相对权重；权重在每个空间点内归一化，色标统一为0–1，比例尺为500 μm。B，三靶基因表达与细胞相对权重的总计数校正秩相关。C，细胞相对权重与相应标志基因平均表达的 Spearman 相关。B、C 数字为相关系数，色标分别为−1至1和0至1。
'''
(out/'Figure_legends.md').write_text('# Updated figure legends\n\n'+s+'\n'+r,encoding='utf-8')
for rel in ['Methods_Results_figure_legends.md','rctd_results/RCTD_results_and_legend.md']:
    p=root/rel;backup=root/'archive'/'before_vertical_merge'/rel
    backup.parent.mkdir(parents=True,exist_ok=True)
    if not backup.exists():shutil.copy2(p,backup)
    t=p.read_text(encoding='utf-8')
    t=t.replace('(Fig. S1)','(Fig. S1A)').replace('(Fig. S2)','(Fig. S1B,C)').replace('(Fig. R1)','(Fig. R1A)').replace('(Fig. R2A)','(Fig. R1B)').replace('(Fig. R2B)','(Fig. R1C)')
    if rel.startswith('rctd_results'):
        a=t.index('## Figure legends');b=t.index('## 复现',a)
        t=t[:a]+'## Figure legends\n\n'+r+'\n'+t[b:]
        t=t.replace('## 图与数据','## 图与数据\n\n**当前论文使用：[Figure R1 上下组合图](../combined_figures/Figure_R1_combined.pdf)，A为空间分布、B为基因关联、C为标志一致性。以下独立图及此前组合仅供单独使用。**')
    else:
        t=t[:t.index('## Figure legends')]+'## Figure legends\n\n'+s+'\n'+r+'\n\n图文件：[Figure S1](combined_figures/Figure_S1_combined.pdf)；[Figure R1](combined_figures/Figure_R1_combined.pdf)。\n'
    p.write_text(t,encoding='utf-8')
for rel in ['README.md','RESULTS_README.md']:
    p=root/rel;t=p.read_text(encoding='utf-8');first,rest=t.split('\n',1)
    t=first+'\n\n**当前上下组图：** [Figure S1（A–C）](combined_figures/Figure_S1_combined.pdf)；[Figure R1（A–C）](combined_figures/Figure_R1_combined.pdf)。[更新图注](combined_figures/Figure_legends.md)。独立图及此前合并图保留，全文采用本次编号。\n'+rest
    p.write_text(t,encoding='utf-8')
(out/'README.md').write_text('''# 当前论文组图

- Figure_S1_combined.pdf/png：A为原S1；B、C分别为原S2的检出率及Moran自相关。
- Figure_R1_combined.pdf/png：A为原R1；B、C分别为原R2合并图的基因关联及标志一致性。
- Figure_legends.md：新版中英文图注。
- PDF保留原图矢量元素；PNG为300 dpi。原始结果及独立图保留。

复现顺序：scripts/11_vertical_panels.R → scripts/12_compose_vertical.py。文稿更新脚本为 scripts/13_update_vertical_legends.py（一次性文稿迁移，不必在重绘时重复运行）。
''',encoding='utf-8')
print('Updated two figure legends and manuscript panel citations.')
