from pathlib import Path
import sqlite3, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# 只读取已完成的统计结果；不重做 MR，不修改原始图。
root = Path('D:/Bioinfor/colorectal/MR')
out = root / '组学MR论文写作'
figdir = out / '附图'
figdir.mkdir(exist_ok=True)
deg = pd.read_csv('D:/Bioinfor/colorectal/GEO_CRC_raw/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv')
deg.GeneSymbol = deg.GeneSymbol.str.strip().str.upper()
deg = deg.drop_duplicates('GeneSymbol')
eq = pd.read_csv(root / 'eQTL/Mr_colorectal_result/eqtl_deg_intersection.csv')
pq = pd.read_csv(root / 'pQTL/MR_IVW_FDR05_DEG_intersection.csv')
scan = pd.read_csv(root / 'eQTL/Mr_colorectal_result/mr_res.tsv', sep='\t')
scan = scan[(scan.method == 'Inverse variance weighted') & (scan.pval < .05)].copy()
scan['ens'] = scan.exposure.str.split('.').str[0]
# 用本地同一注释包恢复完整 MR 集合，避免以已取交集的集合绘图。
con = sqlite3.connect('file:D:/R/R-4.6.1/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite?mode=ro', uri=True)
mapping = pd.read_sql_query('SELECT e.ensembl_id, g.symbol FROM ensembl e JOIN gene_info g ON e._id=g._id ORDER BY e.rowid,g.rowid', con)
mapping = mapping.drop_duplicates('ensembl_id').set_index('ensembl_id').symbol
scan['symbol'] = scan.ens.map(mapping)
eset = set(scan.symbol.fillna(scan.ens).str.upper())
pscan = pd.read_csv(root / 'pQTL/MR_IVW_results.csv')
pset = set(pscan.loc[pscan.pval < .05, 'exposure'].str.upper())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
audit=[]
for kind, data, mrset, color in [('eQTL',eq,eset,'#287C8E'),('pQTL',pq,pset,'#AD5C3E')]:
    # 直接按具名字段排序，避免标签与数值错配。
    data=data.sort_values('pval').reset_index(drop=True)
    fig=plt.figure(figsize=(12,10.5),facecolor='white')
    gs=fig.add_gridspec(2,2,height_ratios=[1,2],hspace=.32,wspace=.28)
    for j, cat in enumerate(['Up','Down']):
        ax=fig.add_subplot(gs[0,j]); ds=set(deg.loc[deg.Category==cat,'GeneSymbol'].dropna())
        both=ds & mrset
        expected=set(data.loc[data.DEG_Category==cat,'GeneSymbol'])
        assert both==expected, (kind,cat,both ^ expected)
        for center,fill in [((.40,.48),'#7299BB'),((.65,.48),color)]:
            ax.add_patch(Circle(center,.25,facecolor=fill,edgecolor=fill,alpha=.35,lw=1.5))
        for x,val in [(.32,len(ds-mrset)),(.525,len(both)),(.73,len(mrset-ds))]:
            ax.text(x,.48,str(val),ha='center',va='center',fontsize=18,fontweight='bold')
        ax.text(.23,.80,cat+'regulated\nDEGs (n = '+str(len(ds))+')',ha='center',fontsize=10)
        ax.text(.82,.80,kind+' IVW P < 0.05\n(n = '+str(len(mrset))+')',ha='center',fontsize=10)
        ax.text(-.06,1.07,chr(65+j),weight='bold',fontsize=17,va='top')
        ax.set(xlim=(0,1.05),ylim=(.1,1));ax.set_aspect('equal');ax.axis('off')
        audit.append(dict(omics=kind,DEG_group=cat,DEG=len(ds),MR=len(mrset),overlap=len(both),genes=';'.join(sorted(both))))
    # 森林图与数值表共用行坐标；横轴为 OR 的对数尺度。
    sub=gs[1,:].subgridspec(1,3,width_ratios=[1.6,3.2,3.1],wspace=.04)
    left,ax,right=[fig.add_subplot(sub[0,i]) for i in range(3)]
    n=len(data);y=np.arange(n)[::-1]
    for a in (left,ax,right): a.set_ylim(-.8,n+.55)
    for row in range(n):
        for a in (left,ax,right):
            if row%2==0:a.axhspan(row-.45,row+.45,color='#F1F4F6',zorder=0)
    left.set_xlim(0,1);right.set_xlim(0,1)
    left.text(0,n+.35,'C',fontsize=17,weight='bold')
    left.text(.0,n-.05,'Gene',weight='bold');left.text(.83,n-.05,'SNPs',ha='center',weight='bold')
    right.text(.02,n-.05,'OR (95% CI)',weight='bold');right.text(.98,n-.05,'P',ha='right',weight='bold')
    for i,row in data.iterrows():
        yy=y[i]
        left.text(0,yy,row.GeneSymbol,va='center',weight='bold' if row.GeneSymbol=='LAMC2' else 'normal')
        left.text(.83,yy,str(int(row.get('nsnp',row.get('nsnp.x')))),ha='center',va='center')
        right.text(.02,yy,f'{row.OR:.3f} ({row.OR_low:.3f}–{row.OR_high:.3f})',va='center',fontsize=10)
        right.text(.98,yy,f'{row.pval:.3g}',ha='right',va='center',fontsize=10)
    ax.errorbar(data.OR,y,xerr=np.array([data.OR-data.OR_low,data.OR_high-data.OR]),fmt='s',color=color,ecolor=color,capsize=3,ms=6,lw=1.5)
    ax.axvline(1,color='#5F6770',ls='--',lw=1)
    ax.set_xscale('log');ax.set_xlim(.25,2.6)
    ax.set_xticks([.25,.5,1,2]);ax.set_xticklabels(['0.25','0.5','1','2'])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_yticks([]);ax.set_xlabel('Odds ratio (95% CI), log scale')
    ax.spines['left'].set_visible(False)
    left.axis('off');right.axis('off')
    fig.suptitle(kind+' associations with colorectal cancer and differential expression',fontsize=16,fontweight='bold',y=.97)
    fig.text(.5,.025,'Candidate selection: unadjusted IVW P < 0.05  |  Venn circle areas are schematic',ha='center',fontsize=10,color='#555555')
    fig.subplots_adjust(left=.055,right=.97,bottom=.09,top=.91)
    for ext in ['png','pdf']:fig.savefig(figdir/f'Figure_{kind}_Venn_forest.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
pd.DataFrame(audit).to_csv(figdir/'韦恩图集合核对.csv',index=False,encoding='utf-8-sig')
print(json.dumps(audit,ensure_ascii=False,indent=2))

# 在结论段落后直接嵌图，同时保留无图版本。
path=out/'方法与结论_中英文.md'
text=path.read_text(encoding='utf-8-sig')
backup=out/'方法与结论_中英文_插图前.md'
if not backup.exists():backup.write_text(text,encoding='utf-8-sig')
zh='''

eQTL 分析中，7 个候选来自上调差异基因，5 个来自下调差异基因（图 1A、B）；对应 MR 效应见图 1C。pQTL 分析分别获得 2 个上调和 6 个下调差异基因候选（图 2A、B），对应 MR 效应见图 2C。

![图1 eQTL韦恩图与森林图](D:/Bioinfor/colorectal/MR/组学MR论文写作/附图/Figure_eQTL_Venn_forest.png)

**图 1. eQTL 关联与差异表达整合筛选。**（A、B）IVW 原始 P < 0.05 的基因与上调、下调差异基因的交集。（C）12 个交集基因的 IVW OR 及 95% CI；虚线表示 OR = 1。P 为未经多重检验校正的数值。韦恩图圆面积不代表集合大小。

![图2 pQTL韦恩图与森林图](D:/Bioinfor/colorectal/MR/组学MR论文写作/附图/Figure_pQTL_Venn_forest.png)

**图 2. pQTL 关联与差异表达整合筛选。**（A、B）IVW 原始 P < 0.05 的蛋白对应基因与上调、下调差异基因的交集。（C）8 个交集候选的 IVW OR 及 95% CI；虚线表示 OR = 1。P 为未经校正的数值；GREM1 和 HHIP 通过全扫描 FDR 校正。韦恩图圆面积不代表集合大小。

'''
assert sum((eq.DEG_Category=='Up'))==7 and sum(pq.DEG_Category=='Up')==2
en='''

The eQTL analysis identified seven candidates among upregulated genes and five among downregulated genes (Fig. 1A,B); their MR estimates are shown in Fig. 1C. The pQTL analysis identified two and six candidates, respectively (Fig. 2A,B), with MR estimates shown in Fig. 2C.

**Figure 1. Integration of eQTL associations with differential expression.** (A,B) Overlap between genes with an unadjusted IVW P < 0.05 and upregulated or downregulated genes. (C) IVW ORs and 95% CIs for the 12 overlapping genes. The dashed line indicates OR = 1. P values are unadjusted. Venn circle areas are schematic.

**Figure 2. Integration of pQTL associations with differential expression.** (A,B) Overlap between genes encoding proteins with an unadjusted IVW P < 0.05 and upregulated or downregulated genes. (C) IVW ORs and 95% CIs for the eight overlapping candidates. The dashed line indicates OR = 1. P values are unadjusted; GREM1 and HHIP passed the full-scan FDR threshold. Venn circle areas are schematic.

'''
if '![图1' not in text:
    text=text.replace('## Methods (English)',zh+'## Methods (English)')
    text=text.replace('## 已核对的数据来源',en+'## 已核对的数据来源')
    path.write_text(text,encoding='utf-8-sig')
