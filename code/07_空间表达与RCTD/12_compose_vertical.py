"""将原始矢量 PDF 上下拼接，保留面板内部尺度及统计值。"""
from pathlib import Path
from io import BytesIO
import json
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
import pypdfium2 as pdfium

root=Path('D:/Bioinfor/colorectal/GSE225857_spatial')
out=root/'combined_figures';out.mkdir(exist_ok=True)
specs=[('Figure_S1_combined',root/'figures/Figure_spatial_expression.pdf',out/'S_bottom_panels.pdf','A  Spatial expression of EPB41L3, KLB and ANGPTL1'),
       ('Figure_R1_combined',root/'rctd_results/Figure_RCTD_spatial.pdf',out/'R_bottom_panels.pdf','A  RCTD-derived cell-type distributions')]
checks=[]
for name,top_path,bottom_path,title in specs:
    top=PdfReader(top_path).pages[0];bottom=PdfReader(bottom_path).pages[0]
    width=864;header=28;gap=14
    s_top=width/float(top.mediabox.width);s_bottom=width/float(bottom.mediabox.width)
    h_top=float(top.mediabox.height)*s_top;h_bottom=float(bottom.mediabox.height)*s_bottom
    height=header+h_top+gap+h_bottom
    writer=PdfWriter();page=writer.add_blank_page(width=width,height=height)
    page.merge_transformed_page(top,Transformation().scale(s_top).translate(0,h_bottom+gap))
    page.merge_transformed_page(bottom,Transformation().scale(s_bottom))
    overlay=BytesIO();c=canvas.Canvas(overlay,pagesize=(width,height));c.setFont('Helvetica-Bold',13);c.drawString(12,height-19,title);c.save();overlay.seek(0)
    page.merge_page(PdfReader(overlay).pages[0]);dest=out/(name+'.pdf')
    with dest.open('wb') as f:writer.write(f)
    doc=pdfium.PdfDocument(str(dest));im=doc[0].render(scale=300/72).to_pil();im.save(out/(name+'.png'),dpi=(300,300))
    checks.append(dict(file=name,pages=len(doc),width=im.width,height=im.height,dpi=300,top_source=str(top_path),bottom_source=str(bottom_path)))
    im.thumbnail((1100,1700));im.save(out/('preview_'+name+'.png'))
(out/'figure_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
