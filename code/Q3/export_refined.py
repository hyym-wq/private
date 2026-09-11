from pathlib import Path
import numpy as np, json
from zipfile import ZipFile

root=Path(__file__).resolve().parents[2]
out=root/'results/Q3/refined_kirchhoff'
p=out/'kirchhoff_N1280_p1_tol1e-10_step120_export.npz'
z=np.load(p)
print(z.files, {k:z[k].shape for k in z.files})
r=np.linspace(0,.02,21); ts=z['times']; C=z['C']; T=z['T']
assert C.shape[0]==ts.size and np.all(np.diff(ts)>0)
dr=r[1]-r[0]; j=np.rint(np.array([0,.005,.01,.015,.02])/dr).astype(int)
target=np.r_[np.arange(21600., ts[-1],21600.), ts[-1]]
rows=[]
for t in target:
    i=int(np.argmin(abs(ts-t))); assert abs(ts[i]-t)<1e-7
    rows.append([float(ts[i]),*map(float,C[i,j])])
head='时间/h,0 cm,0.5 cm,1 cm,1.5 cm,2 cm'
lines=[head]+[f'{row[0]/3600:.6f},'+','.join(f'{v:.8f}' for v in row[1:]) for row in rows]
(out/'table5_refined.csv').write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
def col(n):
    s=''; n+=1
    while n: n,rem=divmod(n-1,26); s=chr(65+rem)+s
    return s
def cell(ref,val):
    if isinstance(val,str): return f'<c r="{ref}" t="inlineStr"><is><t>{val}</t></is></c>'
    return f'<c r="{ref}"><v>{val:.12g}</v></c>'
header=['时间\\到药材中心的距离']+[f'{x:.1f}' for x in np.arange(0,2.0001,.1)]
xml_rows=['<row r="1">'+''.join(cell(col(i)+'1',v) for i,v in enumerate(header))+'</row>']
for ii,(t,crow) in enumerate(zip(ts,C),start=2):
    vals=[float(t)]+[float(crow[k]) for k in range(21)]
    xml_rows.append(f'<row r="{ii}">'+''.join(cell(col(i)+str(ii),v) for i,v in enumerate(vals))+'</row>')
sheet='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'+''.join(xml_rows)+'</sheetData></worksheet>'
ct='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>'
rels='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'
wbxml='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
wbrels='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
path=out/'result3_refined.xlsx'
with ZipFile(path,'w') as f:
    f.writestr('[Content_Types].xml',ct); f.writestr('_rels/.rels',rels); f.writestr('xl/workbook.xml',wbxml); f.writestr('xl/_rels/workbook.xml.rels',wbrels); f.writestr('xl/worksheets/sheet1.xml',sheet)
print('written',path,'frames',len(ts))
for x in rows: print(x)
summary={'method':'Q2 PDE with Kirchhoff moisture flux + BDF','N_internal':1280,'dt_control':'adaptive BDF, max_step=120 s','event_s':205809.80110531385,'first_integer_dry_s':205810.,'first_integer_dry_h':205810/3600,'frames':len(ts),'output_spacing_s':60,'output_spacing_cm':.1,'table':rows,'mesh_times_h':{'320':57.16955889,'640':57.16942356,'1280':57.1693891959},'surface_refined_N640_p1.5_s':205809.81887188723,'convergence_last_s':0.01776657338,'checks':{'finite':bool(np.isfinite(C).all() and np.isfinite(T).all()),'center_max':bool(np.max(np.diff(C,axis=1))<1e-10)}}
(out/'refined_final_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
