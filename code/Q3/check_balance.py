import numpy as np
from pathlib import Path
p=Path(r'private/results/Q3/refined_kirchhoff/kirchhoff_N640_p1.5_tol1e-10_step120.npz')
z=np.load(p)
r=z['r']; t=z['t']; C=z['C']; R=.02; HM=8e-7
faces=.5*(r[:-1]+r[1:]); vol=.5*np.diff(np.r_[0,faces,R]**2)
import zipfile, xml.etree.ElementTree as ET
with zipfile.ZipFile(Path(r'E:/2026数模/CUMCM2026Problems/A题/附件/附件1.xlsx')) as f:
 ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 root=ET.fromstring(f.read('xl/worksheets/sheet1.xml'))
 rows=[]
 for row in root.findall('.//m:sheetData/m:row',ns):
  vals=[float(c.find('m:v',ns).text) for c in row.findall('m:c',ns)[:3]]
  if len(vals)==3: rows.append(vals)
 rows=np.asarray(rows)
M=C@vol; q=R*HM*(C[:,-1]-np.interp(t,rows[:,0],rows[:,2]))
print('N',len(r)-1,'states',len(t),'relative balance',(M[-1]-M[0]+np.trapezoid(q,t))/M[0],'M0,Mend',M[0],M[-1])
