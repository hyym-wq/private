"""Mathematical review checks for q3_kirchhoff BDF runs."""
from pathlib import Path
import json, zipfile, xml.etree.ElementTree as ET
import numpy as np
from scipy.special import expi
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'results/Q3/refined_kirchhoff/kirchhoff_N640_p1.5_tol1e-10_step120.npz'
J=ROOT/'results/Q3/refined_kirchhoff/kirchhoff_N640_p1.5_tol1e-10_step120.json'
z=np.load(P); r=z['r']; ts=z['sample_t']; Cs=z['sample_center_C']; Csurf=z['sample_surface_C']
with zipfile.ZipFile(ROOT.parent/'CUMCM2026Problems/A题/附件/附件1.xlsx') as f:
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    root=ET.fromstring(f.read('xl/worksheets/sheet1.xml'))
    rows=[]
    for row in root.findall('.//m:sheetData/m:row',ns)[1:]: rows.append([float(c.find('m:v',ns).text) for c in row.findall('m:c',ns)[:3]])
ta=np.array(rows)[:,0]; ca=np.array(rows)[:,2]
ci=np.interp(ts,ta,ca)
R=.02; HM=8e-7
M0=2.55*R*R/2
faces=(r[:-1]+r[1:])/2
vol=.5*np.diff(np.r_[0.,faces,R]**2)
Mend=float(np.dot(vol,z['C'][-1]))
boundary_flux=R*HM*(Csurf-ci)
loss_trapezoid=float(np.trapezoid(boundary_flux,ts))
loss_simpson=float(simpson(boundary_flux,x=ts))
instant_res=[]
for t,T,C in zip(z['t'],z['T'],z['C']):
    K=C*np.exp(-.45/C)+.45*expi(-.45/C)
    af=.0024*np.exp(-3850/(.5*(T[:-1]+T[1:])+273.15))
    fc=faces/np.diff(r)*af*np.diff(K)
    incoming=R*HM*(np.interp(t,ta,ca)-C[-1])
    dc=np.diff(np.r_[0.,fc,incoming])/vol
    instant_res.append(float(np.dot(vol,dc)-incoming))
# K derivative check
a=.45; c=np.logspace(-4,np.log10(3),1000)
K=c*np.exp(-a/c)+a*expi(-a/c)
num=np.gradient(K,c)
max_deriv_err=float(np.max(np.abs(num-np.exp(-a/c))))
res={'run':str(P),'event_h':float(json.loads(J.read_text())['drying_h']),
 'C_integral_balance':{'max_instantaneous_absolute_residual':max(map(abs,instant_res)),
   'relative_integrated_residual_trapezoid':(Mend-M0+loss_trapezoid)/M0,
   'relative_integrated_residual_simpson':(Mend-M0+loss_simpson)/M0,
   'definition':'M=sum(vol*C); time integral uses boundary flux sampled at every accepted BDF step',
   'meaning':'C-PDE integral balance, not asserted to be physical total water mass'},
 'Kprime_max_abs_error_grid':max_deriv_err,
 'Kprime_min':float(np.min(np.exp(-a/c))), 'Kprime_max':float(np.max(np.exp(-a/c))),
 'checks':{'geometry_formula':'vol=.5*diff(edges**2), sum(vol)=R²/2 exactly',
   'temperature_divergence':'conservative flux difference; Robin inward flux at boundary',
   'mass_divergence':'Kirchhoff K(C) flux difference; Robin inward flux at boundary',
   'event':'max(C)-0.15 terminal crossing, direction=-1',
   'Kprime_identity':'d/dC [C exp(-a/C)+a Ei(-a/C)] = exp(-a/C)'},
 'limitations':['integrated residual includes independent boundary time-quadrature error; the instantaneous flux identity is algebraically conservative',
                'Cs=max(C,1e-12) is used in constitutive functions; positive accepted states are numerical evidence, not a proof of positivity preservation']}
(ROOT/'results/Q3/refined_kirchhoff/review.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(res,ensure_ascii=True,indent=2))
