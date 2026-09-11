"""Q4 shrinking-radius ALE radial model.

The dry-basis moisture and temperature are attached to material shells.  With
xi=r/R(t), the shrinkage velocity cancels from the material derivative and
the diffusion operator acquires the factor R(t)^-2.  Robin fluxes are kept in
physical units, so the boundary terms scale as 1/R(t).
"""
from pathlib import Path
import argparse, json, zipfile, xml.etree.ElementTree as ET, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expi
from scipy.sparse import lil_matrix

ROOT=Path(__file__).resolve().parents[2]
ATT=ROOT.parent/'CUMCM2026Problems/A题/附件'
OUT=ROOT/'results/Q4'
HH,HM=25.,8.e-7

def read_xlsx3(path):
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    vals=[]
    for row in root.findall('.//m:sheetData/m:row',ns)[1:]:
        rr=[]
        for c in row.findall('m:c',ns)[:3]:
            v=c.find('m:v',ns); rr.append(float(v.text) if v is not None else np.nan)
        if len(rr)>=3: vals.append(rr)
    return np.asarray(vals,float).T

def ambient():
    a=read_xlsx3(ATT/'附件1.xlsx'); return a[0],a[1],a[2]

def radius_data():
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(ATT/'附件2.xlsx') as z:
        root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    vals=[]
    for row in root.findall('.//m:sheetData/m:row',ns)[1:]:
        cs=row.findall('m:c',ns); vv=[]
        for c in cs[:2]:
            v=c.find('m:v',ns); vv.append(float(v.text) if v is not None else np.nan)
        if len(vv)==2: vals.append(vv)
    a=np.asarray(vals,float); return a[:,0],a[:,1]*1e-2

def run(N=320, rtol=1e-8, max_step=300., tmax=400000., export=True):
    ta,te,ca=ambient(); tr,rr=radius_data(); R0=rr[0]
    xi=np.linspace(0.,1.,N+1); face=.5*(xi[:-1]+xi[1:]); dxi=np.diff(xi)
    edges=np.r_[0.,face,1.]; vol=.5*np.diff(edges**2); geom=face/dxi; n=N+1
    def Rfun(t): return float(np.interp(t,tr,rr,left=rr[0],right=rr[-1]))
    def rhs(t,y):
        T=y[0::2]; C=np.maximum(y[1::2],1e-14); Rt=Rfun(t)
        rhocp=(760+90*C)*(1850+2150*C/(C+1))
        k=.12+.20*C/(C+1)
        tf=.5*(T[:-1]+T[1:])+273.15
        af=4.2e-4*np.exp(-3850/tf)
        K=C*np.exp(-.30/C)+.30*expi(-.30/C)
        fc=geom*af*np.diff(K)
        ft=geom*.5*(k[:-1]+k[1:])*np.diff(T)
        Tin=np.interp(t,ta,te); Cin=np.interp(t,ta,ca)
        # face fluxes are r*coefficient*gradient; Robin value at physical R
        ft=np.r_[0.,ft,Rt*HH*(Tin-T[-1])]
        fc=np.r_[0.,fc,Rt*HM*(Cin-C[-1])]
        dt=np.diff(ft)/vol/(Rt*Rt)/rhocp
        dc=np.diff(fc)/vol/(Rt*Rt)
        out=np.empty_like(y); out[0::2]=dt; out[1::2]=dc; return out
    sp=lil_matrix((2*n,2*n),dtype=int)
    for i in range(n): sp[2*i:2*i+2,2*max(i-1,0):2*min(i+2,n)]=1
    def event(t,y): return np.max(y[1::2])-.15
    event.terminal=True; event.direction=-1
    y0=np.tile([28.,2.55],n); t0=time.perf_counter()
    sol=solve_ivp(rhs,(0,tmax),y0,method='BDF',rtol=rtol,
      atol=np.tile([rtol*.01,rtol*.0001],n),max_step=max_step,
      jac_sparsity=sp.tocsc(),events=event,dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    if not len(sol.t_events[0]): raise RuntimeError('drying event absent')
    end=float(sol.t_events[0][0]); elapsed=time.perf_counter()-t0
    # six-hour table through the event; material coordinates are converted to physical distance
    ts=np.r_[np.arange(21600,end,21600),end]; ys=sol.sol(ts)
    targets=np.arange(0,.02001,.005)
    tab=[]
    for j,t in enumerate(ts):
        Rt=Rfun(t); C=ys[1::2,j]; row=[]
        for q in targets: row.append(float(np.interp(q/Rt,xi,C)) if q<=Rt+1e-14 else np.nan)
        row.append(float(C[-1])); tab.append(row)
    result={'N':N,'rtol':rtol,'max_step':max_step,'drying_s':end,'drying_h':end/3600,
      'elapsed_s':elapsed,'steps':len(sol.t),'nfev':sol.nfev,'njev':sol.njev,'nlu':sol.nlu,
      'radius_final_cm':Rfun(end)*100,'table_times_h':(ts/3600).tolist(),'table_C':tab,
      'min_C':float(sol.y[1::2].min()),'max_C':float(sol.y[1::2].max()),
      'radius_data_end_s':float(tr[-1]),'model':'material-ALE xi=r/R(t), Kirchhoff moisture flux'}
    OUT.mkdir(parents=True,exist_ok=True); tag=f'q4_ALE_N{N}_tol{rtol:g}'
    (OUT/(tag+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    if export:
      integer_end=np.floor(end)+1.; tail=solve_ivp(rhs,(end,integer_end),sol.y_events[0][0],method='BDF',rtol=rtol,
        atol=np.tile([rtol*.01,rtol*.0001],n),max_step=.5,dense_output=True,jac_sparsity=sp.tocsc())
      times=np.r_[np.arange(60,integer_end,60),integer_end]; Y=sol.sol(times[:-1]); Y=np.column_stack([Y,tail.y[:,-1]])
      # save xi grid and physical-radius vector per time; callers can remap to 0.1 cm
      np.savez_compressed(OUT/(tag+'_full.npz'),times=times,xi=xi,R=np.array([Rfun(t) for t in times]),T=Y[0::2].T,C=Y[1::2].T)
      result.update(integer_end_s=float(integer_end),previous_max_C=float(np.max(sol.sol(integer_end-1)[1::2])),terminal_max_C=float(np.max(tail.y[1::2,-1])))
      (OUT/(tag+'.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(result,ensure_ascii=False),flush=True); return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=320); p.add_argument('--rtol',type=float,default=1e-8); p.add_argument('--max-step',type=float,default=300.); p.add_argument('--no-export',action='store_true'); a=p.parse_args(); run(a.N,a.rtol,a.max_step,export=not a.no_export)
