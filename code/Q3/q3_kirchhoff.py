"""Independent conservative radial MOL/BDF grid study for Q3."""
from pathlib import Path
import argparse, json, time
import numpy as np
import zipfile
import xml.etree.ElementTree as ET
from scipy.integrate import solve_ivp
from scipy.special import expi
from scipy.sparse import lil_matrix

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/Q3/refined_kirchhoff'
R, HH, HM = .02, 25., 8.e-7

def ambient():
    path=ROOT.parent/'CUMCM2026Problems/A题/附件/附件1.xlsx'
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        doc=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
    rows=[]
    for row in doc.findall('.//m:sheetData/m:row',ns)[1:]:
        rows.append([float(c.find('m:v',ns).text) for c in row.findall('m:c',ns)[:3]])
    return np.asarray(rows,float).T

def run(N=160, power=1., flux='kirchhoff', rtol=1e-8, max_step=300., export=False):
    ta, te, ca=ambient()
    x=np.linspace(0.,1.,N+1)
    # Power > 1 creates small cells at the surface while keeping an axis node.
    r=R*(1-(1-x)**power)
    faces=.5*(r[:-1]+r[1:]); dr=np.diff(r)
    edges=np.r_[0.,faces,R]; vol=.5*np.diff(edges**2)
    geom=faces/dr
    n=N+1
    def rhs(t,y):
        T=y[0::2]; C=y[1::2]; Cs=np.maximum(C,1e-12)
        rhocp=(650+128*Cs)*(1450+2736*Cs/(Cs+1))
        k=.21+.38*Cs/(Cs+1)
        tf=.5*(T[:-1]+T[1:])+273.15
        af=.0024*np.exp(-3850/tf)
        if flux=='kirchhoff':
            K=Cs*np.exp(-.45/Cs)+.45*expi(-.45/Cs)
            fc=geom*af*np.diff(K)
        else:
            D=.0024*np.exp(-.45/Cs-3850/(T+273.15))
            fc=geom*.5*(D[:-1]+D[1:])*np.diff(C)
        ft=geom*.5*(k[:-1]+k[1:])*np.diff(T)
        dt=np.diff(np.r_[0.,ft,R*HH*(np.interp(t,ta,te)-T[-1])])/vol/rhocp
        dc=np.diff(np.r_[0.,fc,R*HM*(np.interp(t,ta,ca)-C[-1])])/vol
        dy=np.empty_like(y);dy[0::2]=dt;dy[1::2]=dc
        return dy
    sparsity=lil_matrix((2*n,2*n),dtype=int)
    for i in range(n):
        sparsity[2*i:2*i+2,2*max(0,i-1):2*min(n,i+2)]=1
    def event(t,y): return np.max(y[1::2])-.15
    event.terminal=True;event.direction=-1
    y0=np.tile([28.,2.55],n)
    before=time.perf_counter()
    sol=solve_ivp(rhs,(0,259200),y0,method='BDF',rtol=rtol,
       atol=np.tile([rtol*.01,rtol*.0001],n),max_step=max_step,
       jac_sparsity=sparsity.tocsc(),events=event,dense_output=True)
    duration=time.perf_counter()-before
    assert sol.success,sol.message
    assert len(sol.t_events[0]),'drying event absent'
    end=float(sol.t_events[0][0])
    ts=np.r_[np.arange(21600,end,21600),end]
    ys=sol.sol(ts); table=np.stack([np.interp([0,.005,.01,.015,.02],r,v) for v in ys[1::2].T])
    result=dict(N=N,power=power,flux=flux,rtol=rtol,max_step=max_step,
      elapsed_s=duration,drying_s=end,drying_h=end/3600,
      nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,steps=len(sol.t),
      table_times_h=(ts/3600).tolist(),table_C=table.tolist(),
      min_C=float(sol.y[1::2].min()),max_C=float(sol.y[1::2].max()),
      max_radial_increase=float(np.diff(sol.y[1::2],axis=0).max()),
      last_surface_C=float(ys[-1,-1]),
      geometry_volume_error=float(abs(vol.sum()-R*R/2)))
    OUT.mkdir(parents=True,exist_ok=True)
    tag=f'{flux}_N{N}_p{power:g}_tol{rtol:g}_step{max_step:g}'
    (OUT/(tag+'.json')).write_text(json.dumps(result,indent=2),encoding='utf8')
    np.savez_compressed(OUT/(tag+'.npz'),r=r,t=ts,T=ys[0::2].T,C=ys[1::2].T,
      sample_t=sol.t,sample_center_C=sol.y[1],sample_surface_C=sol.y[-1])
    if export:
        # Continue after the event with a new IVP; never extrapolate beyond it.
        integer_end=float(np.floor(end)+1)
        tail=solve_ivp(rhs,(end,integer_end),sol.y_events[0][0],method='BDF',
          rtol=rtol,atol=np.tile([rtol*.01,rtol*.0001],n),
          jac_sparsity=sparsity.tocsc(),dense_output=True,max_step=.1)
        assert tail.success
        last=tail.y[:,-1]
        previous=sol.sol(integer_end-1)
        assert np.max(last[1::2])<.15
        assert np.max(previous[1::2])>=.15
        sample_times=np.r_[np.arange(60,integer_end,60),integer_end]
        sample_y=sol.sol(sample_times[:-1])
        sample_y=np.column_stack([sample_y,last])
        targets=np.linspace(0,R,21)
        sample_C=np.stack([np.interp(targets,r,v) for v in sample_y[1::2].T])
        sample_T=np.stack([np.interp(targets,r,v) for v in sample_y[0::2].T])
        table_times=np.r_[np.arange(21600,integer_end,21600),integer_end]
        table_y=np.column_stack([sol.sol(table_times[:-1]),last])
        np.savez_compressed(OUT/(tag+'_export.npz'),
          times=sample_times,radius=targets,C=sample_C,T=sample_T,
          table_times=table_times,table_C=np.stack([np.interp(targets[::5],r,v) for v in table_y[1::2].T]),
          integer_end=integer_end,event_time=end,previous_max_C=np.max(previous[1::2]),
          terminal_max_C=np.max(last[1::2]))
        result.update(integer_end_s=integer_end,previous_max_C=float(np.max(previous[1::2])),
           terminal_max_C=float(np.max(last[1::2])))
        (OUT/(tag+'.json')).write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=160)
    p.add_argument('--power',type=float,default=1.)
    p.add_argument('--flux',default='kirchhoff',choices=['kirchhoff','arithmetic'])
    p.add_argument('--rtol',type=float,default=1e-8)
    p.add_argument('--max-step',type=float,default=300.)
    p.add_argument('--export',action='store_true')
    a=p.parse_args();run(a.N,a.power,a.flux,a.rtol,a.max_step,a.export)
