"""Memory-light Q3 solver for very fine radial grids.

It uses exactly the Q2 operators, but stores only 60-second frames and stops
as soon as the maximum nodal moisture is below 0.15. This permits N>=320
without allocating a 72-hour field history.
"""
from pathlib import Path
import argparse, importlib.util, json, platform, time, hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "数学建模/数学建模/code/p2_model.py"
ATT = ROOT.parent / "CUMCM2026Problems/A题/附件/附件1.xlsx"
OUT = ROOT / "results/Q3"

def run(N, dt, max_iter=3):
    spec=importlib.util.spec_from_file_location('q2',MODEL)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    Tinf,Cinf=m.make_ambient(ATT)
    R=m.R_CYL; dr=R/N; nmax=int(round(259200/dt)); stride=int(round(60/dt))
    T=np.full(N+1,m.T0); C=np.full(N+1,m.C0)
    times=[0.0]; Cs=[C.copy()]; Ts=[T.copy()]; maxima=[float(C.max())]
    mean_iter=0.; start=time.perf_counter(); stop=None
    zero=np.zeros(N+1)
    for n in range(1,nmax+1):
        tn=n*dt; to=(n-1)*dt
        Ti_n,Ti_o=Tinf(tn),Tinf(to); Ci_n,Ci_o=Cinf(tn),Cinf(to)
        Told,Cold=T.copy(),C.copy(); Cstar,Tstar=C.copy(),T.copy()
        aho,bho,cho,gto=m.build_op_heat(dr,N,Cold,R)
        sto=zero.copy(); sto[N]=gto*Ti_o
        used=max_iter
        for k in range(max_iter):
            ahn,bhn,chn,gtn=m.build_op_heat(dr,N,Cstar,R)
            stn=zero.copy(); stn[N]=gtn*Ti_n
            Tnew=m._step(Told,(ahn,bhn,chn),(aho,bho,cho),sto,stn,dt,.5)
            Cnew=m._mass_step(Cold,Cstar,Told,Tnew,Ci_o,Ci_n,dt,.5,dr,N,R,True)
            dT=np.max(np.abs(Tnew-Tstar)); dC=np.max(np.abs(Cnew-Cstar))
            Tstar,Cstar=Tnew,Cnew
            if dT<1e-9 and dC<1e-9:
                used=k+1; break
        T,C=Tnew,Cnew; mean_iter+=used
        mm=float(C.max()); maxima.append(mm)
        if n%stride==0 or mm<.15:
            times.append(tn); Cs.append(C.copy()); Ts.append(T.copy())
        if mm<.15:
            stop=n; break
    if stop is None: raise RuntimeError('threshold not reached')
    times=np.asarray(times); Cs=np.asarray(Cs); Ts=np.asarray(Ts); maxima=np.asarray(maxima)
    # Threshold crossing uses the last two full time levels.
    m0=maxima[stop-1]; m1=maxima[stop]
    cross=(stop-1)*dt+dt*(m0-.15)/(m0-m1)
    out=OUT; out.mkdir(parents=True,exist_ok=True); tag=f'N{N}_dt{dt:g}'
    np.savez_compressed(out/f'fields_{tag}.npz',times=times,C=Cs,T=Ts,radius=np.arange(N+1)*dr)
    summary={'schema_version':1,'question':'Q3','method':'Q2 inherited streaming solve',
      'N':N,'dt_s':dt,'role':'grid_diagnostic','drying_time_s':stop*dt,
      'drying_time_h':stop*dt/3600,'threshold_crossing_interpolated_s':cross,
      'previous_max_C':m0,'terminal_max_C':m1,'terminal_surface_C':float(C[-1]),
      'mean_iterations':mean_iter/stop,'runtime_s':time.perf_counter()-start,
      'frame_count':len(times),'max_radial_increase':float(np.max(np.diff(Cs,axis=1))),
      'finite_fields':bool(np.isfinite(Cs).all() and np.isfinite(Ts).all()),
      'ambient_endpoint':{'time_s':14400,'T_C':Tinf(14400),'C':Cinf(14400)},
      'inputs':{str(MODEL):hashlib.sha256(MODEL.read_bytes()).hexdigest(),str(ATT):hashlib.sha256(ATT.read_bytes()).hexdigest()},
      'environment':{'python':platform.python_version(),'numpy':np.__version__},
      'notes':['Only 60-second frames are stored; threshold is tested every time step.','Same three Picard/Newton iterations as Q2.']}
    (out/f'run_summary_{tag}.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int,required=True); ap.add_argument('--dt',type=float,default=1); x=ap.parse_args(); run(x.N,x.dt)
