"""Independent Q3 numerical experiment: unchanged equations, adaptive CN."""
from pathlib import Path
import argparse, importlib.util, json, time, sys
import numpy as np
from scipy.linalg import solve_banded
sys.path.append('C:/Users/29945/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/Lib/site-packages')

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('inherited', ROOT / '数学建模/数学建模/code/p2_model.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class Solver:
    def __init__(self, n):
        self.n=n; self.dr=m.R_CYL/n; self.r=np.arange(n+1)*self.dr
        self.rf=(self.r[1:]+self.r[:-1])/2
        self.w=m.conservative_weights(self.dr,n)
        self.left=np.zeros(n+1); self.right=np.zeros(n+1)
        self.left[1:]=self.rf/(self.w[1:]*self.dr)
        self.right[:-1]=self.rf/(self.w[:-1]*self.dr)
        self.gT=m.R_CYL*m.H_HEAT/self.w[-1]
        self.gC=m.R_CYL*m.H_MASS/self.w[-1]
        att=ROOT.parent/'CUMCM2026Problems/A题/附件/附件1.xlsx'
        self.ta,self.Ta,self.Ca=m.load_ambient(att)
        self.linear_solves=0; self.max_res_T=0.; self.max_res_C=0.; self.max_iters=0
    def env(self,t):
        return np.interp(t,self.ta,self.Ta),np.interp(t,self.ta,self.Ca)
    def coefficients(self,C):
        inv=1/(m.rho(C)*m.cp(C)); k=m.kcond(C); kh=(k[1:]+k[:-1])/2
        a=self.left*np.r_[0.,kh]*inv; c=self.right*np.r_[kh,0.]*inv
        b=-a-c; b[-1]-=self.gT*inv[-1]
        return a,b,c,self.gT*inv[-1]
    @staticmethod
    def action(a,b,c,u):
        v=b*u; v[1:]+=a[1:]*u[:-1]; v[:-1]+=c[:-1]*u[1:]; return v
    def ft(self,T,C,t):
        a,b,c,g=self.coefficients(C); out=self.action(a,b,c,T); out[-1]+=g*self.env(t)[0]; return out
    def fc(self,T,C,t):
        D=m.diff_coef(C,T); Dh=(D[1:]+D[:-1])/2
        flux=Dh*(C[1:]-C[:-1]); out=np.zeros_like(C)
        out[:-1]+=self.right[:-1]*flux; out[1:]-=self.left[1:]*flux
        out[-1]-=self.gC*(C[-1]-self.env(t)[1]); return out
    def jac(self,T,C):
        D=m.diff_coef(C,T); Dp=m.dD_dC(C,T); Dh=(D[1:]+D[:-1])/2
        diff=C[1:]-C[:-1]
        # flux(C_l,C_r)=D_half*(C_r-C_l)
        ql=-Dh+.5*Dp[:-1]*diff; qr=Dh+.5*Dp[1:]*diff
        a=np.zeros_like(C); b=np.zeros_like(C); c=np.zeros_like(C)
        a[1:]=-self.left[1:]*ql; c[:-1]=self.right[:-1]*qr
        b[:-1]+=self.right[:-1]*ql; b[1:]-=self.left[1:]*qr
        b[-1]-=self.gC; return a,b,c
    def tridiag(self,a,b,c,rhs):
        ab=np.zeros((3,len(b))); ab[0,1:]=c[:-1]; ab[1]=b; ab[2,:-1]=a[1:]
        self.linear_solves+=1
        return solve_banded((1,1),ab,rhs,check_finite=False,overwrite_ab=True,overwrite_b=False)
    def step(self,T0,C0,t,dt):
        ft0=self.ft(T0,C0,t); fc0=self.fc(T0,C0,t)
        T=T0.copy(); C=C0.copy(); hn=dt/2; tn=t+dt
        for it in range(20):
            a,b,c,g=self.coefficients(C)
            rhs=T0+hn*ft0; rhs[-1]+=hn*g*self.env(tn)[0]
            Tn=self.tridiag(-hn*a,1-hn*b,-hn*c,rhs)
            fc=self.fc(Tn,C,tn); a,b,c=self.jac(Tn,C)
            residual=C-C0-hn*(fc0+fc)
            delta=self.tridiag(-hn*a,1-hn*b,-hn*c,-residual)
            lam=1.
            neg=delta < -1e-20
            if np.any(neg): lam=min(1.,.99*float(np.min(-C[neg]/delta[neg])))
            Cn=C+lam*delta
            if not np.isfinite(Tn).all() or not np.isfinite(Cn).all(): raise RuntimeError('nonfinite')
            rt=np.max(np.abs(Tn-T0-hn*(ft0+self.ft(Tn,Cn,tn))))
            rc=np.max(np.abs(Cn-C0-hn*(fc0+self.fc(Tn,Cn,tn))))
            T,C=Tn,Cn
            if rt<2e-10 and rc<2e-12:
                self.max_res_T=max(self.max_res_T,float(rt)); self.max_res_C=max(self.max_res_C,float(rc))
                self.max_iters=max(self.max_iters,it+1)
                return T,C
        raise RuntimeError('nonlinear convergence failure')
    def double(self,T,C,t,h):
        T1,C1=self.step(T,C,t,h)
        Th,Ch=self.step(T,C,t,h/2)
        T2,C2=self.step(Th,Ch,t+h/2,h/2)
        return T1,C1,T2,C2,Th,Ch

def run(n=160,tol=1e-7,maxdt=60.,tend=216000.,output=None):
    s=Solver(n); T=np.full(n+1,m.T0); C=np.full(n+1,m.C0)
    t=0.; h=.25; accepted=0; rejected=0; nextout=60.; records=[]; table=[]
    total_loss=0.; mass0=float(np.dot(s.w,C)); minC=float(C.min()); maxinc=0.
    start=time.perf_counter(); event=None; maxerr=0.
    while t<tend:
        future=s.ta[s.ta>t+1e-7]
        knot=float(future[0]) if len(future) else tend
        h=min(h,maxdt,tend-t,nextout-t,knot-t)
        if h<1e-8:
            if nextout-t<1e-7: nextout+=60.
            continue
        try:
            T1,C1,T2,C2,Th,Ch=s.double(T,C,t,h)
            # Local error estimate of two half steps for order-2 CN.
            et=np.max(np.abs(T2-T1)/(3*(tol*100+tol*np.maximum(np.abs(T),np.abs(T2)))))
            ec=np.max(np.abs(C2-C1)/(3*(tol+tol*np.maximum(np.abs(C),np.abs(C2)))))
            err=max(float(et),float(ec))
            if C2.min()<-1e-12 or T2.min()<min(m.T0,s.env(t+h)[0])-1e-7: err=max(err,2.)
        except RuntimeError:
            err=np.inf
        if err>1:
            rejected+=1; h*=max(.1,.8*err**(-1/3)) if np.isfinite(err) else .5
            if h<1e-5: raise RuntimeError('step size underflow')
            continue
        accepted+=1; maxerr=max(maxerr,err)
        # Both accepted half steps use CN, so their boundary quadrature telescopes.
        flux0=m.R_CYL*m.H_MASS*(C[-1]-s.env(t)[1])
        fluxh=m.R_CYL*m.H_MASS*(Ch[-1]-s.env(t+h/2)[1])
        flux1=m.R_CYL*m.H_MASS*(C2[-1]-s.env(t+h)[1])
        step_loss=h/4*(flux0+2*fluxh+flux1)
        total_loss+=step_loss
        if C.max()>=.15 and C2.max()<.15:
            lo=0.; hi=h; C_hi=C2; T_hi=T2
            for _ in range(30):
                mid=(lo+hi)/2
                _,_,Tm,Cm,_,_=s.double(T,C,t,mid)
                if Cm.max()<.15: hi=mid; C_hi=Cm; T_hi=Tm
                else: lo=mid
                if hi-lo<1e-3: break
            _,_,T_hi,C_hi,T_mid,C_mid=s.double(T,C,t,hi)
            fluxmid=m.R_CYL*m.H_MASS*(C_mid[-1]-s.env(t+hi/2)[1])
            fluxend=m.R_CYL*m.H_MASS*(C_hi[-1]-s.env(t+hi)[1])
            total_loss-=step_loss
            total_loss+=hi/4*(flux0+2*fluxmid+fluxend)
            event={'time_s':t+hi,'time_h':(t+hi)/3600,'bracket_s':[t+lo,t+hi],
                   'C_max':float(C_hi.max()),'C_surface':float(C_hi[-1])}
            records.append((t+hi,T_hi.copy(),C_hi.copy())); table.append((t+hi,C_hi[::n//4].copy()))
            t+=hi; T,C=T_hi,C_hi
            break
        t+=h; T,C=T2,C2; minC=min(minC,float(C.min())); maxinc=max(maxinc,float(np.diff(C).max()))
        if abs(t-nextout)<1e-6:
            records.append((t,T.copy(),C.copy())); nextout+=60.
            if abs(t/21600-round(t/21600))<1e-9: table.append((t,C[::n//4].copy()))
        factor=2. if err<1e-12 else np.clip(.9*err**(-1/3),.5,2.)
        h*=factor
    conservation=float((np.dot(s.w,C)-mass0+total_loss)/mass0)
    result={'N':n,'tol':tol,'maxdt_s':maxdt,'elapsed_s':time.perf_counter()-start,
            'accepted_steps':accepted,'rejected_steps':rejected,'linear_solves':s.linear_solves,
            'max_nonlinear_iterations':s.max_iters,'max_equation_residual_T':s.max_res_T,
            'max_equation_residual_C':s.max_res_C,'max_local_error_ratio':maxerr,
            'event':event,'end_time_s':t,'C_center':float(C[0]),'C_surface':float(C[-1]),
            'min_C':minC,'max_radial_increase':maxinc,'relative_mass_balance_residual':conservation,
            'table_times_h':[x[0]/3600 for x in table],'table_C':[x[1].tolist() for x in table]}
    if output:
        out=Path(output); out.parent.mkdir(parents=True,exist_ok=True)
        out.with_suffix('.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        np.savez_compressed(out.with_suffix('.npz'),times=[x[0] for x in records],r=s.r,
                            T=[x[1] for x in records],C=[x[2] for x in records])
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=160)
    p.add_argument('--tol',type=float,default=1e-7); p.add_argument('--maxdt',type=float,default=60.)
    p.add_argument('--tend',type=float,default=216000.); p.add_argument('--output')
    a=p.parse_args(); run(a.N,a.tol,a.maxdt,a.tend,a.output)
