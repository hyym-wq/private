import os,sys,time, json
import numpy as np
sys.path.insert(0,'code'); import p2_model as m
path=os.path.abspath(os.path.join('A题','附件','附件1.xlsx'))
Ti,Ci=m.make_ambient(path)
N=20; dt=1.; R=m.R_CYL; dr=R/N; theta=.5; pmx=3; tol=1e-9; nsteps=10800
T=np.full(N+1,m.T0); C=np.full(N+1,m.C0); zero=np.zeros(N+1)
iter_counts=[]; dT_hist=[]; dC_hist=[]; converged9=[]; converged8=[]
t0=time.time()
for n in range(1,nsteps+1):
 tn,to=n*dt,(n-1)*dt; Tin,Tio=Ti(tn),Ti(to); Cin,Cio=Ci(tn),Ci(to)
 T_old,C_old=T.copy(),C.copy(); C_star,T_star=C.copy(),T.copy(); it=pmx; ds=[]
 ah_o,bh_o,ch_o,gT_o=m.build_op_heat(dr,N,C_old,R); sT_o=zero.copy(); sT_o[N]=gT_o*Tio
 for k in range(pmx):
  ah_n,bh_n,ch_n,gT_n=m.build_op_heat(dr,N,C_star,R); sT_n=zero.copy(); sT_n[N]=gT_n*Tin
  T_new=m._step(T_old,(ah_n,bh_n,ch_n),(ah_o,bh_o,ch_o),sT_o,sT_n,dt,theta)
  C_new=m._mass_step(C_old,C_star,T_old,T_new,Cio,Cin,dt,theta,dr,N,R,True)
  dT=float(np.max(np.abs(T_new-T_star))); dC=float(np.max(np.abs(C_new-C_star))); ds.append((dT,dC)); T_star,C_star=T_new,C_new
  if dT<tol and dC<tol: it=k+1; break
 iter_counts.append(it); dT_hist.append(ds[-1][0]); dC_hist.append(ds[-1][1]); converged9.append(it<pmx or (ds[-1][0]<tol and ds[-1][1]<tol)); converged8.append(ds[-1][0]<1e-8 and ds[-1][1]<1e-8)
 T,C=T_new,C_new
# summarize
arrT=np.array(dT_hist); arrC=np.array(dC_hist); ic=np.array(iter_counts)
out={'path':path,'t_end':nsteps*dt,'N':N,'dt':dt,'picard_max':pmx,'tol_default':tol,'mean_iters':float(ic.mean()),'counts':{str(i):int(np.sum(ic==i)) for i in range(1,pmx+1)},'not_converged_default':int(np.sum(ic==pmx)), 'last_dT_quantiles':np.quantile(arrT,[0,.5,.9,.99,1]).tolist(),'last_dC_quantiles':np.quantile(arrC,[0,.5,.9,.99,1]).tolist(),'steps_last_change_gt_1e-8':int(np.sum((arrT>=1e-8)|(arrC>=1e-8))), 'steps_last_change_gt_1e-9':int(np.sum((arrT>=1e-9)|(arrC>=1e-9))), 'first10':[(i+1,int(ic[i]),float(arrT[i]),float(arrC[i])) for i in range(min(10,nsteps))], 'elapsed_s':time.time()-t0,'final_center_C':float(C[0]),'final_surface_C':float(C[-1])}
os.makedirs('reports',exist_ok=True)
json.dump(out,open('reports/p2_picard_short.json','w',encoding='utf8'),ensure_ascii=False,indent=2)
print(json.dumps(out,ensure_ascii=False,indent=2))

