import numpy as np
D=np.load('results/p3_fields.npz'); t,C,Cd=D['times'],D['C'],D['C_out'];
for sec in [205740,205771.7578,205800,205860,205920]:
 i=np.where(t==sec)[0][0] if np.any(t==sec) else np.argmin(abs(t-sec)); print(sec,'idx',i,'t',t[i],'Ccenter',C[i,0],'Csurf',C[i,-1],'max',C[i].max(),'Cd',Cd[i,[0,5,10,15,20]])
print('C_end_ds vs t_min diff',D['C_end_ds']-Cd[np.where(t==205800)[0][0]]); print('vs t205860',D['C_end_ds']-Cd[np.where(t==205860)[0][0]])
