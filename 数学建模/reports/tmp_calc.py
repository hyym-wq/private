import numpy as np
# parse grid values manually
seqs={
'uniform':([20,40,80,160,320,640],[202239.7266,204139.9219,205078.5938,205511.6016,205697.9297,205771.7578]),
'q1.5':([20,30,40,60,80],[205137.7734,205468.2422,205603.5938,205711.9922,205753.0078]),
'q2':([20,30,40,60,80,120,160],[205775.2734,205791.0938,205797.5391,205803.9844,205806.5625,205808.0859,205809.0234]),
'q2.5':([40,60,80],[205858.7109,205831.6406,205822.1484]),
'q3':([30,40,60],[205980.8203,205906.6406,None])}
from scipy.optimize import least_squares
for name,(N,y) in seqs.items():
 N=np.array(N,float); y=np.array(y,float); m=np.isfinite(y); N=N[m]; y=y[m]
 print('\n',name)
 if len(N)>=3:
  def fun(x):
   L,a,p=x; return L+a*N**(-p)-y
  # initial L near last, a diff, p1
  for x0 in ([y[-1]+(y[-1]-y[-2]),-1e5,1],[y[-1],1e5,1]):
   try:
    z=least_squares(fun,x0,max_nfev=100000).x; print('fit',z,'rmse',np.sqrt(np.mean(fun(z)**2)))
   except Exception as e: print(e)
 # pair subsequences ratio 2 exact fit 3 pts
 for i in range(len(N)-2):
  rat=N[i+1]/N[i]; rat2=N[i+2]/N[i+1]
  if abs(rat-rat2)<1e-8:
   d1=y[i+1]-y[i]; d2=y[i+2]-y[i+1]; p=np.log(abs(d1/d2))/np.log(rat); L=y[i+2]+d2/(rat**p-1); print('3pt',N[i:i+3],p,L)
# uniform all actual fit and equal ratio
