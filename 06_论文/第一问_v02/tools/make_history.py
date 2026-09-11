from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE=Path(__file__).resolve().parents[1]; ROOT=HERE.parents[1]
src=ROOT/'04_结果/第一问/q1_main_r16.npz'; out=HERE/'figures'; out.mkdir(exist_ok=True)
font=next((n for n in ('SimSun','Microsoft YaHei','SimHei') if n in {f.name for f in font_manager.fontManager.ttflist}), 'DejaVu Sans')
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':[font,'DejaVu Sans'],'font.size':9,'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42})
with np.load(src) as d:
    t=d['times']; T=d['temperature']; C=d['moisture']
fig,ax=plt.subplots(1,2,figsize=(16/2.54,6.0/2.54)); fig.subplots_adjust(left=.09,right=.98,bottom=.18,top=.88,wspace=.32)
for i,lab,col in [(0,'中心 r=0 cm','#0072B2'),(10,'r=1 cm','#E69F00'),(20,'表面 r=2 cm','#D55E00')]:
    ax[0].plot(t,T[:,i],lw=1.15,label=lab,color=col); ax[1].plot(t,C[:,i],lw=1.15,label=lab,color=col)
env=np.loadtxt(ROOT/'02_数据/清洗数据/0910_第一问_烘房边界_v01_Codex.csv',delimiter=',',skiprows=1)
ax[0].plot(env[:,0][env[:,0]<=1800],env[:,1][env[:,0]<=1800],'k--',lw=1,label='烘房温度')
for a,y,yl,title in [(ax[0],T,'温度 T / ℃','(a) 温度时程'),(ax[1],C,'含水率 C / (kg/kg)','(b) 含水率时程')]:
    a.set(xlabel='时间 t / s',ylabel=yl,title=title,xlim=(0,1800)); a.grid(axis='y',alpha=.25); a.spines[['top','right']].set_visible(False); a.tick_params(direction='in')
ax[0].set_ylim(28,43); ax[1].set_ylim(1.45,2.6); ax[0].legend(frameon=False,ncol=2,fontsize=8,loc='upper left'); ax[1].legend(frameon=False,fontsize=8,loc='lower left')
for ext,dpi in [('pdf',None),('png',400)]: fig.savefig(out/f'q1_history.{ext}',dpi=dpi,bbox_inches='tight' if ext=='png' else None)
plt.close(fig)
