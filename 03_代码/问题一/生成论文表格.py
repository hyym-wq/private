from pathlib import Path
import json,csv
root=Path(__file__).resolve().parents[2];out=root/'06_论文/问题一'
r=json.loads((root/'04_结果/问题一/统计与验证.json').read_text(encoding='utf-8'))
def table(caption,label,widths,header,rows):
    col=''.join('P{'+str(w)+'}' for w in widths)
    title=' & '.join(header)+r'\\'
    return '\\begin{longtable}{'+col+'}\n\\caption{'+caption+'}\\label{'+label+'}\\\\\n\\toprule '+title+'\n\\midrule\\endfirsthead\n\\toprule '+title+'\n\\midrule\\endhead\n\\bottomrule\\endfoot\n'+'\n'.join(' & '.join(str(v) for v in row)+r'\\' for row in rows)+'\n\\end{longtable}\n'
rows=[]
for key,name in [('Euler_40_0.25','Euler对照'),('CN_40_0.25','CN同网格对照'),('CN_80_0.25','CN主结果')]:
 c=r['runs'][key];rows.append([name,c['nodes'],c['dt'],f"{c['seconds']:.3f}",c['max_iterations']])
s=table('算法配置与本次运行时间','tab:method',[.28,.14,.14,.14,.14],['方法','节点数','步长/s','耗时/s','最多迭代'],rows)
rows=[]
for key,label in [('time_1_05','时间：1与0.5秒'),('time_05_025','时间：0.5与0.25秒'),('space_20_40','空间：21与41节点'),('space_40_80','空间：41与81节点')]:
 rows.append([label,f"{r[key]['T']:.7f}",f"{r[key]['C']:.7f}"])
s+=table('CN相邻分辨率全时空最大差异','tab:convergence',[.42,.25,.25],['比较','温度差/摄氏度','含水率差/(kg/kg)'],rows)
(out/'数值比较表.tex').write_text(s,encoding='utf-8')
s=''
for name,label,caption in [('温度','temperature','指定时刻药材温度（摄氏度）'),('水分浓度','moisture','指定时刻药材干基含水率（kg/kg）')]:
 with (root/f'04_结果/问题一/指定时刻{name}.csv').open(encoding='utf-8-sig') as f:
  rows=list(csv.reader(f))[1:]
 s+=table(caption,'tab:'+label,[.10,.14,.14,.14,.14,.14],['时间/s','0 cm','0.5 cm','1 cm','1.5 cm','2 cm'],rows)
(out/'指定点结果表.tex').write_text(s,encoding='utf-8')
