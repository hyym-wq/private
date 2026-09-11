from q3_kirchhoff import run, OUT
import json

results=[]
for p in [1.,1.5,2.]:
    for n in [20,40,80,160,320,640]:
        results.append(run(n,p,'kirchhoff',1e-8,300))
for n in [40,160,640]:
    results.append(run(n,1.,'arithmetic',1e-8,300))
for tol in [1e-9,1e-10]:
    results.append(run(640,1.5,'kirchhoff',tol,120))
(OUT/'study.json').write_text(json.dumps(results,indent=2),encoding='utf8')
