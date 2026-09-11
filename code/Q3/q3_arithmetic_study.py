"""Independently refine original arithmetic-flux discretisation via BDF."""
from pathlib import Path
import json
import q3_kirchhoff as k

if __name__ == '__main__':
    outputs=[]
    for n in [20,40,80,160,320,640,1280]:
        outputs.append(k.run(n,1.,'arithmetic',rtol=1e-9,max_step=180.))
    (k.OUT/'arithmetic_study.json').write_text(json.dumps(outputs,indent=2),encoding='utf8')
