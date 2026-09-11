"""Targeted Q1 ablations: same physics, comparable grids, fixed output scope.

Original solver is imported without changing its source. Uniform-grid override
is limited to this process. Results use separate labels, preserving main output.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/'03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py'
spec = importlib.util.spec_from_file_location('q1_fvm_original', SOURCE)
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
OUT = ROOT/'04_结果/第一问'
OriginalGrid = M.Grid


class UniformGrid(OriginalGrid):
    def __init__(self, refinement=16, dimension=1, axial_refinement=None):
        if dimension != 1:
            raise ValueError('Uniform ablation is radial only.')
        self.r = np.linspace(0, M.R, 881)
        self.z = np.array([0.])
        rf = np.r_[0, (self.r[:-1]+self.r[1:])/2, M.R]
        self.volume = np.diff(rf**2)/2
        self.sqrtv = np.sqrt(self.volume)
        self.shape = (1, len(self.r))
        self.size = len(self.r)
        self.p, self.q = np.arange(self.size-1), np.arange(1,self.size)
        self.geom = rf[1:-1]/np.diff(self.r)
        self.side = np.zeros(self.size)
        self.side[-1] = M.R
        self.end = np.zeros(self.size)
        self.boundary = self.side.copy()
        self.rows = np.r_[self.p,self.q,np.arange(self.size)]
        self.cols = np.r_[self.q,self.p,np.arange(self.size)]


def run():
    plan = {'model_change': 'none', 'canonical_main': 'q1_main_r16',
        'predeclared_comparison_scope': 't=1..1800 s, 21 radii; additionally 1..60 s and 7x5 paper table',
        'reference': 'q1_spacecheck_r32: 1729 radial nodes, dt=.25 s, graded startup',
        'numerical_resolution_target': {'temperature_C':5e-5,'moisture_kg_kg':5e-5,
            'meaning':'Half of one 4-decimal display unit as a numerical difference target, not physical uncertainty or a guarantee of identical rounding.'},
        'cases': ['865 graded nodes + graded startup (saved baseline)',
                  '881 uniform nodes + same graded startup',
                  '865 graded nodes + two BE half steps totaling .5 s then fixed .5 s CN'],
        'axial_check': '865 radial nodes fixed, 25 to 49 axial nodes, same .5 s and startup; separate output label',
        'solver_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
    (OUT/'q1_improvement_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
    label = 'q1_uniform881_graded_start'
    if not (OUT/(label+'_summary.json')).exists():
        M.Grid = UniformGrid
        try:
            M.integrate(16,1,.5,label,1,True)
        finally:
            M.Grid = OriginalGrid
    label = 'q1_graded865_standard_start'
    if not (OUT/(label+'_summary.json')).exists():
        M.integrate(16,1,.5,label,1,False)
    summarize()


def diff(a,b):
    result = {}
    table_t = np.array([100,300,600,900,1200,1500,1800])-1
    for k in ['temperature','moisture']:
        d = np.abs(a[k]-b[k])
        idx = np.unravel_index(np.argmax(d),d.shape)
        result[k] = {'all_max_abs':float(d.max()),'early_1_60s_max_abs':float(d[:60].max()),
            'paper_7x5_max_abs':float(d[table_t][:,::5].max()),
            'maximum_time_s_radius_cm':[int(idx[0]+1),float(idx[1]/10)]}
    return result


def summarize():
    ref = np.load(OUT/'q1_spacecheck_r32.npz')
    cases = {}
    for label in ['q1_baseline_r16','q1_uniform881_graded_start','q1_graded865_standard_start']:
        data = np.load(OUT/(label+'.npz'))
        s = json.loads((OUT/(label+'_summary.json')).read_text(encoding='utf-8'))
        cases[label] = {'comparison':diff(data,ref), 'nodes':len(data['r']),
            'steps':s['steps'],'max_nonlinear_residual':s['max_nonlinear_residual_C'],
            'min_moisture':s['moisture_bounds'][0],
            'source_sha256':hashlib.sha256((OUT/(label+'.npz')).read_bytes()).hexdigest()}
    main = np.load(OUT/'q1_main_r16.npz')
    report = {'status':'complete_radial_ablations', 'cases':cases,
        'predeclared_plan':'q1_improvement_plan.json',
        'comparison_is_not_analytic_error_bound':True,
        'axial_status':'pending'}
    axialpath = OUT/'q1_axialcheck_r16_z2.npz'
    if axialpath.exists() and (OUT/'q1_axialcheck_r16_z2_summary.json').exists():
        finer = np.load(axialpath)
        assert np.allclose(main['r'],finer['r']) and np.allclose(main['z'],finer['z'][::2])
        report['axial_status'] = 'complete'
        report['axial'] = {'midplane':diff(main,finer),
            'whole_snapshot_max_abs_T':float(np.abs(main['snapshots_t']-finer['snapshots_t'][:,::2,:]).max()),
            'whole_snapshot_max_abs_C':float(np.abs(main['snapshots_c']-finer['snapshots_c'][:,::2,:]).max()),
            'final_mean_finer':finer['balances'][-1,1:3].tolist(),
            'final_mean_finer_minus_main':(finer['balances'][-1,1:3]-main['balances'][-1,1:3]).tolist(),
            'source_sha256':hashlib.sha256(axialpath.read_bytes()).hexdigest(),
            'limitation':'One axial refinement; supports difference assessment, not full-domain extrapolated accuracy.'}
    target=OUT/'q1_improvement_evidence.json'
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=True,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--summarize',action='store_true')
    args=parser.parse_args()
    summarize() if args.summarize else run()
