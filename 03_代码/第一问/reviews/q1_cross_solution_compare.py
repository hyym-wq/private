"""Read-only comparison of saved Q1 solutions; does not rerun or overwrite solvers."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / '04_结果/第一问'
paths = {
    'ours': OUT / 'q1_main_r16.npz',
    'euler_cn_81_nodes': ROOT / '04_结果/问题一/数值场.npz',
    'dimensionless_21_nodes': ROOT / '数学建模/results/p1_fields.npz',
    'dimensionless_nested_copy': ROOT / '数学建模/数学建模/results/p1_fields.npz',
}
main = np.load(paths['ours'])
times = main['times']
r_cm = np.arange(21) * .1
result = {'scope': 'Saved results compared at t=1..1800 s and r=0..2 cm in steps of 0.1 cm; ours at z=0. Differences are not exact errors.',
          'sources': {k: {'path': str(p.relative_to(ROOT)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},
          'comparisons': {}}
for name, stride in [('euler_cn_81_nodes', 4), ('dimensionless_21_nodes', 1)]:
    other = np.load(paths[name])
    assert np.allclose(other['r'][::stride]*100, r_cm)
    assert np.array_equal(other['time' if name.startswith('euler') else 'times'][1:], times)
    metrics = {}
    for source_key, main_key in [('T', 'temperature'), ('C', 'moisture')]:
        values = other[source_key][1:, ::stride]
        error = values - main[main_key]
        idx = np.unravel_index(np.argmax(np.abs(error)), error.shape)
        indices = np.array([100, 300, 600, 900, 1200, 1500, 1800])-1
        metrics[source_key] = {
            'all_outputs_max_abs_difference': float(np.max(np.abs(error))),
            'max_at_time_s_radius_cm': [int(times[idx[0]]), float(r_cm[idx[1]])],
            'table_7x5_max_abs_difference': float(np.abs(error[indices][:, ::5]).max()),
            'final_surface': float(values[-1,-1]),
            'our_final_surface': float(main[main_key][-1,-1]),
            'surface_t100': float(values[99,-1]),
            'our_surface_t100': float(main[main_key][99,-1]),
        }
    result['comparisons'][name] = metrics
outer = np.load(paths['dimensionless_21_nodes'])
nested = np.load(paths['dimensionless_nested_copy'])
result['nested_solution_equals_outer'] = all(np.array_equal(outer[k],nested[k]) for k in outer.files)
R,L,t,C0,k,rho,cp,h,hm = .02,.25,1800,2.55,.36,820,2600,25,8e-7
alpha = k/(rho*cp)
D0 = 7e-9*np.exp(-.89/C0)
result['dimensionless'] = dict(alpha_m2_s=alpha, initial_D_m2_s=D0,
    diffusivity_ratio=alpha/D0, Bi_heat_radius=h*R/k, Bi_mass_radius=hm*R/D0,
    Fo_heat_radius=alpha*t/R**2, Fo_mass_radius=D0*t/R**2,
    heat_diffusion_length_cm=float(np.sqrt(alpha*t)*100),
    moisture_diffusion_length_cm=float(np.sqrt(D0*t)*100),
    half_length_over_heat_diffusion_length=(L/2)/np.sqrt(alpha*t),
    end_to_side_area_ratio=R/L)
r = main['r']; C = main['snapshots_c'][-1,0,:]
front = {}
for epsilon in [.001,.01,.05]:
    threshold = (1-epsilon)*C0
    radius = float(np.interp(threshold, C[::-1], r[::-1]))
    front[str(epsilon)] = {'radius_cm': radius*100, 'shell_thickness_cm': (R-radius)*100}
result['diagnostic_midplane_threshold_fronts_at_1800s'] = front
result['front_note'] = 'Threshold-defined diagnostic, not a sharp physical interface; piecewise-linear interpolation on saved 865-node radial snapshot.'
result['our_mean_moisture_reduction_percent'] = (1-main['balances'][-1,2]/C0)*100
result['our_center_ambient_gap_C'] = 41.513-main['temperature'][-1,0]
result['our_surface_ambient_gap_C'] = 41.513-main['temperature'][-1,-1]
result['our_radial_temperature_difference_C'] = main['temperature'][-1,-1]-main['temperature'][-1,0]
target = OUT/'0910_第一问_三方案对照指标_v01_Codex.json'
target.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='sources'}, ensure_ascii=True, indent=2))
