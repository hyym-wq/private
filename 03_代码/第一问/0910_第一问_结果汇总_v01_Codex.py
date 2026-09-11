"""2026-09-10 A Q1: canonical result matrices and verification metrics."""
from pathlib import Path
import hashlib
import json
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO/'04_结果/第一问'
LABEL = 'q1_main_r16'


def compare(a, b):
    table_indices = np.array([100, 300, 600, 900, 1200, 1500, 1800])-1
    answer = {}
    for field in ('temperature', 'moisture'):
        difference = np.abs(a[field]-b[field])
        idx = np.unravel_index(difference.argmax(), difference.shape)
        answer[field] = {
            'all_outputs_max_abs': float(difference.max()),
            'table_max_abs': float(difference[table_indices][:, ::5].max()),
            'maximum_at_time_s_radius_cm': [int(idx[0]+1), float(idx[1]/10)],
        }
    return answer


def main():
    main_run = np.load(OUT/f'{LABEL}.npz')
    baseline = np.load(OUT/'q1_baseline_r16.npz')
    timecheck = np.load(OUT/'q1_reference_r16.npz')
    spacecheck = np.load(OUT/'q1_spacecheck_r32.npz')
    indices = np.array([100, 300, 600, 900, 1200, 1500, 1800])-1
    data = {
        'source_run': LABEL,
        'cross_section': 'z=0, finite-cylinder axial midplane',
        'times': main_run['times'].astype(int).tolist(),
        'radii_cm': (np.arange(21)/10).tolist(),
        'table_times': (indices+1).tolist(),
        'table_radii_cm': [0, .5, 1, 1.5, 2],
    }
    for field in ('temperature', 'moisture'):
        rounded = np.round(main_run[field], 4)
        data[field] = rounded.tolist()
        data[f'table_{field}'] = rounded[indices][:, ::5].tolist()
        np.savetxt(OUT/f'0910_第一问_{field}_论文表_v01_Codex.csv',
                   np.c_[indices+1, rounded[indices][:, ::5]], delimiter=',',
                   header='time_s,r_0_cm,r_0.5_cm,r_1_cm,r_1.5_cm,r_2_cm',
                   comments='', fmt=['%d']+['%.4f']*5)
    (OUT/'q1_export.json').write_text(json.dumps(data, ensure_ascii=False), encoding='utf-8')
    metrics = {
        'comparison_scope': '1800 output times x 21 radii, axial midplane',
        'main_vs_1d_identical_radial_time': compare(main_run, baseline),
        'time_refinement_1d_dt_half': compare(baseline, timecheck),
        'space_refinement_1d_radial_double': compare(timecheck, spacecheck),
        'main_vs_finer_1d_combined': compare(main_run, spacecheck),
        'final_volume_means_2d': main_run['balances'][-1, 1:3].tolist(),
        'final_radial_means_1d': baseline['balances'][-1, 1:3].tolist(),
        'final_volume_mean_2d_minus_1d': (
            main_run['balances'][-1, 1:3]-baseline['balances'][-1, 1:3]).tolist(),
        'scope_limit': 'Output convergence is verified at z=0; end-corner field accuracy is not certified to four decimals.',
        'raw_result_sha256': hashlib.sha256((OUT/f'{LABEL}.npz').read_bytes()).hexdigest(),
    }
    evidence = OUT/'experiments/round1'
    (evidence/'metrics').mkdir(parents=True, exist_ok=True)
    (evidence/'metrics/q1_convergence.json').write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2), encoding='utf-8')
    main_summary = json.loads((OUT/f'{LABEL}_summary.json').read_text(encoding='utf-8'))
    baseline_summary = json.loads((OUT/'q1_baseline_r16_summary.json').read_text(encoding='utf-8'))
    run_summary = {
        'schema_version': 1, 'question': 'Q1', 'round': 'round1',
        'implementation_target': 'python', 'random_seed': None,
        'determinism': 'No random algorithm used.',
        'approved_decision_id': 'q1_method_choice',
        'output_scope_confirmation': 'q1_midplane_output_confirmation',
        'methods': [], 'comparison': metrics,
        'fallback_trigger': {'fallback_id': None, 'condition': None, 'observed': False,
                             'evidence': 'Graded startup and linear preconditioning refine the approved solver.'},
        'environment': main_summary['environment'],
        'accuracy_statement': 'Four decimal places are an output format. Convergence differences are numerical estimates, not observational validation.',
    }
    for method_id, role, summary in [
        ('Q1_MAIN_2D_FVM_CN', 'main', main_summary),
        ('Q1_BASELINE_1D_FVM_CN', 'usable_baseline', baseline_summary)]:
        run_summary['methods'].append({
            'method_id': method_id, 'role': role,
            'script': '03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py',
            'status': 'success', 'execution_time_seconds': summary['elapsed_seconds'],
            'input_files': ['02_数据/清洗数据/0910_第一问_烘房边界_v01_Codex.csv'],
            'output_files': [f"04_结果/第一问/{summary['label']}.npz"],
            'figure_files': [], 'metrics_summary': summary,
            'warnings': ['Effective environmental moisture boundary; no latent heat or sorption closure supplied.'],
            'errors': [],
        })
    main_record = run_summary['methods'][0]
    main_record['output_files'] += [
        '04_结果/第一问/result1.xlsx',
        '04_结果/第一问/q1_export.json',
        '04_结果/第一问/0910_第一问_结果与模型说明_v01_Codex.md',
        '04_结果/第一问/0910_第一问_temperature_论文表_v01_Codex.csv',
        '04_结果/第一问/0910_第一问_moisture_论文表_v01_Codex.csv',
    ]
    main_record['figure_files'] = [
        '05_图片/第一问/0910_第一问_中截面径向分布_v01_Codex.png',
        '05_图片/第一问/0910_第一问_1800s二维场诊断_v01_Codex.png',
        '05_图片/第一问/0910_第一问_二维与一维平均量对照_v01_Codex.png',
    ]
    main_record['output_degeneracy_evidence'] = {
        'max_midplane_temperature_change_from_initial': float(np.max(main_run['temperature']-28)),
        'max_midplane_moisture_change_from_initial': float(np.max(2.55-main_run['moisture'])),
        'interpretation': 'Near-constant center moisture is physical early-time diffusion, not solver stagnation; surface values change substantially.',
    }
    (evidence/'run_summary.json').write_text(
        json.dumps(run_summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
