"""Build Q1 LaTeX tables and a traceable writing snapshot, without rerunning the model."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

PAPER = Path(__file__).resolve().parents[1]
REPO = PAPER.parents[1]
RESULT = REPO / '04_结果/第一问'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sci(value, digits=3):
    mantissa, exponent = f'{value:.{digits-1}e}'.split('e')
    return rf'{mantissa}\times10^{{{int(exponent)}}}'


def write_table(name, caption, label, header, rows, widths=None):
    columns = widths or ('l' + 'r' * (len(rows[0]) - 1))
    text = [r'\begin{table}[H]', r'\centering\small',
            rf'\caption{{{caption}}}\label{{{label}}}',
            rf'\begin{{tabular}}{{{columns}}}', r'\toprule', header,
            r'\midrule']
    text.extend(' & '.join(row) + r' \\' for row in rows)
    text.extend([r'\bottomrule', r'\end{tabular}', r'\end{table}', ''])
    (PAPER / 'tables' / name).write_text('\n'.join(text), encoding='utf-8')


def main():
    (PAPER / 'tables').mkdir(exist_ok=True)
    (PAPER / 'evidence').mkdir(exist_ok=True)
    paths = {
        'export': RESULT / 'q1_export.json',
        'summary': RESULT / 'q1_main_r16_summary.json',
        'convergence': RESULT / 'experiments/round1/metrics/q1_convergence.json',
        'independent': REPO / '03_代码/第一问/reviews/q1_independent_checks.json',
    }
    ex, summary, convergence, independent = (read(paths[k]) for k in paths)
    assert ex['source_run'] == summary['label'] == 'q1_main_r16'
    assert len(ex['times']) == 1800 and len(ex['radii_cm']) == 21
    claims = []
    macros = []

    def claim(name, value, unit, source, locator, display):
        claims.append(dict(claim_id=name, value=value, unit=unit,
                           source_file=str(paths[source].relative_to(REPO)).replace('\\', '/'),
                           source_locator=locator, latex_display=display))
        macros.append(rf'\newcommand{{\{name}}}{{{display}}}')

    for key, caption, label in (
        ('temperature', r'30分钟内中截面温度（单位：$\celsius$）', 'tab:temperature'),
        ('moisture', r'30分钟内中截面干基含水率（单位：$\kgkg$）', 'tab:moisture'),
    ):
        rows = []
        for i, (time, values) in enumerate(zip(ex['table_times'], ex['table_' + key])):
            assert len(values) == 5
            # Ensure official table cells also agree with the full exported matrix.
            full = ex[key][time-1]
            assert values == [full[j] for j in [0, 5, 10, 15, 20]]
            rows.append([str(time)] + [f'{v:.4f}' for v in values])
        header = (r'\multirow{2}{*}{时间 / s} & \multicolumn{5}{c}{到药材中心的距离 / cm}\\'
                  '\n' + r'\cmidrule(lr){2-6} & 0 & 0.5 & 1 & 1.5 & 2\\')
        write_table(key + '.tex', caption, label, header, rows, 'c' + 'r'*5)
        claims.append(dict(claim_id='table_' + key, value=ex['table_' + key],
                           unit='Celsius' if key == 'temperature' else 'kg/kg dry basis',
                           source_file=str(paths['export'].relative_to(REPO)).replace('\\', '/'),
                           source_locator='$.table_' + key,
                           time_s=ex['table_times'], radii_cm=ex['table_radii_cm']))

    for name, field, pos, unit in [
        ('FinalCenterT', 'table_temperature', 0, 'Celsius'),
        ('FinalSurfaceT', 'table_temperature', -1, 'Celsius'),
        ('FinalCenterC', 'table_moisture', 0, 'kg/kg'),
        ('FinalSurfaceC', 'table_moisture', -1, 'kg/kg'),
    ]:
        value = ex[field][-1][pos]
        claim(name, value, unit, 'export', f'$.{field}[6][{pos % 5}]', f'{value:.4f}')

    specifications = [
        ('NonlinearResidual', summary['max_nonlinear_residual_C'], 'kg/kg', 'summary', '$.max_nonlinear_residual_C'),
        ('HeatResidual', summary['max_heat_balance_residual_mean_Celsius'], 'Celsius', 'summary', '$.max_heat_balance_residual_mean_Celsius'),
        ('WaterResidual', summary['max_water_balance_residual_mean_kg_kg'], 'kg/kg', 'summary', '$.max_water_balance_residual_mean_kg_kg'),
        ('SeriesDifference', independent['analytic_reference']['comparison_256_512_modes_max_Celsius'], 'Celsius', 'independent', '$.analytic_reference.comparison_256_512_modes_max_Celsius'),
        ('AnalyticError', independent['saved_runs']['q1_main_r16']['heat_error_analytic_max_Celsius'], 'Celsius', 'independent', '$.saved_runs.q1_main_r16.heat_error_analytic_max_Celsius'),
    ]
    for name, value, unit, source, locator in specifications:
        claim(name, value, unit, source, locator, sci(value))
    for index, (name, unit) in enumerate([('MeanDeltaT', 'Celsius'), ('MeanDeltaC', 'kg/kg')]):
        value = convergence['final_volume_mean_2d_minus_1d'][index]
        claim(name, value, unit, 'convergence', f'$.final_volume_mean_2d_minus_1d[{index}]', f'{value:+.4f}')

    checks = [
        ('二维中截面与同网格一维对照', 'main_vs_1d_identical_radial_time'),
        ('一维时间步长减半', 'time_refinement_1d_dt_half'),
        ('一维径向网格加密一倍', 'space_refinement_1d_radial_double'),
        ('二维中截面与更细一维参考', 'main_vs_finer_1d_combined'),
    ]
    rows = []
    for label, key in checks:
        vals = [convergence[key][u]['all_outputs_max_abs'] for u in ['temperature', 'moisture']]
        rows.append([label] + ['$' + sci(v) + '$' for v in vals])
        claims.append(dict(claim_id=key, value=vals, unit=['Celsius', 'kg/kg'],
                           source_file=str(paths['convergence'].relative_to(REPO)).replace('\\', '/'),
                           source_locator=f'$.{key}.{{temperature,moisture}}.all_outputs_max_abs',
                           limitation='z=0 at 1800 x 21 outputs; not an analytic error bound'))
    write_table('verification.tex', '全部中截面输出点的最大绝对差', 'tab:verification',
                r'比较内容 & 温度 / $\celsius$ & 含水率 / $(\kgkg)$\\', rows)
    rows = []
    for label, key in [('二维有限圆柱', 'final_volume_means_2d'),
                       ('一维径向对照', 'final_radial_means_1d'),
                       ('二维减一维', 'final_volume_mean_2d_minus_1d')]:
        rows.append([label] + [f'{v:.4f}' for v in convergence[key]])
        claims.append(dict(claim_id=key, value=convergence[key], unit=['Celsius', 'kg/kg'],
                           source_file=str(paths['convergence'].relative_to(REPO)).replace('\\', '/'),
                           source_locator=f'$.{key}', limitation='Full-volume means used as structural comparison, not four-decimal certification'))
    write_table('means.tex', r'$1800\,\mathrm{s}$ 时的体积平均状态对照', 'tab:means',
                r'模型 & 平均温度 / $\celsius$ & 平均含水率 / $(\kgkg)$\\', rows)
    (PAPER / 'tables/numbers.tex').write_text('\n'.join(macros) + '\n', encoding='utf-8')
    snapshot = {
        'status': 'writing_snapshot_not_human_submission_freeze',
        'created_at': datetime.now(timezone.utc).isoformat(),
        'source_run': ex['source_run'],
        'user_authorization_quote': '根据这个内容写一个论文用latex编译，只写第一问，写法参考mathmodeling skill里面的内容',
        'human_package_signoff': None,
        'scope': 'Q1 manuscript drafted under explicit user request; no external submission',
        'source_hashes': {str(p.relative_to(REPO)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},
        'main_npz_sha256': convergence['raw_result_sha256'],
        'claims': claims,
        'input_parameters': {'R_m': .02, 'L_m': .25, 'T0_Celsius': 28, 'C0_kg_kg': 2.55,
            'rho_kg_m3': 820, 'cp_J_kg_K': 2600, 'k_W_m_K': .36, 'h_W_m2_K': 25,
            'km_m_s': 8e-7, 'D_formula': '7e-9 * exp(-0.89 / C)'},
        'solver_configuration': summary,
    }
    (PAPER / 'evidence/writing_snapshot.json').write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding='utf-8')
    print('Generated 4 tables, numerical macros, and writing snapshot; 70 table cells checked against full export.')


if __name__ == '__main__':
    main()
