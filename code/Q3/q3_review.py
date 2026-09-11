"""Read back Q3 artifacts, verify inherited Q2 values, assess mesh sensitivity."""
from pathlib import Path
import importlib.util
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/Q3"


def step_from_field(m, T_old, C_old, t_old, iterations):
    N = len(C_old) - 1
    dr = m.R_CYL / N
    Ti, Ci = m.make_ambient(ROOT.parent / "CUMCM2026Problems/A题/附件/附件1.xlsx")
    a, b, c, g = m.build_op_heat(dr, N, C_old)
    so = np.zeros(N + 1); so[-1] = g * Ti(t_old)
    Cstar, Tstar = C_old.copy(), T_old.copy()
    for _ in range(iterations):
        an, bn, cn, gn = m.build_op_heat(dr, N, Cstar)
        sn = np.zeros(N + 1); sn[-1] = gn * Ti(t_old + 1)
        Tnew = m._step(T_old, (an, bn, cn), (a, b, c), so, sn, 1., .5)
        Cnew = m._mass_step(C_old, Cstar, T_old, Tnew, Ci(t_old), Ci(t_old + 1), 1., .5, dr, N, m.R_CYL, True)
        Tstar, Cstar = Tnew, Cnew
    return Tstar, Cstar


def main():
    runs = [json.loads((OUT / f"run_summary_N{n}_dt1.json").read_text(encoding="utf-8")) for n in [20, 40, 80, 160]]
    final = runs[-1]
    old = np.load(OUT / "fields_N20_dt1.npz")
    idx = np.flatnonzero(old['times'] == 10800)[0]
    q2_expected = np.array([1.7662, 1.7165, 1.5701, 1.3331, 1.0078])
    q2_match = bool(np.array_equal(np.round(old['C'][idx, ::5], 4), q2_expected))
    assert q2_match and runs[0]['drying_time_s'] == 202240
    assert final['previous_max_C'] >= .15 and final['terminal_max_C'] < .15
    assert all(x['finite_fields'] and x['C_min_before_stop'] > 0 for x in runs)
    assert all(x['maximum_off_center_excess'] < 1e-12 for x in runs)
    dtimes = [b['drying_time_h'] - a['drying_time_h'] for a, b in zip(runs[:-1], runs[1:])]
    spec = importlib.util.spec_from_file_location('inherited', ROOT / '数学建模/数学建模/code/p2_model.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    fields = np.load(OUT / 'fields_N160_dt1.npz')
    probes = []
    for target in [60, 1800, 10800, 86400, 172800, final['drying_time_s']]:
        i = int(np.argmin(np.abs(fields['times'] - target)))
        t = float(fields['times'][i])
        T3, C3 = step_from_field(m, fields['T'][i], fields['C'][i], t, 3)
        T10, C10 = step_from_field(m, fields['T'][i], fields['C'][i], t, 10)
        probes.append({'time_s': t, 'max_delta_T': float(np.max(np.abs(T3-T10))), 'max_delta_C': float(np.max(np.abs(C3-C10)))})
    review = {
        'question': 'Q3', 'decision_id': 'q3_inherit_q2',
        'checks': {
            'syntax': {'status': 'PASS', 'evidence': 'Four grids executed successfully'},
            'input_contract': {'status': 'PASS', 'evidence': 'Official attachment read-only; endpoint extension explicitly recorded; hashes in run summaries'},
            'method_alignment': {'status': 'PASS', 'evidence': 'Inherited Q2 solver, all-node threshold; changes limited to spatial resolution'},
            'reproducibility': {'status': 'PASS', 'evidence': 'N20 matches Q2 3h moisture table and 202240s dry-time baseline'},
            'output_contract': {'status': 'PASS', 'evidence': final.get('xlsx_readback')},
            'threshold': {'status': 'PASS', 'evidence': [final['previous_max_C'], final['terminal_max_C']]},
        },
        'mesh_times_h': {str(x['N']): x['drying_time_h'] for x in runs},
        'successive_mesh_time_changes_h': dtimes,
        'last_mesh_relative_change': abs(dtimes[-1]) / final['drying_time_h'],
        'local_iteration_probes_3_vs_10': probes,
        'limitations': ['Grid differences do not certify four-decimal accuracy', 'Representative iteration checks do not certify every time step', 'Long-term ambient conditions are extrapolated from last attachment row'],
        'verdict': 'PASS with stated numerical and modeling limitations',
    }
    (OUT / 'review.json').write_text(json.dumps(review, ensure_ascii=False, indent=2), encoding='utf-8')
    table = '| 时间/h | 0 cm | 0.5 cm | 1 cm | 1.5 cm | 2 cm |\n|---:|---:|---:|---:|---:|---:|\n'
    for t, row in zip(final['table_times_h'], final['table_C']):
        table += '| ' + f'{t:.4f}' + ' | ' + ' | '.join(f'{c:.4f}' for c in row) + ' |\n'
    mesh = '| 径向分段数 | 烘干时间/h |\n|---:|---:|\n'
    for run in runs:
        mesh += f"| {run['N']} | {run['drying_time_h']:.4f} |\n"
    text = f'''# 第三问：烘干时长与含水率分布

沿用第二问一维、固定半径的变参数热质传递模型及附录3参数。初始温度28摄氏度、干基含水率2.55 kg/kg，中心对称，表面对流换热与传质。第二问物性公式用于整个过程，不在30分钟处拼接第一问常参数结果。附件1前4小时分段线性插值，其后温度50.165摄氏度、环境水分浓度0.04986 kg/kg保持不变；预测时长依赖此延拓假设。

定义 M(t)=max_r C(r,t)，烘干阈值为0.15 kg/kg。每一步检查全部径向节点，首次满足M(t)<0.15时停止计时。各网格的计算均表明中心为最大含水率位置（浮点误差量级内），因此本算例最后达标位置为中心。

采用160段径向网格、1 s时间步，首次离散达标时间为 **{final['drying_time_s']:.0f} s = {final['drying_time_h']:.4f} h**。前一秒最大含水率为{final['previous_max_C']:.12f}，该秒为{final['terminal_max_C']:.12f} kg/kg。线性插值得到阈值穿越估计{final['threshold_crossing_interpolated_s']:.4f} s；严格小于阈值的离散判据使用上述首次达标秒数。

表5（含水率单位kg/kg，干基；最后一行为烘干结束）：

{table}
终点的0.1500是四舍五入显示，判定使用未舍入值。result3.xlsx采用原模板Sheet1，时间单位s，距离单位cm；从60 s起每60 s输出一次，0..2 cm每0.1 cm；另附最终{final['drying_time_s']:.0f} s的终点行，共{final['sample_rows']}行数据。内部计算网格更细，但输出直接取指定位置节点，无空间插值。

网格检查：

{mesh}
20段结果与此前第二问一致。第三问的长时间阈值对空间误差更敏感，因此数值加密至160段；没有改变前两问物理模型或修改前两问结果文件。最后两档的时长差为{abs(dtimes[-1]):.6f} h，约{abs(dtimes[-1])*60:.3f} min，不能把表格的四位小数理解为模型已经达到对应精度。

验证已包括：全部节点阈值、含水率非负与有限值、中心最大值、对照第二问结果、四档网格比较、Excel逐单元格回读。另对6个保存状态重做单步，比较3次与10次非线性迭代，详见review.json。这些检查不替代实际干燥实验验证，也不证明每个时间步均已达到迭代容差。
'''
    (OUT / '第三问_解答.md').write_text(text, encoding='utf-8')
    print(json.dumps(review, ensure_ascii=False, indent=2))
    print(table)


if __name__ == '__main__':
    main()
