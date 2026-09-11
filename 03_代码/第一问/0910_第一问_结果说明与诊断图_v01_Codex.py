"""2026-09-10 Q1: render diagnostics and write the evidence-based report."""
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / '04_结果/第一问'
FIG = REPO / '05_图片/第一问'
CODE = REPO / '03_代码/第一问'
REPORT = OUT / '0910_第一问_结果与模型说明_v01_Codex.md'


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def table(times, radii, values):
    lines = ['| 时间 / s | ' + ' | '.join(f'{r:g} cm' for r in radii) + ' |',
             '|---:|' + '---:|' * len(radii)]
    lines.extend('| ' + str(t) + ' | ' + ' | '.join(f'{v:.4f}' for v in row) + ' |'
                 for t, row in zip(times, values))
    return '\n'.join(lines)


def save_figure(fig, name):
    path = FIG / name
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)
    return path


def figures(main, baseline, times):
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'],
        'axes.unicode_minus': False, 'font.size': 11, 'axes.titlesize': 13,
        'axes.labelsize': 11, 'legend.fontsize': 10, 'figure.titlesize': 15,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.labelcolor': '#333333', 'text.color': '#333333',
        'xtick.color': '#555555', 'ytick.color': '#555555',
    })
    colors = ['#1A6FC4', '#E28E2C', '#7B5FD6', '#33B5A5', '#D9544D', '#666666', '#283B4D']
    styles = ['-', '--', '-.', ':', '-', '--', '-.']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, field, ylabel, title in zip(
            axes, ['snapshots_t', 'snapshots_c'],
            ['温度 T / °C', '干基含水率 C / (kg/kg)'],
            ['(a) 温度径向分布', '(b) 含水率径向分布']):
        for i, (t, color, style) in enumerate(zip(times, colors, styles)):
            ax.plot(main['r'] * 100, main[field][i, 0], color=color,
                    ls=style, lw=1.9, label=f'{t} s')
        ax.set(xlabel='距中心径向距离 r / cm', ylabel=ylabel, title=title, xlim=(0, 2))
        ax.set_xticks(np.arange(0, 2.01, .5))
        ax.grid(color='#D8D8D8', lw=.6, alpha=.65)
    fig.suptitle('第一问：中截面 z = 0 的径向分布', y=.98)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .93), ncol=7,
               frameon=False, columnspacing=1.25)
    fig.subplots_adjust(top=.77, bottom=.14, left=.08, right=.98, wspace=.25)
    radial_path = save_figure(fig, '0910_第一问_中截面径向分布_v01_Codex.png')

    fig, axes = plt.subplots(2, 1, figsize=(11, 6.4), layout='constrained')
    for ax, field, cmap, title, unit in zip(
            axes, ['snapshots_t', 'snapshots_c'], ['inferno', 'viridis'],
            ['(a) 温度场', '(b) 干基含水率场'], ['T / °C', 'C / (kg/kg)']):
        values = main[field][-1].T
        image = ax.pcolormesh(main['z'] * 100, main['r'] * 100, values,
                              shading='gouraud', cmap=cmap, rasterized=True)
        ax.set(xlabel='距中截面的轴向距离 z / cm', ylabel='半径 r / cm', title=title,
               xlim=(0, 12.5), ylim=(0, 2), aspect='equal')
        ax.set_xticks(np.arange(0, 12.51, 2.5))
        ax.set_yticks([0, 1, 2])
        fig.colorbar(image, ax=ax, pad=.025, shrink=.85, label=unit)
    fig.suptitle('第一问：1800 s 的二维轴对称半域（诊断）')
    field_path = save_figure(fig, '0910_第一问_1800s二维场诊断_v01_Codex.png')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), layout='constrained')
    for ax, column, ylabel in zip(axes, [1, 2], ['体积平均温度 / °C', '体积平均含水率 / (kg/kg)']):
        ax.plot(main['times'] / 60, main['balances'][:, column], color='#1A6FC4',
                lw=2, label='二维有限圆柱')
        ax.plot(baseline['times'] / 60, baseline['balances'][:, column], color='#767676',
                lw=1.8, ls='--', label='一维径向对照')
        ax.set(xlabel='时间 / min', ylabel=ylabel, xlim=(0, 30))
        ax.grid(color='#D8D8D8', lw=.6, alpha=.65)
        ax.legend(frameon=False)
    fig.suptitle('第一问：端面传递对整体平均状态的影响（模型内对照）')
    means_path = save_figure(fig, '0910_第一问_二维与一维平均量对照_v01_Codex.png')
    return [radial_path, field_path, means_path]


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    export = read_json(OUT / 'q1_export.json')
    summary = read_json(OUT / 'q1_main_r16_summary.json')
    metrics = read_json(OUT / 'experiments/round1/metrics/q1_convergence.json')
    run = np.load(OUT / 'q1_main_r16.npz')
    baseline = np.load(OUT / 'q1_baseline_r16.npz')
    diagnostic_reference = np.load(OUT / 'q1_fine.npz')
    source_hash = hashlib.sha256((OUT / 'q1_main_r16.npz').read_bytes()).hexdigest()
    assert source_hash == metrics['raw_result_sha256']
    assert export['source_run'] == summary['label'] == 'q1_main_r16'
    for field in ('temperature', 'moisture'):
        assert np.array_equal(np.asarray(export[field]), np.round(run[field], 4))
    assert np.allclose(run['r'][::4], diagnostic_reference['r'], rtol=0, atol=1e-14)
    assert np.allclose(run['z'], diagnostic_reference['z'][::4], rtol=0, atol=1e-14)
    domain_delta = {
        field: float(np.max(np.abs(run[field][:, :, ::4]
                                  - diagnostic_reference[field][:, ::4, :])))
        for field in ('snapshots_t', 'snapshots_c')
    }
    figure_paths = figures(run, baseline, export['table_times'])
    comparison_rows = []
    comparison_labels = {
        'main_vs_1d_identical_radial_time': '二维主解与同径向网格、同步长的一维对照',
        'time_refinement_1d_dt_half': '一维对照时间步长减半',
        'space_refinement_1d_radial_double': '同细时间步长下一维径向网格加密一倍',
        'main_vs_finer_1d_combined': '二维主解与更细一维参考',
    }
    for key, label in comparison_labels.items():
        comparison_rows.append(f"| {label} | {metrics[key]['temperature']['all_outputs_max_abs']:.6g} | "
                               f"{metrics[key]['moisture']['all_outputs_max_abs']:.6g} |")
    template = r'''# 第一问：模型、结果与数值验证

采用有限圆柱二维轴对称守恒有限体积模型，按用户确认的中截面 \(z=0\) 输出径向结果，并附同参数一维径向对照。范围为预热阶段 \(0\le t\le1800\ \mathrm{s}\)。本说明是当前计算成果，不是论文数字冻结或外部提交记录。

## 1. 输入、假设与物性

药材长度 \(L=0.25\ \mathrm m\)、半径 \(R=0.02\ \mathrm m\)，初始温度 \(28^\circ\mathrm C\)、干基含水率 \(2.55\ \mathrm{kg/kg}\)。利用中截面对称性，计算域为 \(0\le r\le R\)、\(0\le z\le L/2\)。输出的“到药材中心距离”按中截面径向距离解释，不能将该表当作任意轴向位置的结果。

| 参数 | 数值及单位 |
|---|---|
| 密度 \(\rho\) | 820 kg/m³ |
| 比热 \(c_p\) | 2600 J/(kg·K) |
| 热导率 \(k\) | 0.36 W/(m·K) |
| 换热系数 \(h\) | 25 W/(m²·K) |
| 传质系数 \(k_m\) | \(8\times10^{-7}\) m/s |
| 扩散系数 \(D(C)\) | \(7\times10^{-9}\exp(-0.89/C)\) m²/s |

环境边界使用附件 1 的 0—1800 s 共 31 个实测时间点，作分段线性插值。该区间无需外推，且没有在1800 s之前人为切换恒温。初始环境为 28 °C、0.01963 kg/kg，1800 s 时为 41.513 °C、0.03307 kg/kg。审计确认附件全表 241 条记录无缺失或重复，数值未作插补、剔除或平滑。

本问采用各向同性有效扩散和有效热物性，忽略体积收缩。侧面与两端面都暴露于同一环境，端面传热、传质系数取侧面系数。空气数据与药材含水率的 kg/kg 分母并未在题面建立物理转换关系，因此 \(C_{\rm eq}=C_\infty\) 是明确的经验边界闭合。本次未加入缺少完整参数的吸附等温线、气膜蒸汽压和潜热机制，不能解释为已经证实潜热影响为零。

## 2. 控制方程与边界

温度与水分分别满足

\[
\rho c_p\frac{\partial T}{\partial t}
=\frac1r\frac{\partial}{\partial r}\left(rk\frac{\partial T}{\partial r}\right)
+\frac{\partial}{\partial z}\left(k\frac{\partial T}{\partial z}\right),
\]

\[
\frac{\partial C}{\partial t}
=\frac1r\frac{\partial}{\partial r}\left(rD(C)\frac{\partial C}{\partial r}\right)
+\frac{\partial}{\partial z}\left(D(C)\frac{\partial C}{\partial z}\right).
\]

第一问的热物性均为常数，扩散系数仅依赖含水率，因而两场可以分别求解；本问没有使用第二问才出现的双向温湿物性耦合。水分方程保持散度形式，不把它替换为 \(D(C)\Delta C\)，后者遗漏 \(D'(C)|\nabla C|^2\)。

中心轴 \(r=0\) 和中截面 \(z=0\) 为对称边界：\(T_r=C_r=0\)、\(T_z=C_z=0\)。在侧面 \(r=R\) 与端面 \(z=L/2\)，以外法向 \(\boldsymbol n\) 定义

\[
-k\nabla T\cdot\boldsymbol n=h(T_s-T_\infty),\qquad
-D(C_s)\nabla C\cdot\boldsymbol n=k_m(C_s-C_{\rm eq}).
\]

表面值由边界控制体求解，不强行设为环境值。端面与侧面的交界控制体同时包含两张外表面的通量。初始含水率均匀而表面立刻发生传质，初值与 Robin 条件在 \(t=0\) 角点不相容；保留题设初值，通过细时间步解析启动层。

## 3. 离散、收敛与守恒

采用轴对称节点中心有限体积，径向 {{NR}} 节点、半轴向 {{NZ}} 节点，共 {{NODES}} 个节点。径向最大间距 {{DRMAX}} mm，最小间距 {{DRMIN}} μm；最后 1 mm 区域加密。轴向网格向端面加密，最小间距 {{DZMIN}} μm。中心通量通过零面积控制面处理，内部相邻控制体使用同一通量且符号相反。

面扩散率取相邻节点扩散率的调和平均。真实边界节点的外通量直接由 Robin 条件计算。温度容量为 \(\rho c_pV_i\)，水分容量为 \(V_i\)。这对应空间均匀、固定干物质密度下约去密度的含水率方程。

主时间步长为 {{DT}} s，最初 {{FIRSTDT}} s 区间采用两个隐式 Euler 半步；随后按启动层分级增大步长，并在每个整数秒输出。其余时间步采用非线性 Crank–Nicolson，当前时间层内更新 \(D(C)\)。Picard 增量阈值为 \(2\times10^{-10}\)，离散方程残差阈值为 \(2\times10^{-9}\)。本次共 {{STEPS}} 个内部时间步，最大 {{MAXITER}} 次 Picard 迭代，记录的最大非线性残差 {{NONLINEAR}} kg/kg。

初期 CN 步长依次受以下上限约束：到0.001 s为 \(\Delta t/5000\)，到0.01 s为 \(\Delta t/500\)，到0.1 s为 \(\Delta t/50\)，到1 s为 \(\Delta t/25\)，到10 s为 \(\Delta t/10\)，到60 s为 \(\Delta t/4\)，之后为 \(\Delta t\)。这属于预设分级启动与加密验证，本次没有实现误差估计驱动的通用自适应步长，不能将其写成已经实现的功能。

令 \(I_C=\int_\Omega C\,dV\)，则检查

\[
I_C(t)-I_C(0)+\int_0^t\int_{\partial\Omega}
k_m(C_s-C_{\rm eq})\,dA\,ds=0.
\]

\(I_C\) 是含水率体积分，不是以 kg 计的水质量；若要转为水质量，还需乘以一致定义的干物质体积密度。热量平衡检查使用 \(\int_\Omega\rho c_p(T-T_0)dV\) 与表面净传热积分。边界时间积分使用与时间推进相同的 Euler/CN 权重。除以域体积及相应容量后的最大累计平衡残差分别为 {{HEATBAL}} °C、{{WATERBAL}} kg/kg。半域计算中共同的 \(2\pi\) 因子在归一化后约去。

## 4. 题目要求的两张结果表

表1：30分钟内中截面温度，单位 °C。

{{TEMPERATURE_TABLE}}

表2：30分钟内中截面干基含水率，单位 kg/kg。

{{MOISTURE_TABLE}}

1800 s 时中截面中心与侧面温度分别为 {{CENTER_T}} °C、{{SURFACE_T}} °C，含水率分别为 {{CENTER_C}}、{{SURFACE_C}} kg/kg。药材外层升温较快、失水较多，中心水分在此时间范围内变化很小。该含水率分布说明第一问结束并不等于达到第三问的最终干燥要求。

## 5. 一维对照与数值误差范围

以下为全部 1800 个输出时刻、21 个径向位置的最大绝对差。二维输出均取 \(z=0\)，一维模型使用相同物性、初值、环境边界和径向有限体积公式。

| 比较内容 | 温度最大差 / °C | 含水率最大差 / (kg/kg) |
|---|---:|---:|
{{COMPARISON_ROWS}}

一维同网格基线的径向细化参数为16、主步长为0.5 s，时间加密为0.25 s；空间加密参考取径向细化参数32及相同的0.25 s步长。上述比较支持当前中截面输出的数值收敛性。四位小数是题目要求的显示格式，网格差不是严格的解析误差上界，也不保证每个临近舍入分界的单元最后一位都不变。不存在内部温度/含水率实测观测，因此这里不报告实验预测精度或统计置信区间。

虽然中截面差很小，端面仍影响全体积平均量。1800 s 的二维体积平均温度和含水率为 {{MEAN_T_2D}} °C、{{MEAN_C_2D}} kg/kg；一维径向对照为 {{MEAN_T_1D}} °C、{{MEAN_C_1D}} kg/kg。两者之差分别为 {{DELTA_T}} °C、{{DELTA_C}} kg/kg。不能由中截面的近似一致推断端面可在全域计算中忽略。

二维图件只作为分布诊断。将当前二维主解与此前另一套二维网格 q1_fine 在7个快照的共同节点比较，温度最大差为 {{DOMAIN_T_DIFF}} °C，含水率最大差为 {{DOMAIN_C_DIFF}} kg/kg；两组还使用不同启动步安排，这是一项差异诊断，并非纯轴向收敛证明。全域端部/角点未认证到四位小数精度，本次精度检验的主要交付范围是官方中截面径向输出。

## 6. 诊断图

![中截面径向温湿分布]({{FIG_RADIAL}})

图1展示7个指定时刻的中截面径向剖面，源自未舍入的二维快照。

![1800秒二维场]({{FIG_FIELD}})

图2展示1800 s时二维轴对称半域；左侧是对称中截面，右侧是暴露端面，上边界是圆柱侧面。色彩表示计算状态；细网格之间的图像平滑不构成额外数值精度。

![二维和一维平均量对照]({{FIG_MEANS}})

图3比较有限圆柱与一维径向模型的平均状态，说明端面影响集中于靠近端部的区域。以上图件属于诊断成果，尚未作为冻结论文图发布。

## 7. 文件与重现

- `result1.xlsx`：两个工作表，每表1800行时间、21列半径结果，共75600个数据值，显示四位小数。
- `q1_main_r16.npz`、`q1_main_r16_summary.json`：主解原始精度数值及参数记录。
- `q1_baseline_r16.npz`、`q1_reference_r16.npz`、`q1_spacecheck_r32.npz`：对照与加密证据。
- `q1_export.json`：由主解生成的官方矩阵与两张论文表的统一数值来源。
- `experiments/round1/metrics/q1_convergence.json`、`experiments/round1/run_summary.json`：比较指标与运行证据。
- `../../02_数据/data_profile.json`：原始来源、输入哈希与清洗记录。

以下在 private 仓库根目录用 PowerShell 执行。`python` 应指向安装 numpy、scipy、matplotlib 的环境；本机为 `D:/python3.10/python.exe`。本次求解环境为 Python {{PYTHON}}、NumPy {{NUMPY}}、SciPy {{SCIPY}}。表格导出另需 bundled Node.js 及已准备好的 artifact-tool 运行目录，详见导出脚本文件头。

```powershell
python "03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py" --dimension 2 --refinement 16 --axial-refinement 1 --dt 0.5 --label q1_main_r16
python "03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py" --dimension 1 --refinement 16 --dt 0.5 --label q1_baseline_r16
python "03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py" --dimension 1 --refinement 16 --dt 0.25 --label q1_reference_r16
python "03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py" --dimension 1 --refinement 32 --dt 0.25 --label q1_spacecheck_r32
python "03_代码/第一问/0910_第一问_热湿模型_v01_Codex.py" --dimension 2 --refinement 4 --axial-refinement 4 --dt 0.5 --legacy-startup --label q1_fine
python "03_代码/第一问/0910_第一问_结果汇总_v01_Codex.py"
python "03_代码/第一问/0910_第一问_结果说明与诊断图_v01_Codex.py"
node "03_代码/第一问/0910_第一问_结果表导出_v01_Codex.mjs"
```

主解 NPZ SHA-256：`{{HASH}}`。本说明与两张表由脚本读取已保存结果生成，不独立手录计算数值。四位显示、数学误差检验、模型物理假设是不同层面的信息，使用结果时应同时保留。
'''
    fields = {
        'NR': summary['grid_shape_z_r'][1], 'NZ': summary['grid_shape_z_r'][0],
        'NODES': int(np.prod(summary['grid_shape_z_r'])),
        'DRMAX': f"{np.diff(run['r']).max()*1000:.6g}",
        'DRMIN': f"{summary['radial_min_spacing_m']*1e6:.6g}",
        'DZMIN': f"{summary['axial_min_spacing_m']*1e6:.6g}",
        'DT': summary['dt_s'], 'FIRSTDT': summary['first_dt_s'], 'STEPS': summary['steps'],
        'MAXITER': summary['max_picard_iterations'],
        'NONLINEAR': f"{summary['max_nonlinear_residual_C']:.4g}",
        'HEATBAL': f"{summary['max_heat_balance_residual_mean_Celsius']:.4g}",
        'WATERBAL': f"{summary['max_water_balance_residual_mean_kg_kg']:.4g}",
        'TEMPERATURE_TABLE': table(export['table_times'], export['table_radii_cm'], export['table_temperature']),
        'MOISTURE_TABLE': table(export['table_times'], export['table_radii_cm'], export['table_moisture']),
        'CENTER_T': f"{summary['final_center_T_C'][0]:.4f}",
        'CENTER_C': f"{summary['final_center_T_C'][1]:.4f}",
        'SURFACE_T': f"{summary['final_midplane_surface_T_C'][0]:.4f}",
        'SURFACE_C': f"{summary['final_midplane_surface_T_C'][1]:.4f}",
        'COMPARISON_ROWS': '\n'.join(comparison_rows),
        'MEAN_T_2D': f"{metrics['final_volume_means_2d'][0]:.4f}",
        'MEAN_C_2D': f"{metrics['final_volume_means_2d'][1]:.4f}",
        'MEAN_T_1D': f"{metrics['final_radial_means_1d'][0]:.4f}",
        'MEAN_C_1D': f"{metrics['final_radial_means_1d'][1]:.4f}",
        'DELTA_T': f"{metrics['final_volume_mean_2d_minus_1d'][0]:+.4f}",
        'DELTA_C': f"{metrics['final_volume_mean_2d_minus_1d'][1]:+.4f}",
        'DOMAIN_T_DIFF': f"{domain_delta['snapshots_t']:.6g}",
        'DOMAIN_C_DIFF': f"{domain_delta['snapshots_c']:.6g}",
        'FIG_RADIAL': figure_paths[0].as_posix(), 'FIG_FIELD': figure_paths[1].as_posix(),
        'FIG_MEANS': figure_paths[2].as_posix(), 'HASH': source_hash,
        'PYTHON': summary['environment']['python'], 'NUMPY': summary['environment']['numpy'],
        'SCIPY': summary['environment']['scipy'],
    }
    for key, value in fields.items():
        template = template.replace('{{' + key + '}}', str(value))
    assert '{{' not in template
    REPORT.write_text(template, encoding='utf-8')
    print(json.dumps({'report': str(REPORT), 'figures': [str(p) for p in figure_paths],
                      'domain_diagnostic_max_abs': domain_delta}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
