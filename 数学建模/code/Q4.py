# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 4（「无对流」简化模型）

为什么会有 51.00 h 这个答案？
    本脚本给出「去掉 Landau 变换对流项」的简化口径，烘干时长 t_dry ≈ 50.83 h，
    四舍五入 / 向上取整到整小时即 51 h。

    与完整模型（code/p4_model.py + code/p4_run.py，t_dry = 52.40 h）的唯一差别是：
    坐标变换（E = r/R(t)）带来的对流项  (E R'/R)·∂/∂E  被去掉，
    收缩只保留 1/R² 的扩散缩放，而不再把水分向中心「浓缩」。
    该对流项是收缩使水分浓缩、从而拖慢干燥的物理来源（占 t_dry 约 3%），
    漏掉它会让烘干时长从 52.40 h 缩到 50.83 h。

输出：
    results/result4_noconv.xlsx           每隔 60 s、0.1 cm 距离 + 表面列
    results/problem4_tables_noconv.txt    表 6、烘干时长、终止剖面
"""
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_model as m
from p4_model import (C_TARGET, make_ambient, make_R, march)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

# 交付口径
N_DELIVER = 400          # E 域均匀网格节点数（与 p4_verify V3 的 50.833 h 口径一致）
DT = 60.0                # s    输出节拍（问题 4 要求每隔 60 s）
OUT_R_CM = np.round(np.arange(20) * 0.1, 1)   # 固定距离 0, 0.1, ..., 1.9 cm
HEADER_T = "时间\\到药材中心的距离"

# ----------------------------------------------------------------------------
# 关键改动：去掉对流项 (E R'/R)·∂/∂E —— 这就是「51.00」的来源口径
# ----------------------------------------------------------------------------
_orig_conv = m.op_conv


def _zero_conv(g, R, Rp):
    a, b, c = _orig_conv(g, R, Rp)
    return a * 0.0, b * 0.0, c * 0.0


m.op_conv = _zero_conv


def interp_E(g, C_row, E):
    """在 E 域网格上取 C(E)。E ∈ [0,1]，np.interp 线性插值。"""
    return round(float(np.interp(E, g["E"], C_row)), 4)


def build_table_and_sheet(res, t_dry):
    """构造表 6 与 result4 的完整数据（固定距离列 + 表面列）。"""
    g = res["g"]
    times = res["times"]
    C = res["C"]

    # 表 6：每隔 6 h，距离 0/0.5/1.0/1.5 cm + 表面
    table6_t = np.arange(6.0, t_dry / 3600.0 + 1e-9, 6.0) * 3600.0
    table6_r = [0.0, 0.5, 1.0, 1.5]
    table6_rows = []
    for tt in table6_t:
        i = int(round(tt / DT))
        R = res["R_t"][i]
        row = [tt / 3600.0]
        for rcm in table6_r:
            r = rcm * 0.01
            row.append(interp_E(g, C[i], r / R) if r < R else np.nan)
        row.append(round(float(C[i, -1]), 4))
        table6_rows.append(row)
    i_end = len(times) - 1
    R_end = res["R_t"][i_end]
    row_end = [t_dry / 3600.0]
    for rcm in table6_r:
        r = rcm * 0.01
        row_end.append(interp_E(g, C[i_end], r / R_end) if r < R_end else np.nan)
    row_end.append(round(float(C[i_end, -1]), 4))
    table6_rows.append(row_end)

    # result4：每隔 60 s，固定距离 0..1.9 cm + 表面
    out_idx = np.arange(1, len(times))
    sheet = []
    for i in out_idx:
        tt = times[i]
        R = res["R_t"][i]
        row = [int(round(tt))]
        for rcm in OUT_R_CM:
            r = rcm * 0.01
            if r < R:
                row.append(interp_E(g, C[i], r / R))
            else:
                row.append("")
        row.append(round(float(C[i, -1]), 4))
        sheet.append(row)
    return table6_rows, sheet


def bisect_t_dry(g, res, R_fn, Rp_fn, T_inf, C_inf, target=C_TARGET, tol=1.0):
    """在达标跨越区间内二分，把 t_dry 精确到 ~1 s。"""
    times = res["times"]
    C = res["C"]
    cmax = C.max(axis=1)
    idx = int(np.argmax(cmax < target))
    if idx <= 0:
        return float(times[-1]), None
    t_a, t_b = float(times[idx - 1]), float(times[idx])
    T_a, C_a = res["T"][idx - 1].copy(), C[idx - 1].copy()
    while t_b - t_a > tol:
        t_m = 0.5 * (t_a + t_b)
        T_m, C_m = march(g, T_a, C_a, t_a, t_m, R_fn, Rp_fn, T_inf, C_inf)
        if C_m.max() < target:
            t_b = t_m
        else:
            t_a = t_m
    return t_b, None


def main():
    np.set_printoptions(precision=6, suppress=True)
    T_inf, C_inf = make_ambient()
    R_fn, Rp_fn, _, _, _ = make_R()

    print("=" * 70)
    print("问题 4（无对流简化模型，N=%d, dt=%.0f s）" % (N_DELIVER, DT))
    print("=" * 70)

    res = m.solve(T_inf, C_inf, R_fn, Rp_fn, N=N_DELIVER, dt=DT)
    g = res["g"]
    times, C = res["times"], res["C"]

    t_dry, _ = bisect_t_dry(g, res, R_fn, Rp_fn, T_inf, C_inf)
    t_min = np.ceil(t_dry / 60.0) * 60.0
    print("烘干时长：t_dry = %.2f s = %.4f h = %.2f min" % (t_dry, t_dry / 3600, t_dry / 60))
    print("分钟向上取整口径：%.0f min = %.4f h（%.0f s）" % (t_min / 60, t_min / 3600, t_min))
    print("→ 整小时口径 ≈ %d h（即「51.00」的来源）" % int(round(t_dry / 3600.0)))

    table6_rows, sheet = build_table_and_sheet(res, t_dry)

    # result4_noconv.xlsx
    header = [HEADER_T] + [float(x) for x in OUT_R_CM] + ["药材表面"]
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(header)
    for row in sheet:
        ws.append(row)
    xlsx_path = os.path.join(RES, "result4_noconv.xlsx")
    wb.save(xlsx_path)
    print("已保存 %s（%d 行 × %d 列）" % (xlsx_path, len(sheet), len(header)))

    # 表 6 文本 + 终止剖面
    lines = []
    lines.append("问题 4 结果（无对流简化模型：Landau 变换但去掉 (E R'/R)·∂/∂E 对流项）")
    lines.append("")
    lines.append("烘干时长 t_dry = %.2f s = %.4f h（整小时口径 ≈ %d h）"
                 % (t_dry, t_dry / 3600, int(round(t_dry / 3600.0))))
    lines.append("")
    lines.append("表 6  药材烘干过程的水分浓度（kg/kg，干基）")
    lines.append("时间/h\t0 cm\t0.5 cm\t1.0 cm\t1.5 cm\t药材表面")
    for row in table6_rows:
        cells = ["%.4f" % v if not np.isnan(v) else "——" for v in row]
        lines.append("\t".join(cells))
    lines.append("")
    lines.append("注：'——' 表示该固定距离处材料已因收缩离开（r ≥ R(t)）。")
    lines.append("")
    i_end = len(times) - 1
    R_end = res["R_t"][i_end]
    lines.append("终止时刻剖面（t=%.0f s，R=%.4f cm）" % (times[i_end], R_end * 100))
    prof = []
    for rcm in np.arange(0, 2.0 + 1e-9, 0.1):
        r = rcm * 0.01
        if r < R_end:
            prof.append((rcm, interp_E(g, C[i_end], r / R_end)))
    lines.append("r/cm\t" + "\t".join("%.1f" % p[0] for p in prof))
    lines.append("C\t" + "\t".join("%.4f" % p[1] for p in prof))

    txt_path = os.path.join(RES, "problem4_tables_noconv.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("已保存 %s" % txt_path)

    print("\n表 6 预览：")
    for row in table6_rows:
        print("  " + "\t".join("%.4f" % v if not np.isnan(v) else "  ——" for v in row))


if __name__ == "__main__":
    main()
