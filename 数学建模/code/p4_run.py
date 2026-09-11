# -*- coding: utf-8 -*-
"""
问题 4 交付脚本：Landau 变换移动边界模型 → 表 6 + result4.xlsx + 烘干时长。

输出：
    results/result4.xlsx        题目要求的交付文件（每隔 60 s、0.1 cm 距离 + 药材表面列）
    results/problem4_tables.txt 表 6、烘干时长、终止剖面
    results/p4_fields.npz       时空场缓存（供绘图复用）
"""
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_model as m
from p4_model import (T0, C0, C_TARGET, H_MASS, make_ambient, make_R,
                      make_grid, march)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

# 交付口径
N_DELIVER = 800          # E 域均匀网格节点数（物理 Δr ≈ R/N）
DT = 60.0                # s    输出节拍（问题 4 要求每隔 60 s）
OUT_R_CM = np.round(np.arange(20) * 0.1, 1)   # 固定距离 0, 0.1, ..., 1.9 cm（表面列另附）
HEADER_T = "时间\\到药材中心的距离"


def interp_E(g, C_row, E):
    """在 E 域网格上取 C(E)。E ∈ [0,1]，用 np.interp（内部点，二阶误差可忽略）。"""
    return round(float(np.interp(E, g["E"], C_row)), 4)


def build_table_and_sheet(res, t_dry):
    """由求解结果构造 表6 与 result4 的完整数据（固定距离列 + 表面列）。"""
    g = res["g"]
    times = res["times"]          # 每个输出时刻（间隔 DT，含 t=0）
    C = res["C"]                  # (n_t, N+1)

    # ---- 表 6：每隔 6 h，距离 0/0.5/1.0/1.5 cm + 表面 ----
    table6_t = np.arange(6.0, t_dry / 3600.0 + 1e-9, 6.0) * 3600.0
    table6_r = [0.0, 0.5, 1.0, 1.5]                       # cm
    table6_rows = []
    for tt in table6_t:
        i = int(round(tt / DT))                            # 输出时刻索引（tt 是 60 s 整数倍）
        R = res["R_t"][i]                                  # m
        row = [tt / 3600.0]
        for rcm in table6_r:
            r = rcm * 0.01
            row.append(interp_E(g, C[i], r / R) if r < R else np.nan)
        row.append(round(float(C[i, -1]), 4))                       # 表面 C(E=1)
        table6_rows.append(row)
    # 结束行（烘干结束时刻）
    i_end = len(times) - 1
    R_end = res["R_t"][i_end]
    row_end = [t_dry / 3600.0]
    for rcm in table6_r:
        r = rcm * 0.01
        row_end.append(interp_E(g, C[i_end], r / R_end) if r < R_end else np.nan)
    row_end.append(round(float(C[i_end, -1]), 4))
    table6_rows.append(row_end)

    # ---- result4：每隔 60 s，固定距离 0..1.9 cm + 表面 ----
    # 时间列自 60 s 起每 60 s（不含 t=0，与模板一致）
    out_idx = np.arange(1, len(times))                     # 跳过 t=0
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
                row.append("")                             # 材料已收缩离开该位置
        row.append(round(float(C[i, -1]), 4))                        # 表面列
        sheet.append(row)
    return table6_rows, sheet


def bisect_t_dry(g, res, R_fn, Rp_fn, T_inf, C_inf, target=C_TARGET, tol=1.0):
    """在达标跨越区间内二分，把 t_dry 精确到 ~1 s。"""
    times = res["times"]
    C = res["C"]
    # 首个 max_r C < target 的输出时刻（C 随时间递减，cmax 单调下降）
    cmax = C.max(axis=1)
    idx = int(np.argmax(cmax < target))
    if idx <= 0:
        return float(times[-1]), None
    t_a, t_b = float(times[idx - 1]), float(times[idx])
    T_a, C_a = res["T"][idx - 1].copy(), C[idx - 1].copy()
    hist = []
    while t_b - t_a > tol:
        t_m = 0.5 * (t_a + t_b)
        T_m, C_m = march(g, T_a, C_a, t_a, t_m, R_fn, Rp_fn, T_inf, C_inf)
        hist.append((t_m, float(C_m.max())))
        if C_m.max() < target:
            t_b = t_m
        else:
            t_a = t_m
    return t_b, hist


def main():
    np.set_printoptions(precision=6, suppress=True)
    T_inf, C_inf = make_ambient()
    R_fn, Rp_fn, tR, R_data, Rp_data = make_R()

    print("=" * 70)
    print("问题 4：Landau 变换移动边界模型（N=%d, dt=%.0f s）" % (N_DELIVER, DT))
    print("=" * 70)

    res = m.solve(T_inf, C_inf, R_fn, Rp_fn, N=N_DELIVER, dt=DT)
    g = res["g"]
    times, C = res["times"], res["C"]

    t_dry, hist = bisect_t_dry(g, res, R_fn, Rp_fn, T_inf, C_inf)
    t_min = np.ceil(t_dry / 60.0) * 60.0
    print("烘干时长：t_dry = %.2f s = %.4f h = %.2f min" % (t_dry, t_dry / 3600, t_dry / 60))
    print("分钟向上取整口径：%.0f min = %.4f h（%.0f s）" % (t_min / 60, t_min / 3600, t_min))

    table6_rows, sheet = build_table_and_sheet(res, t_dry)

    # ---- 保存 result4.xlsx ----
    header = [HEADER_T] + [float(x) for x in OUT_R_CM] + ["药材表面"]
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(header)
    for row in sheet:
        ws.append(row)
    xlsx_path = os.path.join(RES, "result4.xlsx")
    wb.save(xlsx_path)
    print("已保存 %s（%d 行 × %d 列）" % (xlsx_path, len(sheet), len(header)))

    # ---- 表 6 文本 + 终止剖面 ----
    lines = []
    lines.append("问题 4 结果（Landau 变换移动边界模型，附录 4 变参数物性 + 附件 2 收缩半径）")
    lines.append("")
    lines.append("烘干时长 t_dry = %.2f s = %.4f h（分钟向上取整 %.0f min = %.4f h）"
                 % (t_dry, t_dry / 3600, t_min / 60, t_min / 3600))
    lines.append("")
    lines.append("表 6  药材烘干过程的水分浓度（kg/kg，干基）")
    lines.append("时间/h\t0 cm\t0.5 cm\t1.0 cm\t1.5 cm\t药材表面")
    for row in table6_rows:
        cells = ["%.4f" % v if not np.isnan(v) else "——" for v in row]
        lines.append("\t".join(cells))
    lines.append("")
    lines.append("注：'——' 表示该固定距离处材料已因收缩离开（r ≥ R(t)，无药材内部水分）。")
    lines.append("")
    # 终止时刻（分钟口径）径向剖面，降采样到 0.1 cm
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

    txt_path = os.path.join(RES, "problem4_tables.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("已保存 %s" % txt_path)

    # ---- 缓存场 ----
    np.savez_compressed(os.path.join(RES, "p4_fields.npz"),
                        times=times, E=g["E"], T=res["T"], C=C,
                        R_t=res["R_t"], Rp_t=res["Rp_t"],
                        tR=tR, R_data=R_data, Rp_data=Rp_data)
    print("已保存 results/p4_fields.npz")

    print("\n表 6 预览：")
    for row in table6_rows:
        print("  " + "\t".join("%.4f" % v if not np.isnan(v) else "  ——" for v in row))


if __name__ == "__main__":
    main()
