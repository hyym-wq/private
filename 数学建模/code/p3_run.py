# -*- coding: utf-8 -*-
"""
问题 3 求解与结果输出（交付口径）。

交付口径：在细网格上求解（N = 640 均匀网格，Δr = 3.125e-5 m = 0.003125 cm），
按题目要求降采样输出：
  * result3.xlsx —— 每隔 60 s、到中心距离每隔 0.1 cm 的水分浓度；
  * 表 5 —— 每隔 6 h，距离 0/0.5/1/1.5/2 cm + 烘干结束时间行。

降采样用单调三次插值（PCHIP）：分级/细网格节点 → 0.1 cm 输出网格。选用 PCHIP
而非线性插值，是因为干燥末期剖面在表面附近极陡（0.003 cm 内 C 变化 0.03），
线性插值会把表面值系统性抬高；PCHIP 保单调、无过冲，且对细网格的误差是 O(Δr⁴)。

同时给出分级网格（N = 60, q = 2）的独立复算，作为交付结论的交叉校核。
"""
import os
import sys
import time
# 控制台按 UTF-8 输出（Windows GBK 终端下避免 UnicodeEncodeError 与乱码）
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_model as pm
import p3_model as p3

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

GRID_MAIN = dict(N=640, kind="uniform", q=1.0)      # 交付口径（细网格）
GRID_XCHK = dict(N=60, kind="graded", q=2.0)        # 交叉校核（分级网格）
T_MAX = 345600.0
D_TAB = [0.0, 0.5, 1.0, 1.5, 2.0]                   # 表 5 的距离列（cm）
R4 = "{:.4f}"
R6 = "{:.6f}"


def run_case(tag, N, kind, q):
    t0 = time.time()
    res = p3.solve_problem3(N=N, kind=kind, q=q, t_max=T_MAX, verbose=True)
    res["cpu"] = time.time() - t0
    res["tag"] = tag
    print("[%s] N=%d %s q=%.2f  Δr∈[%.5f, %.4f] cm  t_dry=%.4f h  CPU=%.1f s\n"
          % (tag, N, kind, q, np.diff(res["r"]).min() * 100,
             np.diff(res["r"]).max() * 100, res["t_dry"] / 3600, res["cpu"]))
    return res


def qs_table(times, C, w, t_dry, C_inf_late, starts_h):
    """在若干起始时刻用准稳态解析式外推终止时间，与二分结果比较。

    QS-冻结形状  —— 纯解析（只用当前时刻的形状，不依赖历史轨道）；
    QS-漂移修正  —— 用该时刻之前的轨道拟合形状漂移，再积分一维标量方程。
    """
    W = float(w.sum())
    rows = []
    for th in starts_h:
        t0 = th * 3600.0
        i = int(np.searchsorted(times, t0))
        if i >= len(times):
            continue
        qs = p3.qs_state(times[i], C[i], w, C_inf=C_inf_late)
        t_frz, why = p3.qs_extrapolate(qs)
        pb, pn, _, _, _ = p3.qs_shape_fit(times, C, w, t0)
        t_drf, _ = p3.qs_extrapolate_drift(C[i, 0], times[i], pb, pn, W,
                                           target=p3.C_TARGET, C_inf=C_inf_late)
        rows.append(dict(t0=times[i], C0=float(C[i][0]), phibar=qs["phibar"],
                         phiN=qs["phiN"], C_eq=qs["C_eq"],
                         t_frz=t_frz, err_frz=(t_frz - t_dry) / 60.0, why=why,
                         t_drf=t_drf, err_drf=(t_drf - t_dry) / 60.0))
    return rows


def main():
    lines = []

    # ---------------- 交付口径：细网格 ----------------
    res = run_case("交付", **GRID_MAIN)
    times, C, T, r = res["times"], res["C"], res["T"], res["r"]
    t_dry, t_min = res["t_dry"], res["t_dry_min"]
    C_inf_late = float(res["solver"].C_inf(p3.T_ENV_END))

    # 降采样：细网格 -> 0.1 cm 输出网格
    Cd = p3.downsample(r, C, p3.OUT_R)
    Td = p3.downsample(r, T, p3.OUT_R)
    Ce_d = p3.downsample(r, res["C_end"], p3.OUT_R)[0]
    r_cm = p3.OUT_R * 100.0
    d_idx = [int(round(d / 0.1)) for d in D_TAB]

    # 口径核查：max_r C 是否恒在中心
    arg = np.argmax(C, axis=1)
    n_off = int(np.sum(arg != 0))
    lines.append("表 5  药材烘干过程的水分浓度（kg/kg，干基）")
    lines.append("时间/h\t" + "\t".join(f"{d:g} cm" for d in D_TAB))

    th = 6.0
    while th * 3600.0 < t_dry - 1e-9:
        i = int(np.searchsorted(times, th * 3600.0))
        lines.append(f"{th:g}\t" + "\t".join(R4.format(Cd[i][di]) for di in d_idx))
        th += 6.0
    lines.append(f"烘干结束时间（{t_min/3600:.4f} h）\t"
                 + "\t".join(R4.format(Ce_d[di]) for di in d_idx))

    # ---------------- 交叉校核：分级网格 ----------------
    lines.append("")
    lines.append("网格与口径")
    lines.append("  交付口径：N=%d 均匀，Δr=%.6f cm，降采样至 0.1 cm；t_dry=%.4f s=%.4f h"
                 % (GRID_MAIN["N"], np.diff(r).min() * 100, t_dry, t_dry / 3600))
    lines.append("  中心最大值核查：%d/%d 个时刻的 max_r C 取在中心（判据可只用中心值）"
                 % (len(times) - n_off, len(times)))

    xc = run_case("校核", **GRID_XCHK)
    dt_x = (xc["t_dry"] - t_dry) / 60.0
    lines.append("  校核口径：N=%d 分级 q=%.1f，Δr∈[%.5f, %.4f] cm，节点数仅为交付口径的 1/%.0f；"
                 "t_dry=%.4f h" % (GRID_XCHK["N"], GRID_XCHK["q"],
                                   np.diff(xc["r"]).min() * 100,
                                   np.diff(xc["r"]).max() * 100,
                                   (GRID_MAIN["N"] + 1) / (GRID_XCHK["N"] + 1),
                                   xc["t_dry"] / 3600))
    lines.append("  两口径之差 %+.2f min（%.4f h），即 %.2f%%；结论不随网格改变。"
                 % (dt_x, dt_x / 60.0, abs(dt_x / 60.0) / (t_dry / 3600) * 100))

    # ---------------- 收敛性表 ----------------
    grid_txt = os.path.join(RES, "problem3_grid.txt")
    if os.path.exists(grid_txt):
        lines.append("")
        lines.append("网格收敛性（详见 problem3_grid.txt）")
        with open(grid_txt, encoding="utf-8") as f:
            for ln in f.read().splitlines():
                if ln.strip() and not ln.startswith("问题 3 网格收敛性研究"):
                    lines.append("  " + ln)

    # ---------------- 准稳态外推校核 ----------------
    w = res["solver"].w
    rows = qs_table(times, C, w, t_dry, C_inf_late, [30, 36, 42, 45, 48, 50, 52, 54])
    lines.append("")
    lines.append("准稳态（QS）外推校核（基准：二分精定位 t_dry = %.4f s = %.4f h）"
                 % (t_dry, t_dry / 3600))
    lines.append("起始/h\tC_中心\tφ̄\tφ_N\tQS平衡C_中心\t冻结形状外推/h\t误差/min"
                 "\t漂移修正外推/h\t误差/min")
    for q in rows:
        frz = "inf" if not np.isfinite(q["t_frz"]) else "%.4f" % (q["t_frz"] / 3600)
        efrz = "inf" if not np.isfinite(q["err_frz"]) else "%+.2f" % q["err_frz"]
        lines.append("%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%.4f\t%+.2f"
                     % (q["t0"] / 3600, q["C0"], q["phibar"], q["phiN"], q["C_eq"],
                        frz, efrz, q["t_drf"] / 3600, q["err_drf"]))

    # ---------------- 结论 ----------------
    lines.append("")
    lines.append("烘干时长：t = %.1f s = %.4f h = %.1f min" % (t_dry, t_dry / 3600, t_dry / 60))
    lines.append("按分钟向上取整（保守口径）：%.1f min = %.4f h（即 %.0f h %.0f min）"
                 % (t_min / 60, t_min / 3600,
                    np.floor(t_min / 3600), np.floor(t_min / 60) % 60))
    lines.append("终止时刻剖面（降采样到 0.1 cm 网格，kg/kg）：")
    for i in range(len(r_cm)):
        lines.append("  r=%.1f cm  C=%.6f" % (r_cm[i], Ce_d[i]))

    text = "\n".join(lines)
    print("\n" + text)
    out_txt = os.path.join(RES, "problem3_tables.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print("\n已写入", out_txt)

    # ---------------- result3.xlsx ----------------
    t_out = t_min
    n_row = _write_xlsx(times, Cd, r_cm, t_out)
    print("已写入 %s（%d 行 × %d 列）" % (os.path.join(RES, "result3.xlsx"), n_row, len(r_cm)))

    # ---------------- 缓存时空场供绘图复用 ----------------
    np.savez(os.path.join(RES, "p3_fields.npz"),
             times=times, r=r, T=T, C=C, t_dry=t_dry, t_min=t_min,
             r_out=p3.OUT_R, C_out=Cd, T_out=Td, C_end_ds=Ce_d,
             log_h=np.array([e["h"] for e in res["log"]]),
             log_t=np.array([e["t0"] for e in res["log"]]),
             log_m=np.array([e["m"] for e in res["log"]]),
             log_err=np.array([e["err"] for e in res["log"]]),
             bisect_hist=np.array(res["bisect_hist"]) if res["bisect_hist"] else np.zeros((0, 2)),
             qs_starts_h=np.array([q["t0"] / 3600.0 for q in rows]),
             qs_err_frz=np.array([q["err_frz"] for q in rows]),
             qs_err_drf=np.array([q["err_drf"] for q in rows]),
             xc_t_dry=xc["t_dry"], xc_r=xc["r"], xc_q=GRID_XCHK["q"])
    print("已写入", os.path.join(RES, "p3_fields.npz"))


def _write_xlsx(times, Cd, r_cm, t_out, path=None):
    """result3.xlsx：与附件 3 模板同构。

    单张 Sheet1，第 1 行表头「时间\\到药材中心的距离」+ 距离 0, 0.1, ..., 2（cm）；
    A 列为时间（s），自 60 起每 60 s 一行，直到 t_out（分钟向上取整后的烘干完成时刻）。
    """
    import openpyxl

    if path is None:
        path = os.path.join(RES, "result3.xlsx")
    wb = openpyxl.Workbook(write_only=True)
    ws = wb.create_sheet(title="Sheet1")
    ws.append(["时间\\到药材中心的距离"] + [round(float(v), 1) for v in r_cm])

    n_row = 0
    k_end = int(round(t_out / 60.0))
    for k in range(1, k_end + 1):
        t = 60.0 * k
        i = int(round(t / 60.0))
        if i >= len(times) or abs(times[i] - t) > 1e-6:
            break
        ws.append([t] + [round(float(v), 4) for v in Cd[i]])
        n_row += 1
    wb.save(path)
    return n_row


if __name__ == "__main__":
    main()
