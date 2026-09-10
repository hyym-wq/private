# -*- coding: utf-8 -*-
"""
问题 1 输入敏感性分析

S1  D(C) 非线性 vs 常系数 D≡D(C0)
    —— 量化"扩散系数随含水率变化"这一非线性项对结果的实际影响
S2  热物性 k / cp 的扰动
    —— 量化材料参数不确定性对温度场的影响

输出: results/problem1_sensitivity.txt
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p1_model as M
from p1_run import load_ambient, make_linear

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATT = os.path.join(os.path.dirname(BASE), "A题", "附件")
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


def ambient():
    t, T, C = load_ambient(os.path.join(ATT, "附件1.xlsx"))
    return make_linear(t, T), make_linear(t, C)


def run(Tf, Cf, N=20, dt=1.0, t_end=1800.0):
    _, r, Th, Ch = M.solve(Tf, Cf, t_end=t_end, dt=dt, N=N, theta=0.5,
                           n_store=1, bc="halfcell", picard_max=3)
    return r, Th[-1], Ch[-1]


def weighted_mean(u, N=20):
    w = np.asarray(M.conservative_weights(M.R_CYL / N, N), dtype=float)
    return float((w / w.sum()) @ np.asarray(u, dtype=float))


# ----------------------------------------------------------------------------
def s1_constant_D():
    log("=" * 74)
    log("S1  扩散系数的非线性: D(C) vs 常系数 D≡D(C0)")
    log("=" * 74)
    Tf, Cf = ambient()
    r, T_var, C_var = run(Tf, Cf)

    orig = M.diff_coef
    D0 = float(orig(M.C0))
    M.diff_coef = lambda C: np.full_like(np.asarray(C, dtype=float), D0)
    try:
        _, T_con, C_con = run(Tf, Cf)
    finally:
        M.diff_coef = orig          # 必须还原，否则污染后续算例

    log(f"  D(C0) = {D0:.4e} m^2/s")
    log(f"  {'位置':<14}{'D(C) 变系数':>14}{'D≡const':>12}{'相对差':>10}")
    for j, lab in ((0, "中心 r=0"), (10, "中间 r=1.0cm"), (20, "表面 r=2.0cm")):
        a, b = C_var[j], C_con[j]
        log(f"  {lab:<14}{a:>14.4f}{b:>12.4f}{(b - a) / a * 100:>9.2f}%")

    mv, mc = weighted_mean(C_var), weighted_mean(C_con)
    log("")
    log(f"  体积平均 C : 变系数 {mv:.4f}   常系数 {mc:.4f}")
    log(f"  整体失水量 : 变系数 {M.C0 - mv:.4f}   常系数 {M.C0 - mc:.4f} kg/kg")
    log(f"  常系数使整体失水量高估 {((M.C0 - mc) - (M.C0 - mv)) / (M.C0 - mv) * 100:.2f}%")
    log("")
    log("  结论：常系数近似高估失水量，方向与 D 随 C 单调下降(负反馈)一致；")
    log("        但 30 min 内幅度仅 ~1%，因浓度变化集中于表层，主体区域仍接近 C0。")
    log()
    return C_var, C_con


# ----------------------------------------------------------------------------
def s2_thermal_props():
    log("=" * 74)
    log("S2  热物性参数扰动: k / (rho*cp) 变化对温度场的影响")
    log("=" * 74)
    Tf, Cf = ambient()
    alpha0 = M.ALPHA
    log(f"  基准 alpha = {alpha0:.6e} m^2/s  (k={M.K_COND}, rho={M.RHO}, cp={M.CP})")
    # 先算基准，再算扰动，才能给出相对基准的偏差
    # 注意 run() 返回 (r, T, C)，解包顺序不能取错
    M.ALPHA = alpha0
    _, T_base, _ = run(Tf, Cf)

    log(f"  {'alpha 倍率':>11}{'alpha':>16}{'中心T':>10}{'表面T':>10}"
        f"{'Δ中心T':>10}{'Δ表面T':>10}")
    for f in (0.8, 0.9, 1.0, 1.1, 1.2):
        M.ALPHA = alpha0 * f
        try:
            _, T, _ = run(Tf, Cf)
        finally:
            M.ALPHA = alpha0          # 必须还原，否则污染后续算例
        log(f"  {f:>10.1f}x{alpha0 * f:>15.5e}{T[0]:>10.4f}{T[-1]:>10.4f}"
            f"{T[0] - T_base[0]:>+10.4f}{T[-1] - T_base[-1]:>+10.4f}")
    log("")
    log("  结论：alpha 变化 ±20% 对 1800 s 时的温度场影响有限——")
    log("        中心温度约 ∓0.6 degC，表面温度仅约 ±0.2 degC。")
    log("        表面温度对 alpha 不敏感且符号相反：alpha 越小，表面热量越难向内传导，")
    log("        表面温度反而越接近环境温度。故材料热物性的不确定性主要影响内部温度。")
    log("        又 alpha ∝ k/(rho*cp)，k、rho、cp 的相对误差以相同权重传递到 alpha。")
    log()


if __name__ == "__main__":
    s1_constant_D()
    s2_thermal_props()
    with open(os.path.join(OUT, "problem1_sensitivity.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    print(f"\n敏感性报告已写入 {os.path.join(OUT, 'problem1_sensitivity.txt')}")
