# -*- coding: utf-8 -*-
"""
格式对比实验：
  A. 表面边界离散：'ghost'（虚拟节点）vs 'halfcell'（半控制体积守恒）
  B. D(C) 处理：单次滞后 vs Picard 迭代
评价指标：
  1) 与解析解的最大偏差（纯导热算例，验证边界离散精度）
  2) 空间收敛阶（以 N=640 为参考）
  3) 时间收敛阶（T 场与 C 场分别考察）
  4) 全局质量守恒相对偏差
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_model import (solve, RHO, CP, H_HEAT, H_MASS, R_CYL, T0, C0,
                      ALPHA, conservative_weights)
from p1_verify import analytic_cylinder, _radial_integral

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


Tf_bench = lambda t: 50.0
Cf_bench = lambda t: 0.05
Tf_var = lambda t: 28.0 + 22.0 * (1 - np.exp(-t / 400.0))
Cf_var = lambda t: 0.02 + 0.03 * (1 - np.exp(-t / 400.0))


# ---------------------------------------------------------------------------
def A_boundary():
    log("=" * 78)
    log("A. 表面 Robin 边界离散方式对比（纯导热解析解算例, T_inf=50degC, N=20）")
    log("=" * 78)
    log(f"  {'bc':>9} {'t(s)':>7} {'max|T_num-T_exact|':>20} {'表面误差':>12} {'中心误差':>12}")
    for bc in ("ghost", "halfcell"):
        for t in (300.0, 900.0, 1800.0):
            _, r, T_all, C_all = solve(Tf_bench, Cf_bench, t_end=t, dt=0.1, N=20,
                               theta=0.5, n_store=None, bc=bc)
            ex = analytic_cylinder(r, t, 50.0)
            e = np.abs(T_all[-1] - ex)
            log(f"  {bc:>9} {t:>7.0f} {e.max():>20.4e} {e[-1]:>12.4e} {e[0]:>12.4e}")
    log()

    log("  ---- 空间收敛阶（含传质，dt=0.02s，参考 N=640） ----")
    t_end = 600.0
    refs = {}
    for bc in ("ghost", "halfcell"):
        _, _, Tr_h, Cr_h = solve(Tf_var, Cf_var, t_end=t_end, dt=0.02, N=640,
                               theta=0.5, n_store=None, bc=bc)
        refs[bc] = (Tr_h[-1], Cr_h[-1])

    for bc in ("ghost", "halfcell"):
        Tr_ref, Cr_ref = refs[bc]
        log(f"  [{bc}]  {'N':>5} {'‖eT‖∞':>12} {'阶':>6} {'‖eC‖∞':>12} {'阶':>6}")
        prev = None
        for N in (20, 40, 80, 160):
            _, r, T_all, C_all = solve(Tf_var, Cf_var, t_end=t_end, dt=0.02, N=N,
                               theta=0.5, n_store=None, bc=bc)
            s = 640 // N
            eT = np.abs(T_all[-1] - Tr_ref[::s]).max()
            eC = np.abs(C_all[-1] - Cr_ref[::s]).max()
            oT = "" if prev is None else f"{np.log2(prev[0]/eT):.2f}"
            oC = "" if prev is None else f"{np.log2(prev[1]/eC):.2f}"
            log(f"  [{bc}]  {N:>5} {eT:>12.4e} {oT:>6} {eC:>12.4e} {oC:>6}")
            prev = (eT, eC)
    log()


# ---------------------------------------------------------------------------
def B_picard():
    log("=" * 78)
    log("B. D(C) 处理方式对比：单次滞后 vs Picard 迭代（bc=halfcell）")
    log("=" * 78)

    log("  ---- 时间收敛阶（N=40，参考 dt=0.02s，picard_max=30） ----")
    t_end = 600.0
    _, _, Tr_h, Cr_h = solve(Tf_var, Cf_var, t_end=t_end, dt=0.02, N=40, theta=0.5,
                             n_store=None, bc="halfcell", picard_max=30)
    Tr, Cr = Tr_h[-1], Cr_h[-1]
    for pm in (1, 30):
        log(f"  [picard_max={pm}]  {'dt':>6} {'‖eT‖∞':>11} {'阶':>6} {'‖eC‖∞':>11} {'阶':>6}")
        prev = None
        for dt in (4.0, 2.0, 1.0, 0.5):
            _, r, T_all, C_all = solve(Tf_var, Cf_var, t_end=t_end, dt=dt, N=40,
                               theta=0.5, n_store=None, bc="halfcell", picard_max=pm)
            eT = np.abs(T_all[-1] - Tr).max()
            eC = np.abs(C_all[-1] - Cr).max()
            oT = "" if prev is None else f"{np.log2(prev[0]/eT):.2f}"
            oC = "" if prev is None else f"{np.log2(prev[1]/eC):.2f}"
            log(f"  [picard_max={pm}]  {dt:>6.2f} {eT:>11.4e} {oT:>6} {eC:>11.4e} {oC:>6}")
            prev = (eT, eC)

    log()
    log("  ---- Picard 迭代次数对解的改善（dt=1s, N=20） ----")
    _, r, T_ref_h, C_ref_h = solve(Tf_var, Cf_var, t_end=600.0, dt=0.02, N=80,
                                   theta=0.5, n_store=None, bc="halfcell", picard_max=30)
    T_ref, C_ref = T_ref_h[-1], C_ref_h[-1]
    for pm in (1, 2, 3, 5, 10, 30):
        _, r, T_all, C_all = solve(Tf_var, Cf_var, t_end=600.0, dt=1.0, N=20,
                           theta=0.5, n_store=None, bc="halfcell", picard_max=pm)
        eT = np.abs(T_all[-1] - T_ref[::4]).max()
        eC = np.abs(C_all[-1] - C_ref[::4]).max()
        log(f"    picard_max={pm:>3}  ‖eC‖∞={eC:.4e}  ‖eT‖∞={eT:.4e}")
    log()


# ---------------------------------------------------------------------------
def C_conservation():
    log("=" * 78)
    log("C. 全局守恒性（dt=0.5s, N=20, t_end=1800s）")
    log("=" * 78)
    dt = 0.5
    for bc in ("ghost", "halfcell"):
        _, r, T_all, C_all = solve(lambda t: 28.0, Cf_var, t_end=1800.0, dt=dt, N=20,
                           theta=0.5, n_store=1, bc=bc)
        dr = R_CYL / 20
        w = conservative_weights(dr, 20)          # 格式自身的守恒权重

        # 用格式自身权重衡量的质量变化
        dM_fmt = RHO * 2 * np.pi * np.sum(w * (C_all[-1] - C_all[0]))
        # 用真实几何积分
        dM_geo = RHO * 2 * np.pi * (_radial_integral(C_all[-1], dr) - _radial_integral(C_all[0], dr))

        # 边界累计传质通量
        Mf = 0.0
        for n in range(1, len(T_all)):
            t0, t1 = (n - 1) * dt, n * dt
            Mf += (2 * np.pi * R_CYL * RHO * H_MASS * 0.5 *
                   ((C_all[n-1][-1] - Cf_var(t0)) + (C_all[n][-1] - Cf_var(t1))) * dt)

        log(f"  [{bc}]")
        log(f"    格式权重质量变化 ΔM_fmt = {dM_fmt: .6e} kg/m")
        log(f"    真实几何积分 ΔM_geo    = {dM_geo: .6e} kg/m")
        log(f"    边界累计传质   -M_flux = {-Mf: .6e} kg/m")
        log(f"    守恒偏差 |ΔM_fmt + M_flux|/|M_flux| = {abs(dM_fmt + Mf)/abs(Mf):.3e}")
        log(f"    几何积分偏差                          = {abs(dM_geo + Mf)/abs(Mf):.3e}")

    # 热平衡
    log()
    for bc in ("ghost", "halfcell"):
        Tf = lambda t: 28.0 + 22.0 * (1 - np.exp(-t / 400.0))
        _, r, T_all, C_all = solve(Tf, lambda t: 0.02, t_end=1800.0, dt=dt, N=20,
                           theta=0.5, n_store=1, bc=bc)
        dr = R_CYL / 20
        w = conservative_weights(dr, 20)
        dH_fmt = RHO * CP * 2 * np.pi * np.sum(w * (T_all[-1] - T_all[0]))
        Q = 0.0
        for n in range(1, len(T_all)):
            t0, t1 = (n - 1) * dt, n * dt
            Q += (2 * np.pi * R_CYL * H_HEAT *
                  (0.5 * ((Tf(t0) - T_all[n-1][-1]) + (Tf(t1) - T_all[n][-1]))) * dt)
        log(f"  [{bc}] 热平衡: ΔH={dH_fmt:.6e}  Q={Q:.6e}  偏差={abs(dH_fmt-Q)/abs(Q):.3e}")
    log()


if __name__ == "__main__":
    A_boundary()
    B_picard()
    C_conservation()
    with open(os.path.join(OUT, "problem1_format_comparison.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    print(f"\n已写入 {os.path.join(OUT, 'problem1_format_comparison.txt')}")
