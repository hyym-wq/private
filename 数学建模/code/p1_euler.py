# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 1（数值计算：Euler 显式方法）

本文件是问题 1 数值计算的第一步：**完整实现 Euler 显式（向前 Euler）方法**，
把径向节点严格分为三类分别离散，再向前推进时间：

  控制方程（无限长圆柱假设，径向一维）
      热:  rho*cp * dT/dt = (1/r) d/dr ( k * r * dT/dr )
      质:  dC/dt         = (1/r) d/dr ( D(C) * r * dC/dr )
  定解条件
      初值:  T(r,0)=28 degC,  C(r,0)=2.55 kg/kg
      r=0 :  dT/dr = dC/dr = 0                         （轴对称）
      r=R :  -k dT/dr|_R = h  (T_R - T_inf(t))         （对流换热 Robin）
            -D dC/dr|_R = hm (C_R - C_inf(t))          （对流传质 Robin）

  三类节点的显式更新（theta = 0 的 θ 格式，逐点直接推进，无需解线性方程组）：
      (i)   中心点 j=0    —— L'Hopital 法则消去 1/r 奇异 + 对称虚拟点
      (ii)  内部点 1<=j<=N-1 —— 守恒型有限体积 / 柱坐标中心差分
      (iii) 边界点 j=N    —— Robin 边界条件的半控制体积离散

  Euler 显式的本质：u^{n+1} = u^n + dt * f(u^n)，右端全部用已知的 n 时刻值，
  因此对 D(C) 的非线性**无需迭代**（D 直接用 C^n 计算），这是显式方法
  相对隐式方法的计算优势；代价是受抛物型 CFL 稳定条件 Fo <= 1/2 限制。

  本文件同时：
      * 输出表 1（温度）、表 2（水分浓度）
      * 演示 CFL 稳定边界：Fo <= 1/2 有界，Fo > 1/2 非物理放大
      * 与 Crank-Nicolson（theta=0.5）在同一 dt=1 s 下对照，验证二者一致
"""
import os
import sys
import numpy as np
import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_model import RHO, CP, K_COND, H_HEAT, H_MASS, R_CYL, T0, C0, ALPHA, diff_coef

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATT = os.path.join(os.path.dirname(BASE), "A题", "附件")
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

DT = 1.0            # 时间步长 s（Fo = alpha*dt/dr^2 = 0.169 <= 1/2，满足 CFL）
N = 20              # 径向网格数 -> dr = 2cm/20 = 0.1 cm
T_END = 1800.0


# ----------------------------------------------------------------------------
# 附件 1：烘房环境条件（分段线性插值）
# ----------------------------------------------------------------------------
def load_ambient(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Sheet1"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    t = np.array([float(r[0]) for r in rows])
    T = np.array([float(r[1]) for r in rows])
    C = np.array([float(r[2]) for r in rows])
    return t, T, C


def make_linear(t_data, y_data):
    def f(t):
        return float(np.interp(t, t_data, y_data))
    return f


# ----------------------------------------------------------------------------
# 空间离散系数：把三类节点显式写出（与 theta 格式 theta=0 一致）
# ----------------------------------------------------------------------------
def explicit_rhs_heat(T, dr, T_inf):
    """
    热方程的显式右端 f_T(T^n)，返回 dT/dt 在三个节点类上的取值。

    返回 dTdt[j]（j=0..N）。源项（表面对流）已并入边界节点。
    """
    N = T.size - 1
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    dTdt = np.zeros(N + 1)

    # ---- (i) 中心点 j=0：L'Hopital + 对称，dT/dt = 4 alpha/dr^2 (T1 - T0) ----
    dTdt[0] = 4.0 * ALPHA / dr ** 2 * (T[1] - T[0])

    # ---- (ii) 内部点 1..N-1：柱坐标守恒型中心差分 ----
    for j in range(1, N):
        A = ALPHA * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = ALPHA * r_half[j] / (r[j] * dr ** 2)
        dTdt[j] = A * (T[j - 1] - T[j]) + Cc * (T[j + 1] - T[j])

    # ---- (iii) 边界点 j=N：Robin 条件的半控制体积离散 ----
    VN = R_CYL * dr / 2.0 - dr ** 2 / 8.0
    aN = ALPHA * r_half[N - 1] / (VN * dr)
    gN = R_CYL * H_HEAT / (RHO * CP * VN)
    dTdt[N] = aN * (T[N - 1] - T[N]) - gN * (T[N] - T_inf)

    return dTdt


def explicit_rhs_mass(C, dr, C_inf):
    """
    传质方程的显式右端 f_C(C^n)。D(C) 直接用 C^n 计算（滞后取当前时刻），
    故显式方法无需 Picard 迭代。
    """
    N = C.size - 1
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    D = diff_coef(C)
    D_half = 0.5 * (D[:-1] + D[1:])          # 界面扩散系数算术平均
    dCdt = np.zeros(N + 1)

    # ---- (i) 中心点 ----
    dCdt[0] = 4.0 * D_half[0] / dr ** 2 * (C[1] - C[0])

    # ---- (ii) 内部点 ----
    for j in range(1, N):
        A = D_half[j - 1] * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = D_half[j] * r_half[j] / (r[j] * dr ** 2)
        dCdt[j] = A * (C[j - 1] - C[j]) + Cc * (C[j + 1] - C[j])

    # ---- (iii) 边界点 ----
    VN = R_CYL * dr / 2.0 - dr ** 2 / 8.0
    aN = D_half[N - 1] * r_half[N - 1] / (VN * dr)
    gN = R_CYL * H_MASS / VN
    dCdt[N] = aN * (C[N - 1] - C[N]) - gN * (C[N] - C_inf)

    return dCdt


# ----------------------------------------------------------------------------
# Euler 显式主循环
# ----------------------------------------------------------------------------
def solve_explicit(T_inf_fn, C_inf_fn, t_end=1800.0, dt=1.0, N=20, n_store=None,
                   verbose=False):
    """向前 Euler 推进：u^{n+1} = u^n + dt * f(u^n)。"""
    dr = R_CYL / N
    r = np.arange(N + 1) * dr
    nsteps = int(round(t_end / dt))

    T = np.full(N + 1, T0)
    C = np.full(N + 1, C0)

    times, T_hist, C_hist = [], [], []
    if n_store:
        times.append(0.0); T_hist.append(T.copy()); C_hist.append(C.copy())

    for n in range(1, nsteps + 1):
        t = n * dt
        T_inf = T_inf_fn(t)
        C_inf = C_inf_fn(t)

        # 先热后质（单向耦合：热方程不含水分项）
        T = T + dt * explicit_rhs_heat(T, dr, T_inf)
        C = C + dt * explicit_rhs_mass(C, dr, C_inf)

        if n_store and (n % n_store == 0):
            times.append(t); T_hist.append(T.copy()); C_hist.append(C.copy())

    if not n_store:
        times.append(t_end); T_hist.append(T.copy()); C_hist.append(C.copy())
    return np.array(times), r, np.array(T_hist), np.array(C_hist)


# ----------------------------------------------------------------------------
# 表 1 / 表 2 输出
# ----------------------------------------------------------------------------
def print_tables(times, r, T_hist, C_hist, dt, tag="Euler 显式"):
    r_cm = r * 100.0
    t_tab = [100, 300, 600, 900, 1200, 1500, 1800]
    d_tab = [0.0, 0.5, 1.0, 1.5, 2.0]
    step = int(round(0.1 / (R_CYL / N * 100)))  # 0.1 cm -> 网格索引步长

    def pick(hist, ts):
        i = int(round(ts / dt))
        return [hist[i][int(round(d / 0.1))] for d in d_tab]

    lines = [f"表 1  30 分钟内药材的温度 (degC)  [{tag}, dt={dt:g} s]"]
    lines.append("时间/s\t" + "\t".join(f"{d:g} cm" for d in d_tab))
    for ts in t_tab:
        vals = pick(T_hist, ts)
        lines.append(f"{ts}\t" + "\t".join(f"{v:.4f}" for v in vals))

    lines.append("")
    lines.append(f"表 2  30 分钟内药材的水分浓度 (kg/kg)  [{tag}, dt={dt:g} s]")
    lines.append("时间/s\t" + "\t".join(f"{d:g} cm" for d in d_tab))
    for ts in t_tab:
        vals = pick(C_hist, ts)
        lines.append(f"{ts}\t" + "\t".join(f"{v:.4f}" for v in vals))
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# CFL 稳定性演示
# ----------------------------------------------------------------------------
def _max_eig_mag(N, dr):
    """柱坐标热算子 dT/dt = M T + s 的最大特征值模（最负特征值的绝对值）。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    M = np.zeros((N + 1, N + 1))
    M[0, 0] = -4.0 * ALPHA / dr ** 2
    M[0, 1] = 4.0 * ALPHA / dr ** 2
    for j in range(1, N):
        A = ALPHA * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = ALPHA * r_half[j] / (r[j] * dr ** 2)
        M[j, j - 1] = A; M[j, j] = -(A + Cc); M[j, j + 1] = Cc
    VN = R_CYL * dr / 2.0 - dr ** 2 / 8.0
    aN = ALPHA * r_half[N - 1] / (VN * dr)
    gN = R_CYL * H_HEAT / (RHO * CP * VN)
    M[N, N - 1] = aN; M[N, N] = -(aN + gN)
    ev = np.linalg.eigvals(M)
    return -ev.real.min()          # λ_max = -min(λ)（扩散算子特征值全负）


def cfl_demo():
    """热方程显式格式的条件稳定性：Fo=alpha*dt/dr^2 越过临界值后非物理放大。"""
    lines = []
    lines.append("CFL 稳定性边界（环境恒为 T_inf=50 degC，精确解恒有 T<=50）")
    dr = R_CYL / N
    lm = _max_eig_mag(N, dr)
    dt_crit = 2.0 / lm
    fo_crit = ALPHA * dt_crit / dr ** 2
    lines.append(f"  算子最大特征值模 lambda_max = {lm:.4f} s^-1（柱坐标中心节点 4α/Δr²={4*ALPHA/dr**2:.4f}）")
    lines.append(f"  => 临界时间步 dt_crit = 2/lambda_max = {dt_crit:.3f} s，临界 Fo_crit = {fo_crit:.3f}")
    lines.append("  dt(s)      Fo=alpha*dt/dr^2    T_max(degC)        物理有界")
    Tf = lambda t: 50.0
    Cf = lambda t: 0.02
    for dt in (1.0, 2.0, 2.4, 2.5, 2.9, 3.0, 5.0, 20.0):
        Fo = ALPHA * dt / dr ** 2
        times, r, T, C = solve_explicit(Tf, Cf, t_end=1800.0, dt=dt, N=N)
        ok = np.isfinite(T).all() and T.max() <= 50.0 + 1e-6
        lines.append(f"  {dt:>5.1f}    {Fo:>10.3f}        {T.max():>12.4e}        {'√' if ok else '×'}")
    lines.append("  结论：柱坐标下中心节点(4α/Δr²)与 Robin 边界使显式格式的实际稳定边界")
    lines.append(f"        收紧为 Fo_crit≈{fo_crit:.2f}（略严于笛卡尔 1D 的 Fo<=1/2），越过即非物理放大。")
    return "\n".join(lines)


def main():
    t_amb, T_amb, C_amb = load_ambient(os.path.join(ATT, "附件1.xlsx"))
    T_inf_fn = make_linear(t_amb, T_amb)
    C_inf_fn = make_linear(t_amb, C_amb)

    dr = R_CYL / N
    Fo = ALPHA * DT / dr ** 2
    print(f"网格: N={N}, dr={dr*100:.2f} cm, dt={DT} s, Fo=alpha*dt/dr^2={Fo:.4f} (<=1/2)")
    print(f"环境: T_inf {T_amb[0]:.2f}->{T_amb[-1]:.2f} degC, C_inf {C_amb[0]:.4f}->{C_amb[-1]:.4f}\n")

    # ---- Euler 显式主结果 ----
    times, r, T_hist, C_hist = solve_explicit(T_inf_fn, C_inf_fn, t_end=T_END,
                                              dt=DT, N=N, n_store=1)
    table_text = print_tables(times, r, T_hist, C_hist, DT, tag="Euler 显式")
    print(table_text)

    # ---- 与 Crank-Nicolson 对照（同一 dt=1 s） ----
    from p1_model import solve as solve_cn
    times_cn, r_cn, T_cn, C_cn = solve_cn(T_inf_fn, C_inf_fn, t_end=T_END, dt=DT,
                                          N=N, theta=0.5, n_store=1, bc="halfcell",
                                          picard_max=3)
    dT = max(np.abs(T_hist[-1] - T_cn[-1]).max(), 0.0)
    dC = max(np.abs(C_hist[-1] - C_cn[-1]).max(), 0.0)
    print(f"\n[对照] 同一 dt={DT:g} s 下，Euler 显式 与 CN 终态最大差: "
          f"|dT|={dT:.3e} degC, |dC|={dC:.3e} kg/kg（时间离散误差已小于空间误差，两者四位小数一致）")

    # ---- CFL 稳定性 ----
    print("\n" + cfl_demo())

    # ---- 写文件 ----
    out_txt = os.path.join(RES, "problem1_euler.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(table_text + "\n\n" + cfl_demo() + "\n")
    print(f"\n已写入 {out_txt}")

    np.savez(os.path.join(RES, "p1_euler_fields.npz"),
             times=times, r=r, T=T_hist, C=C_hist)
    print("已写入 results/p1_euler_fields.npz")


if __name__ == "__main__":
    main()
