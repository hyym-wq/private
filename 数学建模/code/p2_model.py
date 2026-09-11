# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 2
非线性变参数热质传递模型（线性化隐式 FDM + 迭代）

控制方程（径向一维，无限长圆柱 R=2cm，长 L=25cm）:
    热:  rho(C) cp(C) ∂T/∂t = (1/r) ∂/∂r ( k(C) r ∂T/∂r )
    质:  ∂C/∂t            = (1/r) ∂/∂r ( D(C,T) r ∂C/∂r )

变参数经验公式（附录 3，问题 2/3 统一采用）:
    rho(C) = 650 + 128 C                                    [kg/m^3]
    cp(C)  = 1450 + 2736 C/(C+1)                            [J/(kg K)]
    k(C)   = 0.21 + 0.38 C/(C+1)                            [W/(m K)]
    D(C,T) = 2.4e-3 exp(-0.45/C) exp(-3850/T_K)             [m^2/s], T_K = T + 273.15

定解条件:
    初值: T(r,0)=28 degC,  C(r,0)=2.55 kg/kg
    r=0 : ∂T/∂r = ∂C/∂r = 0                                （轴对称）
    r=R : -k ∂T/∂r|_R = h  (T_R - T_inf)                   （对流换热 Robin）
          -D ∂C/∂r|_R = hm (C_R - C_inf)                   （对流传质 Robin）
    环境: 预热平衡阶段按附件 1（分段线性插值），恒温干燥阶段维持末值。

数值格式（创新处理）:
    * 隐式 θ 格式（θ=1/2 Crank--Nicolson，无条件稳定）
    * D(C,T) 在当前状态处 Taylor 展开一阶线性化（Newton 型线性化，含 ∂D/∂C 敏度项）
    * 每个时间步内 Picard 迭代，2--3 次收敛
    * 温度方程与水分方程弱耦合交替求解
    * Thomas 追赶法求解三对角方程组
"""

import os
import sys
import numpy as np
import openpyxl

# ----------------------------------------------------------------------------
# 物理参数（附录 2 中的 h、h_m 沿用；材料参数用附录 3）
# ----------------------------------------------------------------------------
R_CYL = 0.02         # m   圆柱半径 2 cm
L_CYL = 0.25         # m   圆柱长度 25 cm
H_HEAT = 25.0        # W/(m^2 K)  对流换热系数
H_MASS = 8.0e-7      # m/s        对流传质系数

T0 = 28.0            # degC  初始温度
C0 = 2.55            # kg/kg 初始水分浓度（干基含水率）

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATT = os.path.join(os.path.dirname(BASE), "A题", "附件")


# ----------------------------------------------------------------------------
# 附录 3 变参数经验公式
# ----------------------------------------------------------------------------
def rho(C):
    C = np.asarray(C, dtype=float)
    return 650.0 + 128.0 * C


def cp(C):
    C = np.asarray(C, dtype=float)
    return 1450.0 + 2736.0 * C / (C + 1.0)


def kcond(C):
    C = np.asarray(C, dtype=float)
    return 0.21 + 0.38 * C / (C + 1.0)


def diff_coef(C, T):
    """D(C,T) = 2.4e-3 exp(-0.45/C) exp(-3850/T_K)，T 单位 °C。"""
    C = np.asarray(C, dtype=float)
    T = np.asarray(T, dtype=float)
    Tk = T + 273.15
    return 2.4e-3 * np.exp(-0.45 / np.clip(C, 1e-9, None)) * np.exp(-3850.0 / Tk)


def dD_dC(C, T):
    """∂D/∂C = D · 0.45/C^2（Taylor 一阶线性化所需的敏度系数）。"""
    C = np.asarray(C, dtype=float)
    return diff_coef(C, T) * 0.45 / np.clip(C, 1e-9, None) ** 2


def dD_dT(C, T):
    """∂D/∂T = D · 3850/T_K^2（对摄氏温标的导数 = 对开尔文的导数）。"""
    T = np.asarray(T, dtype=float)
    Tk = T + 273.15
    return diff_coef(C, T) * 3850.0 / Tk ** 2


# ----------------------------------------------------------------------------
# Thomas 追赶法
# ----------------------------------------------------------------------------
def thomas(a, b, c, d):
    """求解三对角方程组 M x = d（a 次对角、b 主对角、c 超对角）。"""
    n = b.size
    psi = np.empty(n)
    phi = np.empty(n)
    psi[0] = c[0] / b[0]
    phi[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * psi[i - 1]
        psi[i] = c[i] / m
        phi[i] = (d[i] - a[i] * phi[i - 1]) / m
    x = np.empty(n)
    x[-1] = phi[-1]
    for i in range(n - 2, -1, -1):
        x[i] = phi[i] - psi[i] * x[i + 1]
    return x


def _shift(u, k):
    v = np.zeros_like(u)
    if k == -1:
        v[1:] = u[:-1]
    else:
        v[:-1] = u[1:]
    return v


def _step(u, op_new, op_old, s_old, s_new, dt, theta):
    """线性 θ 格式单步（允许新旧时间层算子不同）。"""
    an, bn, cn = op_new
    ao, bo, co = op_old
    rhs = u + (1.0 - theta) * dt * (ao * _shift(u, -1) + bo * u + co * _shift(u, +1))
    rhs += dt * (theta * s_new + (1.0 - theta) * s_old)
    lo = -theta * dt * an
    di = 1.0 - theta * dt * bn
    up = -theta * dt * cn
    lo[0] = 0.0
    up[-1] = 0.0
    return thomas(lo, di, up, rhs)


# ----------------------------------------------------------------------------
# 热传导算子（线性于 T，系数冻结于 C_f）
# ----------------------------------------------------------------------------
def build_op_heat(dr, N, C_f, R=R_CYL):
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    inv = 1.0 / (rho(C_f) * cp(C_f))          # 1/(rho cp)，逐节点
    k = kcond(C_f)
    k_half = 0.5 * (k[:-1] + k[1:])           # 界面导热系数 k_{j+1/2}

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    # 中心 r=0：L'Hopital，dT/dt = 4 [k_{1/2}/(rho cp)_0]/dr^2 (T1-T0)
    b[0] = -4.0 * k_half[0] * inv[0] / dr ** 2
    c[0] = 4.0 * k_half[0] * inv[0] / dr ** 2

    # 内部点：守恒型有限体积
    j = np.arange(1, N)
    A = k_half[j - 1] * r_half[j - 1] * inv[j] / (r[j] * dr ** 2)
    Cc = k_half[j] * r_half[j] * inv[j] / (r[j] * dr ** 2)
    a[1:N] = A; c[1:N] = Cc; b[1:N] = -(A + Cc)

    # 边界 r=R：半控制体积精确平衡
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    a[N] = k_half[N - 1] * r_half[N - 1] * inv[N] / (VN * dr)
    gT = R * H_HEAT * inv[N] / VN
    b[N] = -a[N] - gT
    return a, b, c, gT


# ----------------------------------------------------------------------------
# 传质算子：扩散部分 L(C) 与 Newton 雅可比 J*
#   将边界 Robin 项拆为 L(C)_N = [扩散] - gC·C_N 与源项 s_N = gC·C_inf。
# ----------------------------------------------------------------------------
def mass_diffusion(dr, N, C, T, R=R_CYL):
    """L(C)：扩散算子作用于 C（含边界对流流出 -gC·C_N，不含 C_inf 源）。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    D = diff_coef(C, T)
    D_half = 0.5 * (D[:-1] + D[1:])
    out = np.zeros(N + 1)
    out[0] = 4.0 * D_half[0] * (C[1] - C[0]) / dr ** 2
    j = np.arange(1, N)
    out[j] = (D_half[j - 1] * r_half[j - 1] * (C[j - 1] - C[j])
              + D_half[j] * r_half[j] * (C[j + 1] - C[j])) / (r[j] * dr ** 2)
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    gC = R * H_MASS / VN
    out[N] = D_half[N - 1] * r_half[N - 1] * (C[N - 1] - C[N]) / (VN * dr) - gC * C[N]
    return out


def mass_jacobian(dr, N, C_star, T_star, R=R_CYL):
    """Newton 雅可比 J* = ∂L/∂C|_{C*}（Taylor 一阶线性化，含 ∂D/∂C 敏度项）。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    D = diff_coef(C_star, T_star)
    Dp = dD_dC(C_star, T_star)
    D_half = 0.5 * (D[:-1] + D[1:])

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    # 中心 j=0
    d0 = C_star[1] - C_star[0]
    b[0] = 4.0 / dr ** 2 * (-D_half[0] + 0.5 * Dp[0] * d0)
    c[0] = 4.0 / dr ** 2 * (D_half[0] + 0.5 * Dp[1] * d0)

    # 内部点
    j = np.arange(1, N)
    dm = C_star[j - 1] - C_star[j]     # C*_{j-1} - C*_j
    dp = C_star[j + 1] - C_star[j]     # C*_{j+1} - C*_j
    a[1:N] = r_half[j - 1] / (r[j] * dr ** 2) * (D_half[j - 1] + 0.5 * Dp[j - 1] * dm)
    c[1:N] = r_half[j] / (r[j] * dr ** 2) * (D_half[j] + 0.5 * Dp[j + 1] * dp)
    b[1:N] = (-D_half[j - 1] * r_half[j - 1] - D_half[j] * r_half[j]
              + 0.5 * Dp[j] * (r_half[j - 1] * dm + r_half[j] * dp)) / (r[j] * dr ** 2)

    # 边界 j=N
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    gC = R * H_MASS / VN
    dN = C_star[N - 1] - C_star[N]
    a[N] = r_half[N - 1] / (VN * dr) * (D_half[N - 1] + 0.5 * Dp[N - 1] * dN)
    b[N] = r_half[N - 1] / (VN * dr) * (-D_half[N - 1] + 0.5 * Dp[N] * dN) - gC
    return a, b, c


def build_op_mass_lag(dr, N, C_f, T_f, R=R_CYL):
    """冻结系数（0 阶滞后）传质算子，用于对照 Taylor 线性化的收益。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    D = diff_coef(C_f, T_f)
    D_half = 0.5 * (D[:-1] + D[1:])
    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)
    b[0] = -4.0 * D_half[0] / dr ** 2; c[0] = 4.0 * D_half[0] / dr ** 2
    j = np.arange(1, N)
    A = D_half[j - 1] * r_half[j - 1] / (r[j] * dr ** 2)
    Cc = D_half[j] * r_half[j] / (r[j] * dr ** 2)
    a[1:N] = A; c[1:N] = Cc; b[1:N] = -(A + Cc)
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    a[N] = D_half[N - 1] * r_half[N - 1] / (VN * dr)
    gC = R * H_MASS / VN
    b[N] = -a[N] - gC
    return a, b, c, gC


# ----------------------------------------------------------------------------
# 环境条件（附件 1）
# ----------------------------------------------------------------------------
def load_ambient(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Sheet1"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    t = np.array([float(r[0]) for r in rows])
    T = np.array([float(r[1]) for r in rows])
    C = np.array([float(r[2]) for r in rows])
    return t, T, C


def make_ambient(path=None):
    """返回 (T_inf_fn, C_inf_fn)。附件 1 覆盖 0--14400 s，其后维持末值（恒温干燥）。"""
    if path is None:
        path = os.path.join(ATT, "附件1.xlsx")
    t, T, C = load_ambient(path)

    def T_inf(tv):
        return float(np.interp(tv, t, T))

    def C_inf(tv):
        return float(np.interp(tv, t, C))

    return T_inf, C_inf


# ----------------------------------------------------------------------------
# 传质单步（Newton / 滞后两种线性化）
# ----------------------------------------------------------------------------
def _mass_step(C_old, C_star, T_old, T_new, Ci_o, Ci_n, dt, theta, dr, N, R, newton):
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    gC = R * H_MASS / VN
    s_n = np.zeros(N + 1); s_n[N] = gC * Ci_n
    s_o = np.zeros(N + 1); s_o[N] = gC * Ci_o

    if newton:
        a, b, c = mass_jacobian(dr, N, C_star, T_new, R)
        Lstar = mass_diffusion(dr, N, C_star, T_new, R)
        Lold = mass_diffusion(dr, N, C_old, T_old, R)
        JCs = np.zeros(N + 1)
        JCs[0] = b[0] * C_star[0] + c[0] * C_star[1]
        JCs[1:N] = a[1:N] * C_star[0:N - 1] + b[1:N] * C_star[1:N] + c[1:N] * C_star[2:N + 1]
        JCs[N] = a[N] * C_star[N - 1] + b[N] * C_star[N]
        rhs = C_old + dt * (theta * (Lstar - JCs + s_n) + (1.0 - theta) * (Lold + s_o))
        lo = -theta * dt * a; lo[0] = 0.0
        di = 1.0 - theta * dt * b
        up = -theta * dt * c; up[-1] = 0.0
        return thomas(lo, di, up, rhs)
    else:
        am_n, bm_n, cm_n, _ = build_op_mass_lag(dr, N, C_star, T_new, R)
        am_o, bm_o, cm_o, _ = build_op_mass_lag(dr, N, C_old, T_old, R)
        return _step(C_old, (am_n, bm_n, cm_n), (am_o, bm_o, cm_o), s_o, s_n, dt, theta)


# ----------------------------------------------------------------------------
# 主求解器
# ----------------------------------------------------------------------------
def solve(T_inf_fn, C_inf_fn, t_end=259200.0, dt=1.0, N=20, R=R_CYL, theta=0.5,
          n_store=None, picard_max=3, picard_tol=1e-9, newton=True,
          T_init=None, C_init=None, iters_out=None):
    dr = R / N
    r = np.arange(N + 1) * dr
    nsteps = int(round(t_end / dt))

    T = np.full(N + 1, T0 if T_init is None else float(T_init))
    C = np.full(N + 1, C0 if C_init is None else float(C_init))

    times, T_hist, C_hist = [], [], []
    if n_store:
        times.append(0.0); T_hist.append(T.copy()); C_hist.append(C.copy())

    zero = np.zeros(N + 1)
    iter_sum = 0.0

    for n in range(1, nsteps + 1):
        t_new, t_old = n * dt, (n - 1) * dt
        Ti_n, Ti_o = T_inf_fn(t_new), T_inf_fn(t_old)
        Ci_n, Ci_o = C_inf_fn(t_new), C_inf_fn(t_old)

        T_old, C_old = T.copy(), C.copy()
        C_star, T_star = C.copy(), T.copy()

        # 旧时间层热算子（系数取 C^n，只需组装一次）
        ah_o, bh_o, ch_o, gT_o = build_op_heat(dr, N, C_old, R)
        sT_o = zero.copy(); sT_o[N] = gT_o * Ti_o

        iters = picard_max
        for k in range(picard_max):
            # ① 热：冻结 C*，线性求解 T（Thomas）
            ah_n, bh_n, ch_n, gT_n = build_op_heat(dr, N, C_star, R)
            sT_n = zero.copy(); sT_n[N] = gT_n * Ti_n
            T_new = _step(T_old, (ah_n, bh_n, ch_n), (ah_o, bh_o, ch_o), sT_o, sT_n, dt, theta)

            # ② 质：冻结 T_new，D(C,T) 于 C* 处线性化后求解 C
            C_new = _mass_step(C_old, C_star, T_old, T_new, Ci_o, Ci_n, dt, theta, dr, N, R, newton)

            dT = float(np.abs(T_new - T_star).max())
            dC = float(np.abs(C_new - C_star).max())
            T_star, C_star = T_new, C_new
            if dT < picard_tol and dC < picard_tol:
                iters = k + 1
                break
        iter_sum += iters
        T, C = T_new, C_new

        if n_store and (n % n_store == 0):
            times.append(t_new); T_hist.append(T.copy()); C_hist.append(C.copy())

    if iters_out is not None:
        iters_out.append(iter_sum / nsteps)

    if not n_store:
        times.append(t_end); T_hist.append(T.copy()); C_hist.append(C.copy())
    return np.array(times), r, np.array(T_hist), np.array(C_hist)


# ----------------------------------------------------------------------------
# 守恒权重（与问题 1 相同，Σ w_j = R^2/2）
# ----------------------------------------------------------------------------
def conservative_weights(dr, N, R=R_CYL):
    w = np.arange(N + 1) * dr * dr
    w[0] = dr ** 2 / 8.0
    w[N] = R * dr / 2.0 - dr ** 2 / 8.0
    return w


def mass_conservation_residual(dr, N, C, C_inf, R=R_CYL):
    """全局质量平衡回代残差（应处于浮点舍入水平）。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    VN = R * dr / 2.0 - dr ** 2 / 8.0
    gC = R * H_MASS / VN
    w = conservative_weights(dr, N, R)
    L = mass_diffusion(dr, N, C, C, R)  # 温度不参与，仅用于结构；D 用当前 C
    # 源项：表面仅存
    s = np.zeros(N + 1); s[N] = gC * C_inf
    total = np.sum(w * (L + s))
    return total


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    T_inf, C_inf = make_ambient()
    print("环境: T_inf(0)=%.3f  T_inf(14400)=%.3f  C_inf(0)=%.5f  C_inf(14400)=%.5f"
          % (T_inf(0), T_inf(14400), C_inf(0), C_inf(14400)))

    # 冒烟测试：3 小时（10800 s）
    iters = []
    times, r, T, C = solve(T_inf, C_inf, t_end=10800.0, dt=1.0, N=20, n_store=None,
                           picard_max=3, newton=True, iters_out=iters)
    print("\n[t_end=3h, N=20, dt=1s, Newton 线性化]")
    print("平均 Picard 迭代次数 = %.2f" % iters[0])
    print("中心:  T=%.6f  C=%.6f" % (T[-1][0], C[-1][0]))
    print("表面:  T=%.6f  C=%.6f" % (T[-1][-1], C[-1][-1]))
    print("T 范围 [%.4f, %.4f], C 范围 [%.4f, %.4f]" % (T.min(), T.max(), C.min(), C.max()))
