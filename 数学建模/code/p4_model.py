# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 4
坐标变换法（Landau 变换）求解移动边界问题

物理背景：烘干过程中药材因水分流失发生尺寸变化，半径 R(t) 由附件 2 给出
（2 cm → 1.198 cm，随后维持不变）。此时扩散区域 [0, R(t)] 是随时间移动的，
不能用问题 1--3 的固定网格直接求解。

核心思想（Landau 变换）：
    引入归一化坐标  E = r / R(t)，把移动边界 [0, R(t)] 转化为固定域 [0, 1]。
    网格随边界等比例缩放（E 网格固定、物理网格 r_j = E_j · R(t)），
    因而空间分辨率（相对分辨率）在全程保持一致。

变换后的控制方程（径向一维、无限长圆柱、附录 4 变参数物性）：

    热:  ∂T/∂τ = (1/(ρc_p R²)) (1/E) ∂/∂E ( E k ∂T/∂E ) + (E R'/R) ∂T/∂E
    质:  ∂C/∂τ = (1/R²) (1/E) ∂/∂E ( E D ∂C/∂E ) + (E R'/R) ∂C/∂E

    其中对流项 (E R'/R)·∂/∂E 是坐标变换（移动网格）带来的：
        ∂C/∂t|_r = ∂C/∂τ|_E - (E R'/R) ∂C/∂E  （链式法则，E = r/R(t)）
    代入物理方程 ∂C/∂t|_r = (1/r)∂/∂r(r D ∂C/∂r) 得
        ∂C/∂τ|_E = (1/R²)(1/E)∂/∂E(E D ∂C/∂E) + (E R'/R) ∂C/∂E
    注意：对流项符号为「+ E R'/R」。R'<0（收缩）时该对流把水分向中心输送，
    与「收缩使水分浓缩、从而减缓浓度下降」的物理直觉一致（见 p4_verify 校验）。

定解条件（在 E 域）：
    E = 0（中心）:  ∂T/∂E = ∂C/∂E = 0                          （轴对称）
    E = 1（表面 r=R(t)）:  -k ∂T/∂r = h (T - T∞)  →  -(k/R)∂T/∂E = h (T-T∞)
                           -D ∂C/∂r = h_m (C - C∞)  →  -(D/R)∂C/∂E = h_m (C-C∞)
    初值:  T(E,0) = 28 ℃,  C(E,0) = 2.55 kg/kg

创新处理：
    (1) 精确处理坐标变换带来的对流项 (E R'/R)∂/∂E（对 T、C 两方程一致），
        内部用中心差分（二阶）、表面用单侧差分；其符号由链式法则严格导出。
    (2) R(t) 由附件 2 数值微分得到：中心差分（np.gradient，端点单侧），线性插值；
        R'(t) 在数据末端（收缩停止后）自然趋于 0。
    (3) 网格随边界等比例缩放：E 域均匀网格 → 物理网格 r_j(t)=E_j·R(t) 同步缩放，
        全程相对分辨率一致，且守恒型有限体积权重 Σ_j V_j = (1²-0²)/2 = 1/2 恒成立。
    (4) 附录 4 变参数物性（ρ、c_p、k、D 均随 C、T 变化）在移动域上按
        Crank--Nicolson（θ=0.5）+ Picard/Newton 线性化隐式求解。

依赖：numpy、openpyxl（读取附件 2）、scipy（仅输出降采样 PCHIP）。
"""

import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import openpyxl

# ----------------------------------------------------------------------------
# 物理参数（问题 4：附录 4 变参数物性；h、h_m 沿用问题 1/2 口径）
# ----------------------------------------------------------------------------
H_HEAT = 25.0        # W/(m^2 K)  对流换热系数
H_MASS = 8.0e-7      # m/s        对流传质系数

T0 = 28.0            # degC  初始温度
C0 = 2.55            # kg/kg 初始水分浓度（干基含水率）
C_TARGET = 0.15      # kg/kg 烘干判据（各处均低于该值）

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATT = os.path.join(BASE, "A题", "附件")


# ----------------------------------------------------------------------------
# 附录 4 变参数经验公式
# ----------------------------------------------------------------------------
def rho(C):
    C = np.asarray(C, dtype=float)
    return 760.0 + 90.0 * C


def cp(C):
    C = np.asarray(C, dtype=float)
    return 1850.0 + 2150.0 * C / (C + 1.0)


def kcond(C):
    C = np.asarray(C, dtype=float)
    return 0.12 + 0.20 * C / (C + 1.0)


def diff_coef(C, T):
    """D(C,T) = 4.2e-4 exp(-0.30/C) exp(-3850/T_K)，T 单位 ℃。"""
    C = np.asarray(C, dtype=float)
    T = np.asarray(T, dtype=float)
    Tk = T + 273.15
    return 4.2e-4 * np.exp(-0.30 / np.clip(C, 1e-9, None)) * np.exp(-3850.0 / Tk)


def dD_dC(C, T):
    """∂D/∂C = D · 0.30/C^2（Newton 线性化所需敏度系数）。"""
    C = np.asarray(C, dtype=float)
    return diff_coef(C, T) * 0.30 / np.clip(C, 1e-9, None) ** 2


# ----------------------------------------------------------------------------
# Thomas 追赶法
# ----------------------------------------------------------------------------
def thomas(a, b, c, d):
    n = b.size
    psi = np.empty(n); phi = np.empty(n)
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


def _apply(op, u):
    """三对角算子作用：L u = a_j u_{j-1} + b_j u_j + c_j u_{j+1}。"""
    a, b, c = op
    v = b * u
    v[1:] += a[1:] * u[:-1]
    v[:-1] += c[:-1] * u[1:]
    return v


def _step(u, op_new, op_old, s_old, s_new, dt, theta):
    """线性 θ 格式单步（允许新旧时间层算子不同）。"""
    rhs = u + (1.0 - theta) * dt * _apply(op_old, u)
    rhs = rhs + dt * (theta * s_new + (1.0 - theta) * s_old)
    an, bn, cn = op_new
    lo = -theta * dt * an; lo[0] = 0.0
    di = 1.0 - theta * dt * bn
    up = -theta * dt * cn; up[-1] = 0.0
    return thomas(lo, di, up, rhs)


# ----------------------------------------------------------------------------
# 附件 2：半径 R(t) 及其导数 R'(t)
# ----------------------------------------------------------------------------
def load_R(path=None):
    """读取附件 2（时间 s、半径 cm），半径换算为 m。"""
    if path is None:
        path = os.path.join(ATT, "附件2.xlsx")
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Sheet1"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    t = np.array([float(r[0]) for r in rows])
    R = np.array([float(r[1]) for r in rows]) * 0.01   # cm -> m
    return t, R


def make_R(path=None):
    """
    返回 (R_fn, Rp_fn, t_data, R_data, Rp_data)。
    R_fn  : 线性插值的 R(t)（m）；
    Rp_fn : 中心差分（np.gradient，端点单侧）后线性插值的 R'(t)（m/s）。
    """
    t, R = load_R(path)
    Rp = np.gradient(R, t)               # 中心差分，端点单侧

    def R_fn(tv):
        return float(np.interp(tv, t, R))

    def Rp_fn(tv):
        return float(np.interp(tv, t, Rp))

    return R_fn, Rp_fn, t, R, Rp


# ----------------------------------------------------------------------------
# 环境条件（附件 1，复用问题 2 的读取）
# ----------------------------------------------------------------------------
def load_ambient(path=None):
    if path is None:
        path = os.path.join(ATT, "附件1.xlsx")
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Sheet1"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    t = np.array([float(r[0]) for r in rows])
    T = np.array([float(r[1]) for r in rows])
    C = np.array([float(r[2]) for r in rows])
    return t, T, C


def make_ambient(path=None):
    """返回 (T_inf_fn, C_inf_fn)。附件 1 覆盖 0--14400 s，其后维持末值。"""
    t, T, C = load_ambient(path)

    def T_inf(tv):
        return float(np.interp(tv, t, T))

    def C_inf(tv):
        return float(np.interp(tv, t, C))

    return T_inf, C_inf


# ----------------------------------------------------------------------------
# E 域均匀网格（固定域 [0,1]；物理网格随 R(t) 等比例缩放）
# ----------------------------------------------------------------------------
def make_grid(N):
    """E 域有限体积网格。Σ_j V_j = (1² - 0²)/2 = 1/2 恒成立（离散质量守恒恒等式）。"""
    dE = 1.0 / N
    E = np.arange(N + 1) / N
    rl = (np.arange(N + 1) - 0.5) / N      # 左界面
    rl[0] = 0.0
    rr = np.empty(N + 1)                   # 右界面
    rr[:N] = rl[1:]
    rr[N] = 1.0
    V = (rr ** 2 - rl ** 2) / 2.0          # 控制体积权重 ∫ E dE
    dl = np.empty(N + 1); dl[0] = 0.0; dl[1:] = E[1:] - E[:-1]
    dr = np.empty(N + 1); dr[N] = 0.0; dr[:N] = E[1:] - E[:-1]
    return dict(E=E, rl=rl, rr=rr, V=V, dl=dl, dr=dr, N=N, dE=dE)


# ----------------------------------------------------------------------------
# 坐标变换带来的对流项  (E R'/R) ∂/∂E
# ----------------------------------------------------------------------------
def op_conv(g, R, Rp):
    """
    对流算子（三对角），对 T、C 一致（系数均为 E R'/R）。
    内部：中心差分（二阶）；表面 E=1：单侧差分；中心 E=0：系数为 0。
    """
    N = g["N"]; dE = g["dE"]; E = g["E"]
    vE = E * (Rp / R)                     # E 空间的对流速度
    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)
    j = np.arange(1, N)
    a[j] = -vE[j] / (2.0 * dE)
    c[j] = +vE[j] / (2.0 * dE)
    a[N] = -vE[N] / dE                    # 表面单侧：-(R'/R)/dE · C_{N-1}
    b[N] = +vE[N] / dE                    #               +(R'/R)/dE · C_N
    return a, b, c


# ----------------------------------------------------------------------------
# 热传导算子（线性于 T，系数冻结于 C_f；含对流 + 表面 Robin 源 gT）
# ----------------------------------------------------------------------------
def op_heat(g, C_f, R, Rp):
    N = g["N"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    inv = 1.0 / (rho(C_f) * cp(C_f))      # 1/(ρ c_p)
    k = kcond(C_f)
    kh = 0.5 * (k[:-1] + k[1:])           # 界面 k_{j+1/2}
    invR2 = 1.0 / (R * R)

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    # 中心 j=0：L'Hôpital（E→0），dT/dt = (1/(ρc_p R²))·4 k_{1/2}/ΔE² (T1-T0)
    b[0] = -invR2 * inv[0] * kh[0] * rr[0] / (dr[0] * V[0])
    c[0] = +invR2 * inv[0] * kh[0] * rr[0] / (dr[0] * V[0])

    # 内部点
    j = np.arange(1, N)
    A = invR2 * inv[j] * kh[j - 1] * rl[j] / (dl[j] * V[j])
    Cc = invR2 * inv[j] * kh[j] * rr[j] / (dr[j] * V[j])
    a[1:N] = A; c[1:N] = Cc; b[1:N] = -(A + Cc)

    # 表面 j=N：半控制体积 + Robin（-k ∂T/∂r = h(T-T∞)）
    a[N] = invR2 * inv[N] * kh[N - 1] * rl[N] / (dl[N] * V[N])
    gT = H_HEAT * inv[N] / (R * V[N])
    b[N] = -a[N] - gT

    # 对流项
    ac, bc, cc = op_conv(g, R, Rp)
    a += ac; b += bc; c += cc
    return a, b, c, gT


# ----------------------------------------------------------------------------
# 传质算子：扩散（非线性 D(C,T)）+ 对流 + 表面流出
# ----------------------------------------------------------------------------
def op_mass_diff(g, C, T, R, Rp):
    """L(C) = 扩散 + 对流 + 表面流出项 -gC·C_N（不含 C_inf 源，源另加）。"""
    N = g["N"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C, T)
    Dh = 0.5 * (D[:-1] + D[1:])
    invR2 = 1.0 / (R * R)
    gC = H_MASS / (R * V[N])

    out = np.zeros(N + 1)
    out[0] = invR2 * Dh[0] * rr[0] * (C[1] - C[0]) / (dr[0] * V[0])
    j = np.arange(1, N)
    out[j] = invR2 * (Dh[j - 1] * rl[j] * (C[j - 1] - C[j]) / dl[j]
                      + Dh[j] * rr[j] * (C[j + 1] - C[j]) / dr[j]) / V[j]
    out[N] = invR2 * Dh[N - 1] * rl[N] * (C[N - 1] - C[N]) / (dl[N] * V[N]) - gC * C[N]

    ac, bc, cc = op_conv(g, R, Rp)
    out += _apply((ac, bc, cc), C)
    return out


def op_mass_jac(g, C_star, T_star, R, Rp):
    """∂L/∂C 在 C* 处（Taylor 一阶线性化，含 ∂D/∂C 敏度项 + 线性对流项）。"""
    N = g["N"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C_star, T_star)
    Dp = dD_dC(C_star, T_star)
    Dh = 0.5 * (D[:-1] + D[1:])
    invR2 = 1.0 / (R * R)
    gC = H_MASS / (R * V[N])

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    d0 = C_star[1] - C_star[0]
    fac0 = invR2 * rr[0] / (dr[0] * V[0])
    b[0] = fac0 * (-Dh[0] + 0.5 * Dp[0] * d0)
    c[0] = fac0 * (+Dh[0] + 0.5 * Dp[1] * d0)

    j = np.arange(1, N)
    dm = C_star[j - 1] - C_star[j]
    dp = C_star[j + 1] - C_star[j]
    al = invR2 * rl[j] / (dl[j] * V[j])
    ar = invR2 * rr[j] / (dr[j] * V[j])
    a[1:N] = al * (Dh[j - 1] + 0.5 * Dp[j - 1] * dm)
    c[1:N] = ar * (Dh[j] + 0.5 * Dp[j + 1] * dp)
    b[1:N] = al * (-Dh[j - 1] + 0.5 * Dp[j] * dm) + ar * (-Dh[j] + 0.5 * Dp[j] * dp)

    dN = C_star[N - 1] - C_star[N]
    fN = invR2 * rl[N] / (dl[N] * V[N])
    a[N] = fN * (Dh[N - 1] + 0.5 * Dp[N - 1] * dN)
    b[N] = fN * (-Dh[N - 1] + 0.5 * Dp[N] * dN) - gC

    ac, bc, cc = op_conv(g, R, Rp)
    a += ac; b += bc; c += cc
    return a, b, c


# ----------------------------------------------------------------------------
# 单个 Crank--Nicolson 步（弱耦合交替 + Picard 迭代）
# ----------------------------------------------------------------------------
def step_cn(g, T, C, t, dt, R_fn, Rp_fn, T_inf_fn, C_inf_fn, theta=0.5,
            picard_max=4, picard_tol=1e-10):
    N = g["N"]
    R_o, Rp_o = R_fn(t), Rp_fn(t)
    R_n, Rp_n = R_fn(t + dt), Rp_fn(t + dt)
    Ti_n, Ti_o = T_inf_fn(t + dt), T_inf_fn(t)
    Ci_n, Ci_o = C_inf_fn(t + dt), C_inf_fn(t)

    T_old, C_old = T, C
    C_star, T_star = C.copy(), T.copy()
    zero = np.zeros(N + 1)

    opT_o = op_heat(g, C_old, R_o, Rp_o)
    sT_o = zero.copy(); sT_o[N] = opT_o[3] * Ti_o

    T_new, C_new = T_old, C_old
    for _ in range(picard_max):
        opT_n = op_heat(g, C_star, R_n, Rp_n)
        sT_n = zero.copy(); sT_n[N] = opT_n[3] * Ti_n
        T_new = _step(T_old, opT_n[:3], opT_o[:3], sT_o, sT_n, dt, theta)

        gC_n = H_MASS / (R_n * g["V"][N])
        gC_o = H_MASS / (R_o * g["V"][N])
        sC_n = zero.copy(); sC_n[N] = gC_n * Ci_n
        sC_o = zero.copy(); sC_o[N] = gC_o * Ci_o

        J = op_mass_jac(g, C_star, T_new, R_n, Rp_n)
        Lstar = op_mass_diff(g, C_star, T_new, R_n, Rp_n)
        Lold = op_mass_diff(g, C_old, T_old, R_o, Rp_o)
        JC = _apply(J, C_star)
        rhs = C_old + dt * (theta * (Lstar - JC + sC_n)
                            + (1.0 - theta) * (Lold + sC_o))
        a, b, c = J
        lo = -theta * dt * a; lo[0] = 0.0
        di = 1.0 - theta * dt * b
        up = -theta * dt * c; up[-1] = 0.0
        C_new = thomas(lo, di, up, rhs)

        dT = float(np.abs(T_new - T_star).max())
        dC = float(np.abs(C_new - C_star).max())
        T_star, C_star = T_new, C_new
        if dT < picard_tol and dC < picard_tol:
            break
    return T_new, C_new


# ----------------------------------------------------------------------------
# 主求解器
# ----------------------------------------------------------------------------
def solve(T_inf_fn, C_inf_fn, R_fn, Rp_fn, N=800, dt=60.0, t_end=300000.0,
          theta=0.5, picard_max=4, picard_tol=1e-10, stop_cmax=C_TARGET,
          t_start=0.0, T_init=None, C_init=None,
          dt_startup=2.0, t_startup=3600.0):
    """
    推进求解。返回 dict(times, g, T, C, R_t, Rp_t)。
    T/C 的每一行对应一个输出时刻（间隔 dt，含 t_start 的初值行）。

    初始瞬变（t < t_startup）的每个输出步内部用 dt_startup 细分：
    表面边界层的形成在数十秒内完成，Crank--Nicolson（θ=0.5）用大步长 Δt
    会在此产生非物理振荡（表面 C 在首分钟内振荡 ~0.12 kg/kg）；用 2 s 子步
    启动即可消除，中心浓度不受影响（t_dry 由中心控制）。
    """
    g = make_grid(N)
    T = np.full(N + 1, T0 if T_init is None else float(T_init))
    C = np.full(N + 1, C0 if C_init is None else float(C_init))

    times = [t_start]
    T_hist = [T.copy()]
    C_hist = [C.copy()]
    R_hist = [R_fn(t_start)]
    Rp_hist = [Rp_fn(t_start)]

    t = float(t_start)
    while t < t_end - 1e-9:
        t_new = t + dt
        if t_new <= t_startup + 1e-9:
            T, C = march(g, T, C, t, t_new, R_fn, Rp_fn, T_inf_fn, C_inf_fn,
                         h_max=dt_startup, theta=theta, picard_max=picard_max,
                         picard_tol=picard_tol)
        else:
            T, C = step_cn(g, T, C, t, dt, R_fn, Rp_fn, T_inf_fn, C_inf_fn,
                           theta, picard_max, picard_tol)
        t = t_new
        times.append(t)
        T_hist.append(T.copy())
        C_hist.append(C.copy())
        R_hist.append(R_fn(t))
        Rp_hist.append(Rp_fn(t))
        if stop_cmax is not None and C.max() < stop_cmax:
            break
    return dict(times=np.array(times), g=g, T=np.array(T_hist), C=np.array(C_hist),
                R_t=np.array(R_hist), Rp_t=np.array(Rp_hist))


def march(g, T, C, t, t_to, R_fn, Rp_fn, T_inf_fn, C_inf_fn,
          h_max=60.0, theta=0.5, picard_max=4, picard_tol=1e-10):
    """从状态 (T,C,t) 小步推进到 t_to，用于终止时间二分。"""
    t = float(t)
    T, C = T.copy(), C.copy()
    while t < t_to - 1e-9:
        dt = min(h_max, t_to - t)
        T, C = step_cn(g, T, C, t, dt, R_fn, Rp_fn, T_inf_fn, C_inf_fn,
                       theta, picard_max, picard_tol)
        t += dt
    return T, C


# ----------------------------------------------------------------------------
# 顶层接口
# ----------------------------------------------------------------------------
def solve_problem4(N=800, dt=60.0, t_end=300000.0, verbose=True,
                   T_inf_fn=None, C_inf_fn=None, R_fn=None, Rp_fn=None):
    """完整求解问题 4，返回 dict。"""
    if T_inf_fn is None or C_inf_fn is None:
        T_inf_fn, C_inf_fn = make_ambient()
    if R_fn is None or Rp_fn is None:
        R_fn, Rp_fn, _, _, _ = make_R()

    res = solve(T_inf_fn, C_inf_fn, R_fn, Rp_fn, N=N, dt=dt, t_end=t_end)
    t_dry = float(res["times"][-1])          # 首个 max C < 0.15 的输出时刻
    if verbose:
        print("[问题4] N=%d, dt=%.0f s, 达标时刻 t_dry = %.1f s = %.4f h"
              % (N, dt, t_dry, t_dry / 3600.0))
    return dict(res=res, t_dry=t_dry, T_inf_fn=T_inf_fn, C_inf_fn=C_inf_fn,
                R_fn=R_fn, Rp_fn=Rp_fn)


if __name__ == "__main__":
    import time
    np.set_printoptions(precision=6, suppress=True)
    t0 = time.time()
    out = solve_problem4(N=200, dt=60.0)
    res = out["res"]
    print("最终 中心 C=%.6f  表面 C=%.6f  表面 R=%.4f cm"
          % (res["C"][-1, 0], res["C"][-1, -1], res["R_t"][-1] * 100))
    print("耗时 %.1f s" % (time.time() - t0))
