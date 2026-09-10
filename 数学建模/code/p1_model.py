# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 1
一维非稳态热质传递耦合模型 (Crank-Nicolson 隐式有限差分)

控制方程（无限长圆柱假设，径向一维）:
    热:  rho*cp * dT/dt = (1/r) * d/dr ( k * r * dT/dr )
    质:  dC/dt          = (1/r) * d/dr ( D(C) * r * dC/dr )

定解条件:
    初值: T(r,0)=28 degC,  C(r,0)=2.55 kg/kg
    r=0 : dT/dr = dC/dr = 0                       （轴对称）
    r=R : -k dT/dr|_R = h  (T_R - T_inf(t))       （对流换热 Robin）
          -D dC/dr|_R = hm (C_R - C_inf(t))       （对流传质 Robin）

数值格式:
    * Crank-Nicolson 隐式格式（时间二阶、无条件稳定）
    * r=0 处用 L'Hopital 法则消去 1/r 奇异项
    * r=R 处 Robin 边界精确离散（二阶精度、离散守恒）
    * D(C) 采用滞后迭代（不动点 / Picard 迭代，用上一轮 C 计算扩散系数）
    * Thomas 追赶法求解三对角方程组

边界离散方式 `bc`:
    'halfcell' —— 表面用半控制体积精确平衡（默认）：Robin 通量直接入账，
                  离散格式严格守恒，权 Σw_j = R^2/2 精确等于圆柱截面积。
    'ghost'    —— 经典虚拟节点法：在 r=R 之外引入虚拟点，用中心差分表达
                  Robin 条件。实现简单，但表面控制体积含虚构区域，不守恒。
"""

import numpy as np

# ----------------------------------------------------------------------------
# 物理参数（附录 2）
# ----------------------------------------------------------------------------
RHO = 820.0          # kg/m^3      药材密度
CP = 2600.0          # J/(kg K)    比热容
K_COND = 0.36        # W/(m K)     热传导系数
H_HEAT = 25.0        # W/(m^2 K)   对流换热系数
H_MASS = 8.0e-7      # m/s         对流传质系数

R_CYL = 0.02         # m           圆柱半径 2 cm
L_CYL = 0.25         # m           圆柱长度 25 cm

T0 = 28.0            # degC        初始温度
C0 = 2.55            # kg/kg       初始水分浓度（干基含水率）

ALPHA = K_COND / (RHO * CP)      # 热扩散系数 m^2/s = 1.6886e-7

BC_DEFAULT = "halfcell"


def diff_coef(C):
    """水分扩散系数经验公式（附录 2）: D = 7e-9 * exp(-0.89 / C)   [m^2/s]"""
    C = np.asarray(C, dtype=float)
    return 7.0e-9 * np.exp(-0.89 / np.clip(C, 1e-12, None))


# ----------------------------------------------------------------------------
# Thomas 追赶法
# ----------------------------------------------------------------------------
def thomas(a, b, c, d):
    """
    求解三对角方程组 M x = d。
    a: 次对角 (a[0] 不用)；b: 主对角；c: 超对角 (c[-1] 不用)；d: 右端项。
    """
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


# ----------------------------------------------------------------------------
# 空间算子的三对角系数  dT/dt = L_T T + gT*T_inf ;  dC/dt = L_C C + gC*C_inf
# ----------------------------------------------------------------------------
def build_op_heat(dr, N, R=R_CYL, bc=BC_DEFAULT):
    """热传导算子：返回 (a, b, c, gT)。源项仅作用于第 N 行：src_N = gT*T_inf。"""
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr

    a = np.zeros(N + 1)
    b = np.zeros(N + 1)
    c = np.zeros(N + 1)

    # ---- 中心 r=0：L'Hopital 法则 ----
    #   lim_{r->0} (1/r) d/dr(alpha*r*dT/dr) = 2*alpha*T''(0)
    #   配合对称虚拟点 T_{-1}=T_1 的二阶差分：
    #   dT/dt|_0 = 4*alpha*(T_1 - T_0)/dr^2
    b[0] = -4.0 * ALPHA / dr ** 2
    c[0] = 4.0 * ALPHA / dr ** 2

    # ---- 内点 j=1..N-1：守恒型有限体积（二阶中心差分） ----
    for j in range(1, N):
        A = ALPHA * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = ALPHA * r_half[j] / (r[j] * dr ** 2)
        a[j] = A
        c[j] = Cc
        b[j] = -(A + Cc)

    if bc == "halfcell":
        # ---- 表面 r=R：半控制体积精确平衡 ----
        # ∫_{r_{N-1/2}}^{R} T_t r dr = [alpha r T_r]_{r_{N-1/2}}^{R}
        # 且 alpha*T_r|_R = -h(T_N - T_inf)/(rho cp)，控制体积 V_N = ∫_{r_{N-1/2}}^R r dr
        VN = R * dr / 2.0 - dr ** 2 / 8.0
        a[N] = ALPHA * r_half[N - 1] / (VN * dr)
        b[N] = -a[N] - R * H_HEAT / (RHO * CP * VN)
        c[N] = 0.0
        gT = R * H_HEAT / (RHO * CP * VN)
    elif bc == "ghost":
        # ---- 虚拟节点法（对照用） ----
        beta = dr * H_HEAT / K_COND
        r_half_N = R + 0.5 * dr
        a[N] = 2.0 * ALPHA / dr ** 2
        b[N] = -ALPHA * (r_half_N * (1.0 + 2.0 * beta) + r_half[N - 1]) / (R * dr ** 2)
        c[N] = 0.0
        gT = 2.0 * H_HEAT * r_half_N / (RHO * CP * R * dr)
    else:
        raise ValueError(bc)
    return a, b, c, gT


def build_op_mass(dr, N, C_lag, R=R_CYL, bc=BC_DEFAULT):
    """
    传质算子 L_C：返回 (a, b, c, gC)，src_N = gC*C_inf。
    D(C) 滞后：全部扩散系数由上一轮迭代的 C_lag 计算，使算子在本轮内为常系数。
    """
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr

    D = diff_coef(C_lag)                    # 逐节点扩散系数（滞后）
    D_half = 0.5 * (D[:-1] + D[1:])         # 界面系数取算术平均

    a = np.zeros(N + 1)
    b = np.zeros(N + 1)
    c = np.zeros(N + 1)

    # ---- 中心 r=0：L'Hopital 法则 ----
    b[0] = -4.0 * D_half[0] / dr ** 2
    c[0] = 4.0 * D_half[0] / dr ** 2

    # ---- 内点 ----
    for j in range(1, N):
        A = D_half[j - 1] * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = D_half[j] * r_half[j] / (r[j] * dr ** 2)
        a[j] = A
        c[j] = Cc
        b[j] = -(A + Cc)

    if bc == "halfcell":
        VN = R * dr / 2.0 - dr ** 2 / 8.0
        a[N] = D_half[N - 1] * r_half[N - 1] / (VN * dr)
        b[N] = -a[N] - R * H_MASS / VN
        c[N] = 0.0
        gC = R * H_MASS / VN
    elif bc == "ghost":
        DN = D[N]
        gamma = dr * H_MASS / DN
        r_half_N = R + 0.5 * dr
        a[N] = 2.0 * DN / dr ** 2
        b[N] = -DN * (r_half_N * (1.0 + 2.0 * gamma) + r_half[N - 1]) / (R * dr ** 2)
        c[N] = 0.0
        gC = 2.0 * H_MASS * r_half_N / (R * dr)
    else:
        raise ValueError(bc)
    return a, b, c, gC


# ----------------------------------------------------------------------------
# 守恒权重：使得 Σ_j w_j * dC_j/dt 只剩边界通量项（严格离散守恒）
# ----------------------------------------------------------------------------
def conservative_weights(dr, N, R=R_CYL):
    """w_0 = dr^2/8, w_j = r_j dr (1<=j<=N-1), w_N = R dr/2 - dr^2/8。Σw_j = R^2/2。"""
    w = np.arange(N + 1) * dr * dr
    w[0] = dr ** 2 / 8.0
    w[N] = R * dr / 2.0 - dr ** 2 / 8.0
    return w


# ----------------------------------------------------------------------------
# theta 格式单步
# ----------------------------------------------------------------------------
def _shift(u, k):
    v = np.zeros_like(u)
    if k == -1:
        v[1:] = u[:-1]
    else:
        v[:-1] = u[1:]
    return v


def _step(u, op_new, op_old, s_old, s_new, dt, theta):
    """
    theta 格式单步（允许新旧时间层使用不同的空间算子）：

        (I - theta*dt*L_new) u^{n+1} = (I + (1-theta)*dt*L_old) u^n
                                       + dt*(theta*s_new + (1-theta)*s_old)

    对常系数问题 L_old = L_new；对变扩散系数 D(C) 的传质方程，
    必须 L_old = L(D^n)、L_new = L(D^{n+1}) —— 否则 O(dt^2) 的局部误差
    会累积为 O(dt) 的全局误差，使 CN 退化为时间一阶。
    """
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
# 主求解器
# ----------------------------------------------------------------------------
def solve(T_inf_fn, C_inf_fn, t_end=1800.0, dt=1.0, N=20, R=R_CYL,
          theta=0.5, n_store=None, bc=BC_DEFAULT,
          picard_max=3, picard_tol=1e-10, T_init=None, C_init=None):
    """
    求解一维非稳态热质耦合模型。

    参数
    ----
    T_inf_fn, C_inf_fn : callable(t)->float  烘房温度 / 水分浓度
    t_end              : 终止时间 (s)
    dt                 : 时间步长 (s)
    N                  : 径向网格数（dr = R/N）
    theta              : 0.5 = Crank-Nicolson；1.0 = 向后 Euler；0.0 = 显式
    n_store            : 每隔多少步记录一次；None 表示只记录终态
    bc                 : 'halfcell'（默认，守恒）或 'ghost'
    picard_max         : D(C) 滞后迭代的最大次数。
                         1 = 单次滞后（显式滞后）；>1 = Picard 不动点迭代直至收敛。
    picard_tol         : Picard 迭代收敛容差（K 或 kg/kg 的无穷范数）

    返回
    ----
    times, r, T_hist, C_hist
    """
    dr = R / N
    r = np.arange(N + 1) * dr
    nsteps = int(round(t_end / dt))

    T = np.full(N + 1, T0 if T_init is None else float(T_init))
    C = np.full(N + 1, C0 if C_init is None else float(C_init))

    ah, bh, ch, gT = build_op_heat(dr, N, R, bc)

    times, T_hist, C_hist = [], [], []
    if n_store:
        times.append(0.0); T_hist.append(T.copy()); C_hist.append(C.copy())

    zero = np.zeros(N + 1)

    for n in range(1, nsteps + 1):
        t_new, t_old = n * dt, (n - 1) * dt
        Ti_n, Ti_o = T_inf_fn(t_new), T_inf_fn(t_old)
        Ci_n, Ci_o = C_inf_fn(t_new), C_inf_fn(t_old)

        # ---------------- 温度场（常系数，L_old = L_new） ----------------
        sT_n = zero.copy(); sT_n[N] = gT * Ti_n
        sT_o = zero.copy(); sT_o[N] = gT * Ti_o
        T = _step(T, (ah, bh, ch), (ah, bh, ch), sT_o, sT_n, dt, theta)

        # ---------------- 水分浓度场（D(C) 滞后 / Picard 迭代） ----------------
        C_old = C.copy()                 # 时间步初值 C^n
        C_it = C.copy()                  # 迭代解，兼作 D^{n+1} 的滞后取值

        # 旧时间层算子：D 取 C^n（不随迭代变化，只需组装一次）
        am_o, bm_o, cm_o, gC_o = build_op_mass(dr, N, C_old, R, bc)
        sC_o = zero.copy(); sC_o[N] = gC_o * Ci_o

        for _ in range(picard_max):
            # 新时间层算子：D 取当前迭代解 C_it
            am_n, bm_n, cm_n, gC_n = build_op_mass(dr, N, C_it, R, bc)
            sC_n = zero.copy(); sC_n[N] = gC_n * Ci_n
            C_new = _step(C_old, (am_n, bm_n, cm_n), (am_o, bm_o, cm_o),
                          sC_o, sC_n, dt, theta)
            done = np.abs(C_new - C_it).max() < picard_tol
            C_it = C_new
            if done:
                break
        C = C_it

        if n_store and (n % n_store == 0):
            times.append(t_new); T_hist.append(T.copy()); C_hist.append(C.copy())

    if not n_store:
        times.append(t_end); T_hist.append(T.copy()); C_hist.append(C.copy())

    return np.array(times), r, np.array(T_hist), np.array(C_hist)
