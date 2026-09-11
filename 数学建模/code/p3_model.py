# -*- coding: utf-8 -*-
"""
2026 高教社杯 A 题「药材的烘干问题」——问题 3
长时间烘干过程模拟 + 自适应步长 + 终止时间精确判定

问题 3 要求：药材**各处**的水分浓度低于 0.15 kg/kg，确定烘干所需时间（h）。

控制方程、物性经验公式、定解条件与问题 2 完全一致（附录 3），本模块在其之上
新增五项数值处理（前四项对应题解要求，第五项为解决长时间模拟特有的网格刚性问题）：

  (1) 多时间尺度自适应步长
      宏观步长 h = k·Δt_out，Δt_out = 60 s（输出节拍）。预热平衡阶段 k≡1（h = 60 s），
      进入恒温干燥后随误差反馈逐步放大到 k = K_MAX（默认 5，即 h = 300 s）。
      误差用**步长加倍法**（Richardson：m 个子步 vs m/2 个子步）估计；超差时先加密
      子步 m ← 2m（上限 M_MAX），仍不达标再回退宏观步长 k ← ⌊k/2⌋。
      宏观步长恒为 Δt_out 的整数倍，故宏观步端点精确落在输出网格上。

  (2) 稠密输出（dense output）
      宏观步内部由各子步锚点的值与端导数作分段三次 Hermite 插值，补齐落在步内的
      60 s 输出点。误差 O(h_panel⁴)，在 p3_verify.py V2 中以"从同一起点直接积分到
      该时刻"作定向隔离测试。

  (3) 终止时间二分搜索
      粗定位：宏观步推进，记录 max_r C 首次跌破 0.15 的区间 [t_a, t_b]（宽度 ≤ 300 s）；
      细定位：自 t_a 的检查点起二分积分，收敛到 < 1 s；报出时按"分钟向上取整"给出
      保守口径（第一个确实各处都达标的分钟）。

  (4) 准稳态（quasi-steady, QS）近似与外推
      后段剖面趋于自相似 C(r,t) ≈ φ(r)·C(0,t)。由质量平衡 dM/dt = -R h_m (C_N - C_inf)
      得两条外推：
        QS-0 冻结形状：M(t) = M_eq + (M_0 - M_eq)e^{-a(t-t_0)}，a = R h_m φ_N/(W φ̄)；
        QS-1 形状漂移修正：把 φ̄(C_0)、φ_N(C_0) 由已算轨道拟合成多项式后，把 PDE 降为
             一维标量 ODE  dC_0/dt = -(R h_m/W)(φ_N C_0 - C_inf)/(φ̄ + C_0 φ̄′)。
      QS-0 是纯解析式；QS-1 需要历史轨道标定，二者共同构成对二分的独立校核。

  (5) r = R 附近的分级（边界层加密）网格 —— 长时间模拟的网格刚性问题
      干燥末期 D(C) 在径向上可变化两个量级以上（本例 272 倍），使表面出现极薄的
      扩散边界层：t = 50 h 时 C 在最后 0.003 cm 内由 0.078 掉到 0.053。均匀网格必须
      把 Δr 降到 1e-5 m 量级才能分辨，收敛缓慢（t_dry 的表观空间阶仅约 1.5）。
      改用分级坐标 ξ → r = R[1-(1-ξ)^q]（q > 1，单元尺寸在 r = R 处按 (1-ξ)^{q-1}
      衰减），可用很小的 N 分辨该薄层。q = 1 退化为均匀网格，故本实现与问题 1/2 的
      格式严格兼容（见 p3_verify.py V10 的退化一致性检查）。

所有守恒型有限体积算子（中心 L'Hôpital、表面半控制体积 Robin、界面系数）
在分级网格上重新组装，且在 q = 1 时逐位退化为问题 2 的均匀网格算子。

依赖：numpy、scipy（仅用于输出降采样的 PCHIP 插值）。
"""

import os
import sys
# 控制台按 UTF-8 输出（Windows GBK 终端下避免 UnicodeEncodeError 与乱码）
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_model as pm
from p2_model import (R_CYL, H_HEAT, H_MASS, T0, C0, thomas,
                      rho, cp, kcond, diff_coef, dD_dC, make_ambient)

# ----------------------------------------------------------------------------
# 全局口径
# ----------------------------------------------------------------------------
OUT_DT = 60.0          # s     输出节拍（问题 3 要求每隔 60 s）
OUT_R = np.arange(21) * 1e-3   # m     输出距离网格：0, 0.1, ..., 2.0 cm
T_ENV_END = 14400.0    # s     附件 1 覆盖到 14400 s，其后恒温干燥（维持末值）
C_TARGET = 0.15        # kg/kg 烘干判据（各处均低于该值）
K_MAX = 5              # 宏观步长上限倍数（5×60 = 300 s）
M_MAX = 32             # 子步加密上限
TOL = 1e-8             # 步长加倍法的接受容差


# ============================================================================
# 网格
# ============================================================================
def make_grid(N, R=R_CYL, kind="uniform", q=1.0):
    """
    构造径向网格及其有限体积度量。

    kind = "uniform" : r_j = j·R/N（与问题 1/2 相同）
    kind = "graded"  : ξ_j = j/N,  r_j = R[1-(1-ξ_j)^q]，q > 1 时单元尺寸在 r=R 处
                       按 (1-ξ)^{q-1} 收缩；q = 1 时退化为均匀网格。

    有限体积构造：节点 r_j 取为单元中心，单元 j 的左/右界面取在 ξ 的中点
        ξ^f_j = (ξ_{j-1}+ξ_j)/2,  r^f_j = R[1-(1-ξ^f_j)^q],   r^f_0 = 0,  r^f_N = R。
    于是 Σ_j V_j = (R² - 0)/2 = R²/2 恒成立（离散质量守恒是恒等式）。
    q = 1 时 r^f_j = (j-1/2)·R/N，与问题 2 的 r_half 完全一致。

    返回 dict：r 节点、rl 左界面、rr 右界面、V 控制体积权重、dl/dr 节点间距。
    """
    if kind == "uniform" or (kind == "graded" and q == 1.0):
        r = np.arange(N + 1) * (R / N)
        rl = R * (np.arange(N + 1) - 0.5) / N
        rl[0] = 0.0
    elif kind == "graded":
        xi = np.arange(N + 1) / N
        r = R * (1.0 - (1.0 - xi) ** q)
        xif = (np.arange(N + 1) - 0.5) / N
        rl = R * (1.0 - (1.0 - xif) ** q)
        rl[0] = 0.0
    else:
        raise ValueError("未知网格类型 %r" % kind)

    rr = np.empty(N + 1)          # 右界面与左界面共用（界面唯一）
    rr[:N] = rl[1:]
    rr[N] = R

    V = (rr ** 2 - rl ** 2) / 2.0

    dl = np.empty(N + 1)      # dl[j] = r_j - r_{j-1}（到左邻点距离）
    dl[0] = 0.0
    dl[1:] = r[1:] - r[:-1]
    dr = np.empty(N + 1)      # dr[j] = r_{j+1} - r_j（到右邻点距离）
    dr[N] = 0.0
    dr[:N] = r[1:] - r[:-1]

    return dict(r=r, rl=rl, rr=rr, V=V, dl=dl, dr=dr, N=N, R=R, kind=kind, q=q)


# ============================================================================
# 守恒型有限体积算子（分级网格通用；q = 1 时与 p2_model 的均匀算子逐位相同）
# ============================================================================
def op_heat(g, C_f):
    """热传导算子 + 表面 Robin 源系数 gT（s_N = gT·T_inf）。线性于 T，系数冻结于 C_f。"""
    N, R = g["N"], g["R"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    inv = 1.0 / (rho(C_f) * cp(C_f))
    k = kcond(C_f)
    kh = 0.5 * (k[:-1] + k[1:])          # 界面 k_{j+1/2}

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    # 中心 j=0：半控制体积 + 对称条件，dT/dt = k_{1/2} r_{1/2}(T_1-T_0)/(dr_0·V_0·ρc_p)
    b[0] = -kh[0] * rr[0] * inv[0] / (dr[0] * V[0])
    c[0] = +kh[0] * rr[0] * inv[0] / (dr[0] * V[0])

    # 内部点
    j = np.arange(1, N)
    A = kh[j - 1] * rl[j] * inv[j] / (dl[j] * V[j])
    Cc = kh[j] * rr[j] * inv[j] / (dr[j] * V[j])
    a[1:N] = A; c[1:N] = Cc; b[1:N] = -(A + Cc)

    # 表面 j=N：半控制体积，外边界面（半径 R）直接代入 Robin 条件
    a[N] = kh[N - 1] * rl[N] * inv[N] / (dl[N] * V[N])
    gT = R * H_HEAT * inv[N] / V[N]
    b[N] = -a[N] - gT
    return a, b, c, gT


def op_mass_diff(g, C, T):
    """传质扩散算子 L(C)（含表面流出项 -gC·C_N，不含 C_inf 源）。"""
    N, R = g["N"], g["R"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C, T)
    Dh = 0.5 * (D[:-1] + D[1:])
    gC = R * H_MASS / V[N]

    out = np.zeros(N + 1)
    out[0] = Dh[0] * rr[0] * (C[1] - C[0]) / (dr[0] * V[0])
    j = np.arange(1, N)
    out[j] = (Dh[j - 1] * rl[j] * (C[j - 1] - C[j]) / dl[j]
              + Dh[j] * rr[j] * (C[j + 1] - C[j]) / dr[j]) / V[j]
    out[N] = Dh[N - 1] * rl[N] * (C[N - 1] - C[N]) / (dl[N] * V[N]) - gC * C[N]
    return out


def op_mass_jac(g, C_star, T_star):
    """∂L/∂C 在 C*=C_star 处的三对角（Taylor 一阶线性化，含 ∂D/∂C 敏度项）。"""
    N, R = g["N"], g["R"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C_star, T_star)
    Dp = dD_dC(C_star, T_star)
    Dh = 0.5 * (D[:-1] + D[1:])
    gC = R * H_MASS / V[N]

    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)

    # j = 0
    d0 = C_star[1] - C_star[0]
    fac0 = rr[0] / (dr[0] * V[0])
    b[0] = fac0 * (-Dh[0] + 0.5 * Dp[0] * d0)
    c[0] = fac0 * (+Dh[0] + 0.5 * Dp[1] * d0)

    # 内部
    j = np.arange(1, N)
    dm = C_star[j - 1] - C_star[j]
    dp = C_star[j + 1] - C_star[j]
    al = rl[j] / (dl[j] * V[j])
    ar = rr[j] / (dr[j] * V[j])
    a[1:N] = al * (Dh[j - 1] + 0.5 * Dp[j - 1] * dm)
    c[1:N] = ar * (Dh[j] + 0.5 * Dp[j + 1] * dp)
    b[1:N] = al * (-Dh[j - 1] + 0.5 * Dp[j] * dm) + ar * (-Dh[j] + 0.5 * Dp[j] * dp)

    # 表面
    dN = C_star[N - 1] - C_star[N]
    fN = rl[N] / (dl[N] * V[N])
    a[N] = fN * (Dh[N - 1] + 0.5 * Dp[N - 1] * dN)
    b[N] = fN * (-Dh[N - 1] + 0.5 * Dp[N] * dN) - gC
    return a, b, c


def _apply(op, u):
    """三对角算子作用：L u = a_j u_{j-1} + b_j u_j + c_j u_{j+1}。"""
    a, b, c = op
    v = b * u
    v[1:] += a[1:] * u[:-1]
    v[:-1] += c[:-1] * u[1:]
    return v


def _step(u, op_new, op_old, s_old, s_new, dt, theta):
    """θ 格式单步（允许新旧时间层算子不同）。"""
    rhs = u + (1.0 - theta) * dt * _apply(op_old, u)
    rhs = rhs + dt * (theta * s_new + (1.0 - theta) * s_old)
    a_n, b_n, c_n = op_new
    lo = -theta * dt * a_n; lo[0] = 0.0
    di = 1.0 - theta * dt * b_n
    up = -theta * dt * c_n; up[-1] = 0.0
    return thomas(lo, di, up, rhs)


# ============================================================================
# 右端项（供 Hermite 稠密输出的端点导数）
# ============================================================================
def rhs(g, T, C, t, T_inf_fn, C_inf_fn):
    N = g["N"]
    Ti, Ci = T_inf_fn(t), C_inf_fn(t)
    ah, bh, ch, gT = op_heat(g, C)
    sT = np.zeros(N + 1); sT[N] = gT * Ti
    dT = _apply((ah, bh, ch), T) + sT
    dC = op_mass_diff(g, C, T)
    dC[N] += (g["R"] * H_MASS / g["V"][N]) * Ci
    return dT, dC


# ============================================================================
# 单个 Crank--Nicolson 步（与问题 2 求解器的内部逻辑一致）
# ============================================================================
def step_cn(g, T, C, t, dt, T_inf_fn, C_inf_fn, theta=0.5,
            picard_max=3, picard_tol=1e-10, newton=True):
    N = g["N"]
    Ti_n, Ti_o = T_inf_fn(t + dt), T_inf_fn(t)
    Ci_n, Ci_o = C_inf_fn(t + dt), C_inf_fn(t)
    T_old, C_old = T, C
    C_star, T_star = C.copy(), T.copy()
    zero = np.zeros(N + 1)

    opT_o = op_heat(g, C_old)
    sT_o = zero.copy(); sT_o[N] = opT_o[3] * Ti_o

    T_new, C_new = T_old, C_old
    for _ in range(picard_max):
        opT_n = op_heat(g, C_star)
        sT_n = zero.copy(); sT_n[N] = opT_n[3] * Ti_n
        T_new = _step(T_old, opT_n[:3], opT_o[:3], sT_o, sT_n, dt, theta)

        gC = g["R"] * H_MASS / g["V"][N]
        sC_n = zero.copy(); sC_n[N] = gC * Ci_n
        sC_o = zero.copy(); sC_o[N] = gC * Ci_o

        if newton:
            J = op_mass_jac(g, C_star, T_new)
            Lstar = op_mass_diff(g, C_star, T_new)
            Lold = op_mass_diff(g, C_old, T_old)
            JC = _apply(J, C_star)
            rhs_v = C_old + dt * (theta * (Lstar - JC + sC_n)
                                  + (1.0 - theta) * (Lold + sC_o))
            a, b, c = J
            lo = -theta * dt * a; lo[0] = 0.0
            di = 1.0 - theta * dt * b
            up = -theta * dt * c; up[-1] = 0.0
            C_new = thomas(lo, di, up, rhs_v)
        else:
            opC_n = _op_mass_lag(g, C_star, T_new)
            opC_o = _op_mass_lag(g, C_old, T_old)
            C_new = _step(C_old, opC_n, opC_o, sC_o, sC_n, dt, theta)

        dT = float(np.abs(T_new - T_star).max())
        dC = float(np.abs(C_new - C_star).max())
        T_star, C_star = T_new, C_new
        if dT < picard_tol and dC < picard_tol:
            break
    return T_new, C_new


def _op_mass_lag(g, C_f, T_f):
    """冻结系数（0 阶滞后）传质算子，用于对照 Taylor 线性化的收益。"""
    N, R = g["N"], g["R"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C_f, T_f)
    Dh = 0.5 * (D[:-1] + D[1:])
    gC = R * H_MASS / V[N]
    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)
    b[0] = -Dh[0] * rr[0] / (dr[0] * V[0])
    c[0] = +Dh[0] * rr[0] / (dr[0] * V[0])
    j = np.arange(1, N)
    A = Dh[j - 1] * rl[j] / (dl[j] * V[j])
    Cc = Dh[j] * rr[j] / (dr[j] * V[j])
    a[1:N] = A; c[1:N] = Cc; b[1:N] = -(A + Cc)
    a[N] = Dh[N - 1] * rl[N] / (dl[N] * V[N])
    b[N] = -a[N] - gC
    return (a, b, c)


# ============================================================================
# 三次 Hermite 稠密输出
# ============================================================================
def hermite(u0, d0, u1, d1, h, s):
    """区间 [t, t+h] 上由端点值与端点导数确定的三次 Hermite 多项式在 s∈[0,1] 处的值。"""
    h00 = 2 * s ** 3 - 3 * s ** 2 + 1
    h10 = s ** 3 - 2 * s ** 2 + s
    h01 = -2 * s ** 3 + 3 * s ** 2
    h11 = s ** 3 - s ** 2
    return h00 * u0 + h10 * h * d0 + h01 * u1 + h11 * h * d1


# ============================================================================
# 自适应求解器
# ============================================================================
class DryingSolver:
    """可变步长 Crank--Nicolson 求解器（自适应宏观步 + 三次 Hermite 稠密输出）。"""

    def __init__(self, T_inf_fn, C_inf_fn, N=20, R=R_CYL, out_dt=OUT_DT,
                 theta=0.5, newton=True, tol=TOL, k_max=K_MAX, m_max=M_MAX,
                 kind="uniform", q=1.0):
        self.T_inf, self.C_inf = T_inf_fn, C_inf_fn
        self.N, self.R, self.out_dt, self.theta = N, R, out_dt, theta
        self.newton, self.tol, self.k_max, self.m_max = newton, tol, k_max, m_max
        self.g = make_grid(N, R, kind, q)
        self.dr = R / N
        self.r = self.g["r"]
        self.w = self.g["V"]

        self.t = 0.0
        self.T = np.full(N + 1, T0)
        self.C = np.full(N + 1, C0)

        self.k = 1          # 宏观步长倍数（h = k·out_dt）
        self.m = 2          # 宏观步内子步数
        self.nstep = 0      # 累计 CN 步数
        self.log = []       # [(t0, h, m, err, k)]

    # -- 状态检查点 ---------------------------------------------------------
    def checkpoint(self):
        return dict(t=self.t, T=self.T.copy(), C=self.C.copy())

    def restore(self, cp):
        self.t = float(cp["t"])
        self.T = cp["T"].copy()
        self.C = cp["C"].copy()

    def cmax(self):
        return float(self.C.max())

    # -- 子步推进（不改变自身状态） -----------------------------------------
    def _substep(self, h, m, cp=None):
        t = self.t if cp is None else cp["t"]
        T = (self.T if cp is None else cp["T"]).copy()
        C = (self.C if cp is None else cp["C"]).copy()
        dd = h / m
        anchors = [(t, T.copy(), C.copy())]
        for _ in range(m):
            T, C = step_cn(self.g, T, C, t, dd, self.T_inf, self.C_inf,
                           self.theta, newton=self.newton)
            t += dd
            anchors.append((t, T.copy(), C.copy()))
        self.nstep += m
        return T, C, anchors

    # -- 宏观单步：误差控制 --------------------------------------------------
    def macro_step(self, h, cp=None):
        """执行一个宏观步（长度 h）。返回 (T, C, anchors, err, m)。"""
        m = max(2, self.m)
        while True:
            T_f, C_f, anc_f = self._substep(h, m, cp)
            T_c, C_c, _ = self._substep(h, max(1, m // 2), cp)
            err = float(max(np.abs(T_f - T_c).max(), np.abs(C_f - C_c).max()))
            if err <= self.tol or m >= self.m_max:
                break
            m *= 2
        self.m = m
        return T_f, C_f, anc_f, err, m

    # -- 稠密输出：把宏观步内部的 60 s 节拍补出来 ----------------------------
    def _dense(self, anchors, t_from, t_to):
        marks = np.arange(np.floor(t_from / self.out_dt) * self.out_dt + self.out_dt,
                          t_to + 1e-9, self.out_dt)
        marks = marks[(marks > t_from + 1e-9) & (marks < t_to - 1e-9)]
        if marks.size == 0:
            return [], [], []
        der = [rhs(self.g, Ta, Ca, ta, self.T_inf, self.C_inf)
               for (ta, Ta, Ca) in anchors]
        acoord = np.array([a[0] for a in anchors])
        ts, Ts, Cs = [], [], []
        for tm in marks:
            i = int(np.searchsorted(acoord, tm, side="right")) - 1
            ta, Ta, Ca = anchors[i]
            tb, Tb, Cb = anchors[i + 1]
            hh = tb - ta
            s = (tm - ta) / hh
            dTa, dCa = der[i]
            dTb, dCb = der[i + 1]
            ts.append(tm)
            Ts.append(hermite(Ta, dTa, Tb, dTb, hh, s))
            Cs.append(hermite(Ca, dCa, Cb, dCb, hh, s))
        return np.array(ts), np.array(Ts), np.array(Cs)

    # -- 主推进 --------------------------------------------------------------
    def run(self, t_end, record=True, stop_cmax=None, k_max=None, tol=None):
        if k_max is not None:
            self.k_max = k_max
        if tol is not None:
            self.tol = tol
        t_end = float(t_end)
        rec_t, rec_T, rec_C = [self.t], [self.T.copy()], [self.C.copy()]
        stop = None

        while self.t < t_end - 1e-9:
            cp = self.checkpoint()
            k = 1 if self.t < T_ENV_END else self.k
            h = k * self.out_dt
            if self.t + h > t_end + 1e-9:
                k = int(round((t_end - self.t) / self.out_dt))
                h = k * self.out_dt
                if h <= 0:
                    break
            T_f, C_f, anchors, err, m = self.macro_step(h, cp)
            if not (np.isfinite(T_f).all() and np.isfinite(C_f).all()):
                raise RuntimeError(
                    "解在 t = %.1f s 处失去有限性（网格过度加密导致三对角矩阵退化）"
                    % cp["t"])

            self.T, self.C = T_f, C_f
            self.t = self.t + h

            if record:
                ts, Ts, Cs = self._dense(anchors, cp["t"], self.t)
                for i in range(len(ts)):
                    rec_t.append(ts[i]); rec_T.append(Ts[i]); rec_C.append(Cs[i])
                rec_t.append(self.t); rec_T.append(self.T.copy()); rec_C.append(self.C.copy())

            self.log.append(dict(t0=cp["t"], h=h, m=m, err=err, k=k))

            if err > self.tol:
                self.k = max(1, k // 2)
                self.m = max(2, m)
            else:
                if err < self.tol / 64.0 and self.m > 2:
                    self.m = max(2, self.m // 2)
                if self.t >= T_ENV_END:
                    self.k = min(self.k_max, k + 1)

            if stop_cmax is not None and self.cmax() < stop_cmax:
                stop = dict(t_a=cp["t"], t_b=self.t, cp=cp,
                            T_a=cp["T"], C_a=cp["C"],
                            T_b=self.T.copy(), C_b=self.C.copy())
                break

        return dict(t=np.array(rec_t), T=np.array(rec_T), C=np.array(rec_C),
                    r=self.r, stop=stop, log=self.log)

    # -- 从任意检查点小步推进（用于二分细定位） ------------------------------
    def march(self, cp, t_to, h_max=OUT_DT):
        t = float(cp["t"])
        T, C = cp["T"].copy(), cp["C"].copy()
        while t < t_to - 1e-9:
            dt = min(h_max, t_to - t)
            T, C = step_cn(self.g, T, C, t, dt, self.T_inf, self.C_inf,
                           self.theta, newton=self.newton)
            t += dt
            self.nstep += 1
        return T, C


# ============================================================================
# 输出降采样：把（可能是分级网格的）节点解插值到题目要求的 0.1 cm 网格
# ============================================================================
def downsample(r_src, field, r_out=OUT_R):
    """沿半径用单调三次（PCHIP）插值到输出网格；field 形状 (n, len(r_src))。"""
    from scipy.interpolate import PchipInterpolator
    field = np.atleast_2d(field)
    out = np.empty((field.shape[0], len(r_out)))
    for i in range(field.shape[0]):
        out[i] = PchipInterpolator(r_src, field[i])(r_out)
    return out


# ============================================================================
# 终止时间：二分细定位
# ============================================================================
def bisect_termination(solver, stop, target=C_TARGET, tol=1.0, h_max=OUT_DT):
    """在 [t_a, t_b] 内二分定位 max_r C(t) = target；返回 (t*, 过程记录)。"""
    t_a, t_b = float(stop["t_a"]), float(stop["t_b"])
    cp_a = dict(t=t_a, T=stop["T_a"].copy(), C=stop["C_a"].copy())
    hist = []
    while t_b - t_a > tol:
        t_m = 0.5 * (t_a + t_b)
        T_m, C_m = solver.march(cp_a, t_m, h_max=h_max)
        cmax = float(C_m.max())
        hist.append((t_m, cmax))
        if cmax < target:
            t_b = t_m
        else:
            t_a = t_m
    return t_b, hist


# ============================================================================
# 准稳态（QS）分析与外推
# ============================================================================
def qs_state(t, C, w, R=R_CYL, h_m=H_MASS, C_inf=0.0):
    """
    由守恒权重 w（Σw_j = R²/2）定义
        M = Σ w_j C_j,  W = Σ w_j,  φ_j = C_j/C_0,  φ̄ = Σ w_j φ_j / W。
    质量平衡（守恒格式的恒等式）：dM/dt = -R h_m (C_N - C_inf)。
    形状冻结时 C_N = φ_N C_0 = φ_N M/(W φ̄)，得 dM/dt = -aM + b，
        a = R h_m φ_N/(W φ̄),  b = R h_m C_inf,  M_eq = b/a。
    """
    W = float(w.sum())
    phi = C / C[0]
    phibar = float(np.dot(w, phi) / W)
    phiN = float(phi[-1])
    M = float(np.dot(w, C))
    a = R * h_m * phiN / (W * phibar)
    b = R * h_m * C_inf
    M_eq = b / a
    return dict(t=t, M=M, phi=phi, phibar=phibar, phiN=phiN, a=a, b=b,
                M_eq=M_eq, C_eq=(M_eq / (W * phibar)), W=W)


def qs_extrapolate(qs, target=C_TARGET):
    """QS-0 冻结形状解析外推到 C_0 = target；返回 (t*, 说明)。"""
    W, phibar = qs["W"], qs["phibar"]
    M_star = target * W * phibar
    M0, M_eq, a = qs["M"], qs["M_eq"], qs["a"]
    if M0 <= M_star:
        return qs["t"], "已达标"
    if M_star <= M_eq:
        return np.inf, "QS 平衡浓度高于判据，指数外推失效"
    return qs["t"] + float(np.log((M0 - M_eq) / (M_star - M_eq)) / a), "ok"


def qs_shape_fit(times, C, w, t_fit_from, deg=2):
    """把 φ̄(C_0)、φ_N(C_0) 拟合成多项式（形状随 C_0 缓慢漂移，故以 C_0 为自变量）。"""
    idx = np.where(times >= t_fit_from - 1e-9)[0]
    idx = idx[idx < len(C)]
    C0 = C[idx, 0]
    phi = C[idx] / C0[:, None]
    pb = np.polyfit(C0, phi @ w / w.sum(), deg)
    pn = np.polyfit(C0, phi[:, -1], deg)
    return pb, pn, float(C0.min()), float(C0.max()), float(times[idx[0]])


def qs_drift_rhs(C0, pb, pn, W, R=R_CYL, h_m=H_MASS, C_inf=0.0):
    """把 PDE 降成一维标量 ODE：dC_0/dt = -(R h_m/W)(φ_N C_0 - C_inf)/(φ̄ + C_0 φ̄′)。"""
    P = np.polyval(pb, C0)
    Q = np.polyval(pn, C0)
    dP = np.polyval(np.polyder(pb), C0)
    return -(R * h_m / W) * (Q * C0 - C_inf) / (P + C0 * dP)


def qs_extrapolate_drift(C0_now, t_now, pb, pn, W, target=C_TARGET,
                         R=R_CYL, h_m=H_MASS, C_inf=0.0, dt=60.0, t_max=6e5):
    """RK4 积分形状漂移修正的一维 QS 方程直到 C_0 = target；返回 (t*, 轨迹)。"""
    C0, t = float(C0_now), float(t_now)
    tr = [(t, C0)]
    if C0 <= target:
        return t, tr
    while C0 > target and t < t_max:
        k1 = qs_drift_rhs(C0, pb, pn, W, R, h_m, C_inf)
        k2 = qs_drift_rhs(C0 + 0.5 * dt * k1, pb, pn, W, R, h_m, C_inf)
        k3 = qs_drift_rhs(C0 + 0.5 * dt * k2, pb, pn, W, R, h_m, C_inf)
        k4 = qs_drift_rhs(C0 + dt * k3, pb, pn, W, R, h_m, C_inf)
        C0 += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        tr.append((t, C0))
    return t, tr


# ============================================================================
# 顶层接口
# ============================================================================
def solve_problem3(t_max=345600.0, N=20, out_dt=OUT_DT, verbose=True,
                   T_inf_fn=None, C_inf_fn=None, t_qs_start=36 * 3600.0,
                   kind="uniform", q=1.0):
    """完整求解问题 3。返回 dict（times/T/C 为 60 s 节拍、节点网格上的时空场）。"""
    if T_inf_fn is None or C_inf_fn is None:
        T_inf_fn, C_inf_fn = make_ambient()

    s = DryingSolver(T_inf_fn, C_inf_fn, N=N, out_dt=out_dt, kind=kind, q=q)
    res = s.run(t_max, record=True, stop_cmax=C_TARGET)
    stop = res["stop"]
    if stop is None:
        raise RuntimeError("在 t_max=%.0f s 内未达到烘干判据" % t_max)

    if verbose:
        print("[粗定位] 达标区间 [%.1f s, %.1f s]，宽度 %.1f s"
              % (stop["t_a"], stop["t_b"], stop["t_b"] - stop["t_a"]))

    t_dry, hist = bisect_termination(s, stop, target=C_TARGET, tol=1.0)

    cp_a = dict(t=stop["t_a"], T=stop["T_a"].copy(), C=stop["C_a"].copy())
    T_end, C_end = s.march(cp_a, t_dry, h_max=out_dt)
    w, W = s.w, float(s.w.sum())
    C_inf_late = float(C_inf_fn(T_ENV_END))
    qs = qs_state(t_dry, C_end, w, C_inf=C_inf_late)
    t_qs_frozen, why = qs_extrapolate(qs)

    times_, Cc = res["t"], res["C"]
    pb, pn, C0lo, C0hi, _ = qs_shape_fit(times_, Cc, w, t_qs_start)
    i_qs = min(int(np.searchsorted(times_, t_qs_start)), len(times_) - 1)
    t_qs_drift, _ = qs_extrapolate_drift(Cc[i_qs, 0], times_[i_qs], pb, pn, W,
                                         target=C_TARGET, C_inf=C_inf_late)

    t_min = np.ceil(t_dry / 60.0) * 60.0

    if verbose:
        print("[细定位] t_dry = %.4f s = %.4f h（%.1f min）" % (t_dry, t_dry / 3600, t_dry / 60))
        print("[分钟口径] ceil(t/60) = %.1f min = %.4f h" % (t_min / 60, t_min / 3600))
        print("[准稳态] φ̄=%.4f  φ_N=%.4f  冻结形状速率 a=%.3e /s"
              % (qs["phibar"], qs["phiN"], qs["a"]))
        print("[QS 冻结形状外推] t* = %.1f s（%s）" % (t_qs_frozen, why))
        print("[QS 漂移修正外推] 起始 t=%.1f h → t* = %.4f h（与二分相差 %+.2f min）"
              % (t_qs_start / 3600, t_qs_drift / 3600, (t_qs_drift - t_dry) / 60))
        print("[步数] CN 步总数 %d；宏观步数 %d" % (s.nstep, len(s.log)))

    return dict(times=times_, T=res["T"], C=Cc, r=s.r, grid=s.g,
                t_dry=t_dry, t_dry_min=t_min, bisect_hist=hist,
                T_end=T_end, C_end=C_end, qs=qs, t_qs_extrap=t_qs_frozen,
                qs_reason=why, t_qs_drift=t_qs_drift, qs_fit=(pb, pn, C0lo, C0hi),
                log=s.log, nstep=s.nstep, solver=s, stop=stop)


if __name__ == "__main__":
    import time
    np.set_printoptions(precision=6, suppress=True)
    t0 = time.time()
    r = solve_problem3()
    print("\n总耗时 %.1f s" % (time.time() - t0))
