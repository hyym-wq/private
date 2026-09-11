# -*- coding: utf-8 -*-
"""
问题 2 配图：
  fig9_变参数物性.pdf     —— 附录 3 变参数 ρ, cp, k, D 随含水率/温度的变化
  fig10_全烘干过程剖面.pdf —— 温度（预热期）与水分浓度（全烘干）的径向剖面
  fig11_收敛性.pdf        —— Taylor(Newton) 线性化 vs 0 阶滞后的迭代收敛

配色沿用 p1_figures.py 设计令牌。
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_model as m

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
RES = os.path.join(BASE, "results")
os.makedirs(FIG, exist_ok=True)

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SEQ7 = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#184f95"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CRIT = "#d03b3b"

plt.rcParams.update({
    "font.family": ["Microsoft YaHei", "SimHei"],
    "font.size": 9,
    "axes.unicode_minus": False,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 2.0,
    "figure.dpi": 120,
})


def style(ax, xlabel="", ylabel="", title=""):
    ax.grid(True, color=GRID, linewidth=0.6, alpha=1.0, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, color=INK, fontsize=10, loc="left", pad=8)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("  ->", path)


def load_or_run(n_store=1800):
    """全 3 天模拟（每 0.5 h 存一帧），优先读取已存 npz。"""
    npz = os.path.join(RES, "p2_fields.npz")
    if os.path.exists(npz):
        d = np.load(npz)
        return d["times"], d["r"], d["T"], d["C"]
    T_inf, C_inf = m.make_ambient()
    times, r, T, C = m.solve(T_inf, C_inf, t_end=259200.0, dt=1.0, N=20,
                             n_store=n_store, picard_max=3, newton=True)
    np.savez(npz, times=times, r=r, T=T, C=C)
    return times, r, T, C


# ---------------------------------------------------------------------------
# 图 9：变参数物性
# ---------------------------------------------------------------------------
def fig_param():
    Cc = np.linspace(0.05, 2.55, 400)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))

    ax = axes[0]
    ax.plot(Cc, m.rho(Cc) / m.rho(2.55), color=CAT[0], label=r"$\rho(C)/\rho_0$")
    ax.plot(Cc, m.cp(Cc) / m.cp(2.55), color=CAT[1], label=r"$c_p(C)/c_{p0}$")
    ax.plot(Cc, m.kcond(Cc) / m.kcond(2.55), color=CAT[2], label=r"$k(C)/k_0$")
    ax.axhline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    style(ax, "含水率 $C$ / (kg·kg$^{-1}$)", "相对初值的归一化取值",
          "(a) 密度/比热/导热系数随含水率变化")
    ax.legend(loc="lower left")
    ax.set_xlim(0, 2.55); ax.set_ylim(0.2, 1.7)

    ax = axes[1]
    for Tlab, col, ls in ((28.0, CAT[0], "-"), (50.0, CAT[1], "--")):
        ax.semilogy(Cc, m.diff_coef(Cc, Tlab), color=col, linestyle=ls,
                    label=f"固定温度 $T$ = {Tlab:.0f} °C")
    style(ax, "含水率 $C$ / (kg·kg$^{-1}$)", "扩散系数 $D(C,T)$ / (m$^2$·s$^{-1}$，对数)",
          "(b) 扩散系数随含水率与温度变化")
    ax.legend(loc="upper left")
    ax.set_xlim(0, 2.55)

    fig.subplots_adjust(wspace=0.30)
    save(fig, "fig9_变参数物性.pdf")


# ---------------------------------------------------------------------------
# 图 10：全烘干过程剖面
# ---------------------------------------------------------------------------
def fig_dry(times, r, T, C):
    r_cm = r * 100.0
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))

    # (a) 温度：预热期
    for k, th in enumerate([0.0, 0.5, 1.0, 2.0, 3.0]):
        i = int(np.argmin(np.abs(times - th * 3600.0)))
        axes[0].plot(r_cm, T[i], color=SEQ7[k + 1], lw=2.0, zorder=3,
                     linestyle="-" if k % 2 == 0 else "--", label=f"t = {times[i]/3600:g} h")
    style(axes[0], "到药材中心的距离 $r$ / cm", "温度 $T$ / °C", "(a) 预热平衡阶段温度剖面")
    axes[0].set_xlim(0, 2.1); axes[0].set_ylim(26, 52)
    axes[0].xaxis.set_major_locator(MultipleLocator(0.5))
    axes[0].legend(loc="center left", ncol=2, handlelength=1.6)

    # (b) 水分：全烘干
    for k, th in enumerate([0.0, 6.0, 12.0, 24.0, 48.0, 72.0]):
        i = int(np.argmin(np.abs(times - th * 3600.0)))
        axes[1].plot(r_cm, C[i], color=SEQ7[k], lw=2.0, zorder=3,
                     linestyle="-" if k % 2 == 0 else "--", label=f"t = {times[i]/3600:g} h")
    axes[1].axhline(0.15, color=CRIT, lw=1.2, ls=(0, (4, 3)), zorder=2)
    axes[1].text(0.07, 0.18, "烘干标准 $C$=0.15", color=CRIT, fontsize=8,
                 va="bottom", ha="left")
    style(axes[1], "到药材中心的距离 $r$ / cm", "水分浓度 $C$ / (kg·kg$^{-1}$)",
          "(b) 全烘干过程水分浓度剖面")
    axes[1].set_xlim(0, 2.1); axes[1].set_ylim(0, 2.7)
    axes[1].xaxis.set_major_locator(MultipleLocator(0.5))
    axes[1].legend(loc="upper right", ncol=2, handlelength=1.6)

    fig.subplots_adjust(wspace=0.30)
    save(fig, "fig10_全烘干过程剖面.pdf")


# ---------------------------------------------------------------------------
# 图 11：收敛性（Taylor/Newton vs 0 阶滞后）
# ---------------------------------------------------------------------------
def _one_step_iterate_change(newton, T_old, C_old, Ti_o, Ti_n, Ci_o, Ci_n, dt, dr, N, theta):
    C_star, T_star = C_old.copy(), T_old.copy()
    zero = np.zeros(N + 1)
    ah_o, bh_o, ch_o, gT_o = m.build_op_heat(dr, N, C_old)
    sT_o = zero.copy(); sT_o[N] = gT_o * Ti_o
    res = []
    for _ in range(8):
        ah, bh, ch, gT = m.build_op_heat(dr, N, C_star)
        sT_n = zero.copy(); sT_n[N] = gT * Ti_n
        T_new = m._step(T_old, (ah, bh, ch), (ah_o, bh_o, ch_o), sT_o, sT_n, dt, theta)
        C_new = m._mass_step(C_old, C_star, T_old, T_new, Ci_o, Ci_n, dt, theta, dr, N, m.R_CYL, newton)
        res.append(float(np.abs(C_new - C_star).max()))
        C_star, T_star = C_new, T_new
    return res


def fig_conv():
    T_inf, C_inf = m.make_ambient()
    N, dt, theta = 20, 1.0, 0.5
    dr = m.R_CYL / N
    # 取预热期 t=1799 -> 1800 s 的一个 dt=1 s 步作局部示例。
    times, r, T, C = m.solve(T_inf, C_inf, t_end=1799.0, dt=dt, N=N,
                             n_store=None, picard_max=3, newton=True)
    T_old, C_old = T[-1].copy(), C[-1].copy()
    Ti_o, Ci_o = T_inf(1799.0), C_inf(1799.0)
    Ti_n, Ci_n = T_inf(1800.0), C_inf(1800.0)

    resN = _one_step_iterate_change(True, T_old, C_old, Ti_o, Ti_n, Ci_o, Ci_n, dt, dr, N, theta)
    resL = _one_step_iterate_change(False, T_old, C_old, Ti_o, Ti_n, Ci_o, Ci_n, dt, dr, N, theta)

    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    k = np.arange(1, len(resN) + 1)
    ax.semilogy(k, np.array(resN) + 1e-16, color=CAT[0], lw=2.0, marker="o", ms=4,
                label="Taylor(Newton) 线性化", zorder=3)
    ax.semilogy(k, np.array(resL) + 1e-16, color=CAT[1], lw=2.0, marker="s", ms=4,
                label="0 阶滞后 Picard", zorder=3)
    ax.axhline(1e-9, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(1.02, 1e-9, "浓度迭代差阈值 $10^{-9}$", color=MUTED, fontsize=8, va="bottom")
    style(ax, "非线性迭代次数 $k$", r"$\Vert C^{(k)}-C^{(k-1)}\Vert_\infty$ / (kg/kg)",
          "单步迭代差：$t=1799$ 至 $1800$ s，$\\Delta t=1$ s")
    ax.text(0.04, 0.05, "浓度相邻迭代差；单步示例", transform=ax.transAxes,
            fontsize=7.5, color=INK2)
    leg = ax.legend(loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_xlim(0.9, 8.1)
    ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])
    save(fig, "fig11_收敛性.pdf")


if __name__ == "__main__":
    print("生成问题 2 配图：")
    fig_param()
    times, r, T, C = load_or_run()
    fig_dry(times, r, T, C)
    fig_conv()
    print("完成")
