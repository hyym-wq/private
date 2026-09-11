# -*- coding: utf-8 -*-
"""
问题 1 新增配图（配合"先 Euler 显式、后 CN"的叙事）：
  fig7_节点分类离散.pdf   —— 中心点 / 内部点 / 边界点三类节点的离散示意图
  fig8_稳定性对比.pdf     —— Euler 显式（条件稳定） vs Crank-Nicolson（无条件稳定）

配色沿用 p1_figures.py 的设计令牌。
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_model import ALPHA, R_CYL, T0, C0
from p1_euler import solve_explicit
from p1_model import solve as solve_cn

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- 设计令牌（与 p1_figures.py 一致） ------------------------------------
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # 分类：蓝/橙/青/黄
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
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=6)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight")
    png = os.path.join(FIG, "_preview_" + os.path.splitext(name)[0] + ".png")
    fig.savefig(png, bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("  ->", path)


def _node(ax, x, y, color, r=0.10, ring=False, label=None, dylab=-0.42):
    ax.add_patch(Circle((x, y), r, facecolor=color, edgecolor="none", zorder=4))
    if ring:
        ax.add_patch(Circle((x, y), r + 0.07, facecolor="none",
                            edgecolor=CRIT, linewidth=1.6, zorder=5))
    if label is not None:
        ax.text(x, y + dylab, label, ha="center", va="top", color=INK2, fontsize=9)


def _flux(ax, x0, x1, y, color, both=False, label=None, dylab=0.0):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                 mutation_scale=11, color=color, lw=1.5, zorder=3))
    if both:
        ax.add_patch(FancyArrowPatch((x1, y), (x0, y), arrowstyle="-|>",
                     mutation_scale=11, color=color, lw=1.5, zorder=3))
    if label is not None:
        ax.text((x0 + x1) / 2, y + 0.22 + dylab, label, ha="center", va="bottom",
                color=color, fontsize=8.5)


# ---------------------------------------------------------------------------
# 图 7：三类节点离散示意
# ---------------------------------------------------------------------------
def fig_nodes():
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.4))
    for ax in axes:
        style(ax)
        ax.set_xlim(-0.3, 3.2)
        ax.set_ylim(-1.35, 1.55)

    y = 0.55

    # ---- (a) 中心点 j=0 ----
    ax = axes[0]
    ax.set_title("(a) 中心点 $j=0$（L'H\\^opital + 对称）")
    _node(ax, 0.4, y, CAT[0], ring=True, label="$T_{-1}=T_1$（虚拟）")
    _node(ax, 1.4, y, CAT[0], label="$T_1$")
    # 中心节点（目标）
    ax.add_patch(Circle((0.4, y), 0.10, facecolor=CAT[0], edgecolor="none", zorder=4))
    ax.text(0.4, y, "$T_0$", ha="center", va="center", color="white", fontsize=8, zorder=6)
    _flux(ax, 0.52, 1.28, y, MUTED, label="$\\partial_r u|_0=0$", dylab=0.02)
    ax.text(0.4, -0.18, "$\\frac{\\mathrm{d}T_0}{\\mathrm{d}t}=\\frac{4\\alpha}{\\Delta r^2}(T_1-T_0)$",
            ha="center", color=INK, fontsize=10)
    ax.text(0.4, -0.72, "对称面（无通量）", ha="center", color=MUTED, fontsize=8)

    # ---- (b) 内部点 j ----
    ax = axes[1]
    ax.set_title("(b) 内部点 $1\\leq j\\leq N-1$（守恒中心差分）")
    _node(ax, 0.5, y, CAT[1], label="$T_{j-1}$")
    _node(ax, 1.5, y, CAT[1], label="$T_{j+1}$")
    ax.add_patch(Circle((1.0, y), 0.10, facecolor=CAT[1], edgecolor="none", zorder=4))
    ax.add_patch(Circle((1.0, y), 0.17, facecolor="none", edgecolor=CRIT, linewidth=1.6, zorder=5))
    ax.text(1.0, y, "$T_j$", ha="center", va="center", color="white", fontsize=8, zorder=6)
    _flux(ax, 0.62, 0.9, y, MUTED, both=True, label="$r_{j-1/2}$")
    _flux(ax, 1.1, 1.38, y, MUTED, both=True, label="$r_{j+1/2}$")
    ax.text(1.0, -0.18, "$\\frac{\\mathrm{d}T_j}{\\mathrm{d}t}=\\frac{\\alpha}{\\Delta r^2}\\left[\\frac{j-1/2}{j}(T_{j-1}-T_j)+\\frac{j+1/2}{j}(T_{j+1}-T_j)\\right]$",
            ha="center", color=INK, fontsize=8.6)
    ax.text(1.0, -0.72, "两侧净导入的导热热流", ha="center", color=MUTED, fontsize=8)

    # ---- (c) 边界点 j=N ----
    ax = axes[2]
    ax.set_title("(c) 边界点 $j=N$（Robin 半控制体）")
    _node(ax, 0.5, y, CAT[2], label="$T_{N-1}$")
    ax.add_patch(Circle((1.5, y), 0.10, facecolor=CAT[2], edgecolor="none", zorder=4))
    ax.add_patch(Circle((1.5, y), 0.17, facecolor="none", edgecolor=CRIT, linewidth=1.6, zorder=5))
    ax.text(1.5, y, "$T_N$", ha="center", va="center", color="white", fontsize=8, zorder=6)
    _flux(ax, 0.62, 1.38, y, MUTED, both=True, label="$r_{N-1/2}$")
    _flux(ax, 1.62, 2.5, y, CAT[3], label="$h(T_N-T_\\infty)$", dylab=0.02)
    ax.text(2.72, y, "$T_\\infty$", ha="center", va="center", color=CAT[3], fontsize=9)
    ax.text(1.5, -0.18, "$\\frac{\\mathrm{d}T_N}{\\mathrm{d}t}=\\frac{\\alpha r_{N-1/2}}{V_N\\Delta r}(T_{N-1}-T_N)-\\frac{Rh}{\\rho c_p V_N}(T_N-T_\\infty)$",
            ha="center", color=INK, fontsize=8.6)
    ax.text(1.5, -0.72, "内部导热 $+$ 表面对流", ha="center", color=MUTED, fontsize=8)

    fig.subplots_adjust(wspace=0.34)
    save(fig, "fig7_节点分类离散.pdf")


# ---------------------------------------------------------------------------
# 图 8：稳定性对比
# ---------------------------------------------------------------------------
def fig_stability():
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.5))
    Tf = lambda t: 50.0
    Cf = lambda t: 0.02
    N = 20
    dr = R_CYL / N

    # ---- (a) Euler 显式：log 刻度显示 CFL 临界 ----
    ax = axes[0]
    dts = np.array([0.5, 1.0, 1.5, 2.0, 2.3, 2.4, 2.44, 2.45, 2.46, 2.5, 2.6, 2.8, 3.0, 4.0])
    tmax = []
    for dt in dts:
        _, _, T, C = solve_explicit(Tf, Cf, t_end=1800.0, dt=dt, N=N)
        tmax.append(float(np.abs(T).max()))
    tmax = np.array(tmax)
    fo = ALPHA * dts / dr ** 2
    # 临界线（特征值导出）
    fo_crit = 0.413
    dt_crit = fo_crit * dr ** 2 / ALPHA
    ax.axvline(dt_crit, color=CRIT, lw=1.4, ls=(0, (4, 3)), zorder=2)
    ax.annotate(f"临界 $\\Delta t_{{\\mathrm{{crit}}}}$={dt_crit:.2f} s\n($Fo$={fo_crit:.2f})",
                xy=(dt_crit, 1e4), xytext=(dt_crit + 0.12, 2e6),
                color=CRIT, fontsize=8, va="bottom")
    ax.semilogy(dts, tmax, color=CAT[0], lw=2.0, marker="o", ms=4, zorder=3)
    ax.axhline(50.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(3.2, 50.0, "物理上界 $T$=50 °C", color=MUTED, fontsize=8, va="bottom", ha="right")
    style(ax, "时间步长 $\\Delta t$ / s", "最大温度 $T_{\\max}$ / °C（对数刻度）",
          "(a) Euler 显式：$Fo>0.41$ 即发散")
    ax.set_xlim(0, 4.2)
    ax.set_ylim(20, 1e200)

    # ---- (b) Crank-Nicolson：无条件稳定 ----
    ax = axes[1]
    dts_cn = np.array([1.0, 10.0, 60.0, 120.0, 240.0, 300.0])
    tmax_cn = []
    for dt in dts_cn:
        _, _, T, C = solve_cn(Tf, Cf, t_end=1800.0, dt=dt, N=N, theta=0.5,
                              n_store=None, bc="halfcell", picard_max=3)
        tmax_cn.append(float(T.max()))
    ax.plot(dts_cn, tmax_cn, color=CAT[2], lw=2.0, marker="o", ms=4, zorder=3)
    ax.axhline(50.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.text(30, 50.2, "物理上界 $T$=50 °C", color=MUTED, fontsize=8, va="bottom", ha="center")
    style(ax, "时间步长 $\\Delta t$ / s", "最大温度 $T_{\\max}$ / °C（线性刻度）",
          "(b) Crank--Nicolson：任意 $\\Delta t$ 有界")
    ax.set_xlim(0, 320)
    ax.set_ylim(44, 51)

    fig.subplots_adjust(wspace=0.30)
    save(fig, "fig8_稳定性对比.pdf")


if __name__ == "__main__":
    print("生成新增配图：")
    fig_nodes()
    fig_stability()
    print("完成")
