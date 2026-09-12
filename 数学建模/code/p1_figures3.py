# -*- coding: utf-8 -*-
"""
问题 1 新增配图（§1.5 开头，作为有限差分"差商"概念铺垫）：
  fig20_时空网格剖分.pdf  —— r–t 平面的离散网格，标注 Δr、Δt 与代表节点 u_j^n
  fig21_差商模板.pdf      —— 一阶向前差商 / 二阶中心差商 / 一阶中心差商 三种模板

配色沿用 p1_figures2.py 的设计令牌，与 fig7/fig8 成套。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- 设计令牌（与 p1_figures.py / p1_figures2.py 一致） --------------------
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # 蓝/橙/青/黄
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


def style(ax, title=""):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if title:
        ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=6)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight")
    png = os.path.join(FIG, "_preview_" + os.path.splitext(name)[0] + ".png")
    fig.savefig(png, bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("  ->", path)


def _node(ax, x, y, color, r=0.11, ring=False, label=None, dylab=-0.45, white=True):
    ax.add_patch(Circle((x, y), r, facecolor=color, edgecolor="none", zorder=4))
    if ring:
        ax.add_patch(Circle((x, y), r + 0.06, facecolor="none",
                            edgecolor=CRIT, linewidth=1.5, zorder=5))
    if label is not None:
        ax.text(x, y + dylab, label, ha="center", va="top", color=INK2, fontsize=9, zorder=6)
    return x, y


def _arrow(ax, x0, y0, x1, y1, color, label=None, lx=0, ly=0.16):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                 mutation_scale=11, color=color, lw=1.5, zorder=3))
    if label is not None:
        ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, ha="center",
                va="center", color=color, fontsize=8.5)


# ---------------------------------------------------------------------------
# 图 20：时空网格剖分（r–t 平面）
# ---------------------------------------------------------------------------
def fig_grid():
    NX, NY = 7, 5           # 空间节点 r_0..r_6（r_N=R 表面），时间层 t^0..t^4
    jx, ny = 3, 2           # 代表节点 (r_j, t^n)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    # 网格线
    for i in range(NX):
        ax.plot([i, i], [0, NY - 1], color=GRID, lw=0.9, zorder=1)
    for n in range(NY):
        ax.plot([0, NX - 1], [n, n], color=GRID, lw=0.9, zorder=1)

    # 节点（小灰点）
    for i in range(NX):
        for n in range(NY):
            ax.plot(i, n, "o", ms=3.0, color=MUTED, zorder=3)

    # 代表节点（蓝点 + 红环）
    ax.add_patch(Circle((jx, ny), 0.15, facecolor=CAT[0], edgecolor="none", zorder=6))
    ax.add_patch(Circle((jx, ny), 0.21, facecolor="none", edgecolor=CRIT,
                        linewidth=1.5, zorder=7))
    ax.text(jx, ny - 0.48, "$u_j^{\\,n}$", ha="center", va="top", color=INK, fontsize=11, zorder=6)

    # 空间刻度 r_j（r_0=轴心，r_N=表面）
    ax.set_xticks(range(NX))
    ax.set_xticklabels(["$r_0$", "$r_1$", "$\\cdots$", "$r_j$", "$\\cdots$",
                        "$r_{N-1}$", "$r_N$"], fontsize=9)
    # 时间刻度 t^n
    ax.set_yticks(range(NY))
    ax.set_yticklabels(["$t^0$", "$t^1$", "$t^n$", "$\\cdots$", "$t^M$"], fontsize=9)

    # Δr 标注（r_0 与 r_1 之间，下方）
    ax.annotate("", xy=(1, -0.28), xytext=(0, -0.28),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2))
    ax.text(0.5, -0.56, "$\\Delta r$", ha="center", va="top", color=INK2, fontsize=10)

    # Δt 标注（t^0 与 t^1 之间，左侧）
    ax.annotate("", xy=(-0.32, 1), xytext=(-0.32, 0),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2))
    ax.text(-0.72, 0.5, "$\\Delta t$", ha="center", va="center", color=INK2, fontsize=10)

    # 轴心 / 表面标注
    ax.text(0, -1.0, "轴心 $r=0$", ha="center", va="top", color=MUTED, fontsize=8)
    ax.text(NX - 1, -1.0, "表面 $r=R$", ha="center", va="top", color=MUTED, fontsize=8)

    ax.set_xlim(-1.05, NX - 0.4)
    ax.set_ylim(-1.35, NY - 0.25)
    style(ax, "时空网格剖分（径向 $r$ × 时间 $t$，节点 $u_j^{\\,n}=u(r_j,t^n)$）")
    ax.set_xlabel("空间坐标 $r$（径向）", fontsize=9, color=INK2, labelpad=4)
    ax.set_ylabel("时间 $t$", fontsize=9, color=INK2, labelpad=4)

    save(fig, "fig20_时空网格剖分.pdf")


# ---------------------------------------------------------------------------
# 图 21：差商模板（向前差商 / 二阶中心差商 / 一阶中心差商）
# ---------------------------------------------------------------------------
def fig_difference():
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.3))
    for ax in axes:
        style(ax)
        ax.set_xlim(-0.4, 2.6)
        ax.set_ylim(-1.6, 1.7)

    # ---- (a) 一阶向前差商（时间方向） ----
    ax = axes[0]
    ax.set_title("(a) 一阶向前差商（时间 $\\partial u/\\partial t$）")
    _node(ax, 0.6, 1.1, CAT[0], label=None, dylab=0.18)
    ax.text(0.6, 1.32, "$u_j^{n+1}$", ha="center", va="bottom", color=INK, fontsize=10)
    _node(ax, 0.6, -0.7, CAT[0], label=None, dylab=0.18)
    ax.text(0.6, -0.48, "$u_j^{n}$", ha="center", va="bottom", color=INK, fontsize=10)
    _arrow(ax, 0.6, -0.55, 0.6, 0.95, MUTED)
    ax.text(1.0, 0.2, "$\\Delta t$", ha="left", va="center", color=INK2, fontsize=9)
    ax.text(1.1, -1.42,
            "$\\dfrac{\\partial u}{\\partial t}\\approx\\dfrac{u_j^{n+1}-u_j^{n}}{\\Delta t}$",
            ha="center", color=INK, fontsize=11)

    # ---- (b) 二阶中心差商（空间方向） ----
    ax = axes[1]
    ax.set_title("(b) 二阶中心差商（空间 $\\partial^2 u/\\partial r^2$）")
    _node(ax, 0.35, 0.35, CAT[1], label=None, dylab=0.0)
    ax.text(0.35, 0.0, "$u_{j-1}$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(0.35, 0.78, "$+1$", ha="center", color=CRIT, fontsize=9)
    _node(ax, 1.3, 0.35, CAT[1], ring=True, label=None, dylab=0.0)
    ax.text(1.3, 0.0, "$u_j$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(1.3, 0.78, "$-2$", ha="center", color=CRIT, fontsize=9)
    _node(ax, 2.25, 0.35, CAT[1], label=None, dylab=0.0)
    ax.text(2.25, 0.0, "$u_{j+1}$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(2.25, 0.78, "$+1$", ha="center", color=CRIT, fontsize=9)
    _arrow(ax, 0.5, 0.35, 1.15, 0.35, MUTED)
    _arrow(ax, 1.45, 0.35, 2.1, 0.35, MUTED)
    ax.text(0.825, 0.02, "$\\Delta r$", ha="center", color=INK2, fontsize=8.5)
    ax.text(1.775, 0.02, "$\\Delta r$", ha="center", color=INK2, fontsize=8.5)
    ax.text(1.3, -1.42,
            "$\\dfrac{\\partial^2 u}{\\partial r^2}\\approx\\dfrac{u_{j-1}-2u_j+u_{j+1}}{\\Delta r^2}$",
            ha="center", color=INK, fontsize=11)

    # ---- (c) 一阶中心差商（空间方向） ----
    ax = axes[2]
    ax.set_title("(c) 一阶中心差商（空间 $\\partial u/\\partial r$）")
    _node(ax, 0.35, 0.35, CAT[2], label=None, dylab=0.0)
    ax.text(0.35, 0.0, "$u_{j-1}$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(0.35, 0.78, "$-1$", ha="center", color=CRIT, fontsize=9)
    _node(ax, 1.3, 0.35, CAT[2], ring=True, label=None, dylab=0.0)
    ax.text(1.3, 0.0, "$u_j$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(1.3, 0.78, "$0$", ha="center", color=CRIT, fontsize=9)
    _node(ax, 2.25, 0.35, CAT[2], label=None, dylab=0.0)
    ax.text(2.25, 0.0, "$u_{j+1}$", ha="center", va="top", color=INK, fontsize=10)
    ax.text(2.25, 0.78, "$+1$", ha="center", color=CRIT, fontsize=9)
    _arrow(ax, 0.5, 0.35, 1.15, 0.35, MUTED)
    _arrow(ax, 1.45, 0.35, 2.1, 0.35, MUTED)
    ax.text(0.825, 0.02, "$\\Delta r$", ha="center", color=INK2, fontsize=8.5)
    ax.text(1.775, 0.02, "$\\Delta r$", ha="center", color=INK2, fontsize=8.5)
    ax.text(1.3, -1.42,
            "$\\dfrac{\\partial u}{\\partial r}\\approx\\dfrac{u_{j+1}-u_{j-1}}{2\\Delta r}$",
            ha="center", color=INK, fontsize=11)

    fig.subplots_adjust(wspace=0.42)
    save(fig, "fig21_差商模板.pdf")


if __name__ == "__main__":
    print("生成差商概念配图：")
    fig_grid()
    fig_difference()
    print("完成")
