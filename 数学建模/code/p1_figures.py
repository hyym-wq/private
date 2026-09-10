# -*- coding: utf-8 -*-
"""
问题 1 配图生成：输出 PDF 到 figures/

配色遵循已验证的数据可视化规范：
  * 有序序列（时刻系列）用单一蓝色阶的 4 个步长
        #86b6ef → #5598e7 → #2a78d6 → #184f95
    该 4 步组合通过顺序色阶的全部校验（单色相、亮度单调、步距 >= 0.06、浅端对比 >= 2:1）
  * 分类序列（中心/中间/表面、数值/解析、T/C）用前 3 个分类槽位
        #2a78d6 (蓝) · #eb6834 (橙) · #1baf7a (青)
    通过 CVD 分离度与常视觉分离度校验
  * 坐标轴/网格使用中性墨色，文字一律使用墨色而非系列色
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1_model import solve, R_CYL, conservative_weights
from p1_verify import analytic_cylinder
from p1_run import load_ambient, make_linear

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
RES = os.path.join(BASE, "results")
os.makedirs(FIG, exist_ok=True)

# ---- 设计令牌 --------------------------------------------------------------
SEQ4 = ["#86b6ef", "#5598e7", "#2a78d6", "#184f95"]      # 有序：时刻系列
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]        # 分类：实体系列
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
    png = os.path.join(FIG, "_preview_" + os.path.splitext(name)[0] + ".png")
    fig.savefig(png, bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("  ->", path)


# ---------------------------------------------------------------------------
# 载入结果
# ---------------------------------------------------------------------------
d = np.load(os.path.join(RES, "p1_fields.npz"))
times, r, T, C = d["times"], d["r"], d["T"], d["C"]
r_cm = r * 100.0

T_TIMES = [100.0, 300.0, 900.0, 1800.0]


def idx(t):
    return int(round(t / 1.0))


# ---------------------------------------------------------------------------
# 图 1/2：温度与水分浓度剖面
# ---------------------------------------------------------------------------
def fig_profiles():
    panels = [("温度", "温度 $T$ / °C", T, "#eef4fb"),
              ("水分浓度", "水分浓度 $C$ / (kg·kg$^{-1}$)", C, "#eef4fb")]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))

    for ax, (tag, ylab, field, _) in zip(axes, panels):
        for k, t in enumerate(T_TIMES):
            ax.plot(r_cm, field[idx(t)], color=SEQ4[k], lw=2.0,
                    zorder=3, solid_capstyle="round")
            ax.annotate(f"{int(t)} s", xy=(r_cm[-1], field[idx(t)][-1]),
                        xytext=(4, 0), textcoords="offset points",
                        color=INK2, fontsize=8, va="center")
        style(ax, "到药材中心的距离 $r$ / cm", ylab)
        ax.set_xlim(0, 2.28)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
    axes[0].set_title("(a) 温度剖面", color=INK, fontsize=10, loc="left", pad=8)
    axes[1].set_title("(b) 水分浓度剖面", color=INK, fontsize=10, loc="left", pad=8)
    fig.subplots_adjust(wspace=0.32)
    save(fig, "fig1_预热平衡阶段剖面.pdf")


# ---------------------------------------------------------------------------
# 图 3：特征位置时程曲线（上下两面板，避免双纵轴）
# ---------------------------------------------------------------------------
def fig_history():
    pos = [(0, "中心 $r$=0"), (10, "中间 $r$=1.0 cm"), (20, "表面 $r$=2.0 cm")]
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 5.0), sharex=True)

    for k, (j, lab) in enumerate(pos):
        axes[0].plot(times, T[:, j], color=CAT[k], lw=2.0, label=lab, zorder=3)
        axes[1].plot(times, C[:, j], color=CAT[k], lw=2.0, label=lab, zorder=3)

    style(axes[0], "", "温度 $T$ / °C", "(a) 特征位置的温度时程")
    style(axes[1], "时间 $t$ / s", "水分浓度 $C$ / (kg·kg$^{-1}$)", "(b) 特征位置的水分浓度时程")

    # 烘房环境（附件 1 实测数据；中性虚线 = 参照，不作为系列）
    ATT = os.path.join(os.path.dirname(BASE), "A题", "附件", "附件1.xlsx")
    t_amb, T_amb, C_amb = load_ambient(ATT)
    m = t_amb <= 1800
    axes[0].plot(t_amb[m], T_amb[m], color=MUTED, lw=1.4, ls=(0, (4, 3)), zorder=2)
    axes[0].annotate("烘房环境温度", xy=(t_amb[m][-1], T_amb[m][-1]), xytext=(-6, -14),
                     textcoords="offset points", ha="right", color=MUTED, fontsize=8)
    scl = 100.0
    axes[1].plot(t_amb[m], C_amb[m] * scl, color=MUTED, lw=1.4, ls=(0, (4, 3)), zorder=2)
    axes[1].annotate(f"烘房水汽浓度（×{int(scl)}）", xy=(t_amb[m][-1], C_amb[m][-1] * scl),
                     xytext=(-6, -13), textcoords="offset points", ha="right",
                     color=MUTED, fontsize=8)

    axes[0].set_ylim(27.2, 44.0)
    axes[1].set_ylim(1.30, 3.00)
    for ax, loc in zip(axes, ("upper left", "lower left")):
        ax.set_xlim(0, 1800)
        leg = ax.legend(loc=loc)
        for txt in leg.get_texts():
            txt.set_color(INK2)
    fig.subplots_adjust(hspace=0.28)
    save(fig, "fig3_特征位置时程.pdf")


# ---------------------------------------------------------------------------
# 图 4：与解析解对比（验证 Robin 边界与中心处理）
# ---------------------------------------------------------------------------
def fig_analytic():
    T_inf = 50.0
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.4))

    for k, t in enumerate([300.0, 1800.0]):
        _, rr, Th, Ch = solve(lambda s: T_inf, lambda s: 0.0, t_end=t, dt=0.1,
                              N=20, theta=0.5, n_store=None, bc="halfcell",
                              picard_max=3)
        num = Th[-1]
        ex = analytic_cylinder(rr, t, T_inf)
        axes[0].plot(rr * 100, ex, color=SEQ4[3 - 2 * k], lw=3.4, alpha=0.35,
                     zorder=2)
        axes[0].plot(rr * 100, num, color=CAT[k], lw=1.6, ls=(0, (5, 2)),
                     marker="o", ms=3.2, markevery=2, zorder=3,
                     label=f"t={int(t)} s")
        axes[1].semilogy(rr * 100, np.abs(num - ex), color=CAT[k], lw=2.0,
                         marker="o", ms=3.4, markevery=2, zorder=3,
                         label=f"t={int(t)} s")

    # 左图：实线=解析，虚线+圆点=数值
    axes[0].plot([], [], color=MUTED, lw=1.6, ls=(0, (5, 2)), marker="o", ms=3.2,
                 label="Crank–Nicolson 数值解")
    axes[0].plot([], [], color=MUTED, lw=3.4, alpha=0.35, label="级数解析解")
    style(axes[0], "到药材中心的距离 $r$ / cm", "温度 $T$ / °C",
          "(a) 数值解与解析解剖面（$T_\\infty$=50 °C）")
    axes[0].legend(loc="lower right")
    axes[0].set_xlim(0, 2)

    style(axes[1], "到药材中心的距离 $r$ / cm", "数值解与解析解的温度偏差 / °C",
          "(b) 逐点绝对偏差")
    axes[1].legend(loc="upper left")
    axes[1].set_xlim(0, 2)
    axes[1].yaxis.grid(True, which="minor", color=GRID, linewidth=0.4)
    for ax in axes:
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
    fig.subplots_adjust(wspace=0.30)
    save(fig, "fig4_解析解验证.pdf")


# ---------------------------------------------------------------------------
# 图 5：收敛阶
# ---------------------------------------------------------------------------
def fig_convergence():
    Tf = lambda t: 28.0 + 22.0 * (1 - np.exp(-t / 500.0))
    Cf = lambda t: 2.55 - 2.53 * (1 - np.exp(-t / 500.0))
    t_end = 600.0

    # 空间
    _, _, Tr_h, Cr_h = solve(Tf, Cf, t_end=t_end, dt=0.02, N=640, theta=0.5,
                             n_store=None, bc="halfcell", picard_max=3)
    Tr, Cr = Tr_h[-1], Cr_h[-1]
    Ns = np.array([20, 40, 80, 160])
    eT_s, eC_s = [], []
    for N in Ns:
        _, _, Th, Ch = solve(Tf, Cf, t_end=t_end, dt=0.02, N=N, theta=0.5,
                             n_store=None, bc="halfcell", picard_max=3)
        s = 640 // N
        eT_s.append(np.abs(Th[-1] - Tr[::s]).max())
        eC_s.append(np.abs(Ch[-1] - Cr[::s]).max())

    # 时间
    _, _, Tr_h, Cr_h = solve(Tf, Cf, t_end=t_end, dt=0.02, N=40, theta=0.5,
                             n_store=None, bc="halfcell", picard_max=3)
    Tr2, Cr2 = Tr_h[-1], Cr_h[-1]
    dts = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    eT_t, eC_t = [], []
    for dt in dts:
        _, _, Th, Ch = solve(Tf, Cf, t_end=t_end, dt=dt, N=40, theta=0.5,
                             n_store=None, bc="halfcell", picard_max=3)
        eT_t.append(np.abs(Th[-1] - Tr2).max())
        eC_t.append(np.abs(Ch[-1] - Cr2).max())

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))

    def plot_panel(ax, x, ys, labels, xlabel, title, ticks, legend_loc, anchor=0):
        # 参考斜率线：斜率恒为 2，整体下移 0.25 倍使其与数据分离，仅作视觉参照
        xa, ya = x[anchor], ys[0][anchor]
        ref = 0.25 * ya * (x / xa) ** 2
        ax.plot(x, ref, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=2,
                label="参考斜率 2")
        for k, (y, lab) in enumerate(zip(ys, labels)):
            ax.loglog(x, y, color=CAT[k], lw=2.0, marker="o", ms=5,
                      markeredgecolor=SURFACE, markeredgewidth=1.4,
                      label=lab, zorder=3)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v:g}" for v in ticks])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        style(ax, xlabel, "与参考解的最大偏差（对数刻度）", title)
        leg = ax.legend(loc=legend_loc)
        for t in leg.get_texts():
            t.set_color(INK2)

    plot_panel(axes[0], Ns, [eT_s, eC_s], ["温度 $T$", "水分浓度 $C$"],
               "径向网格数 $N$（$\\Delta r=R/N$）", "(a) 空间收敛：固定 $\\Delta t$=0.02 s",
               [20, 40, 80, 160], "upper right", anchor=-1)
    plot_panel(axes[1], dts, [eT_t, eC_t], ["温度 $T$", "水分浓度 $C$"],
               "时间步长 $\\Delta t$ / s", "(b) 时间收敛：固定 $N$=40",
               [0.25, 0.5, 1, 2, 4], "lower right", anchor=0)

    fig.subplots_adjust(wspace=0.30)
    save(fig, "fig5_收敛阶验证.pdf")


# ---------------------------------------------------------------------------
# 图 6：格式对比
# ---------------------------------------------------------------------------
def fig_format():
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.4))

    # (a) 守恒偏差
    labels = ["虚拟节点法\n(ghost)", "半控制体积法\n(halfcell)"]
    vals = [3.061e-2, 1.788e-14]
    bars = axes[0].bar(labels, vals, width=0.42,
                       color=[CAT[1], CAT[0]], zorder=3)
    for b, v in zip(bars, vals):
        axes[0].annotate(f"{v:.2e}", xy=(b.get_x() + b.get_width() / 2, v),
                         xytext=(0, 5), textcoords="offset points",
                         ha="center", color=INK, fontsize=8.5)
    axes[0].set_yscale("log")
    axes[0].set_ylim(1e-15, 1)
    style(axes[0], "", "质量守恒相对偏差 $|\\Delta M+M_{flux}|/|M_{flux}|$",
          "(a) 表面 Robin 离散的守恒性")
    axes[0].tick_params(axis="x", length=0)

    # (b) 收敛阶对比
    groups = ["单次滞后\n$k_{iter}$=1", "Picard 迭代\n$k_{iter}\\geq$2"]
    x = np.arange(2)
    w = 0.34
    order_C = [1.00, 2.00]
    axes[1].bar(x - w / 2, [2.00, 2.00], width=w, color=CAT[0],
                label="温度 $T$", zorder=3)
    axes[1].bar(x + w / 2, order_C, width=w, color=CAT[1],
                label="水分浓度 $C$", zorder=3)
    for xi, v in zip(x - w / 2, [2.00, 2.00]):
        axes[1].annotate(f"{v:.2f}", xy=(xi, v), xytext=(0, 4),
                         textcoords="offset points", ha="center", color=INK, fontsize=8.5)
    for xi, v in zip(x + w / 2, order_C):
        axes[1].annotate(f"{v:.2f}", xy=(xi, v), xytext=(0, 4),
                         textcoords="offset points", ha="center", color=INK, fontsize=8.5)
    axes[1].axhline(2.0, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=2)
    axes[1].annotate("虚线 = 二阶精度", xy=(1.48, 2.30), ha="right", va="center",
                     color=MUTED, fontsize=8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(groups)
    axes[1].set_ylim(0, 2.6)
    axes[1].tick_params(axis="x", length=0)
    style(axes[1], "", "时间收敛阶（观测值）", "(b) $D(C)$ 处理方式对时间精度的影响")
    leg = axes[1].legend(loc="upper left")
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.subplots_adjust(wspace=0.32)
    save(fig, "fig6_离散格式对比.pdf")


if __name__ == "__main__":
    print("生成图表：")
    fig_profiles()
    fig_history()
    fig_analytic()
    fig_convergence()
    fig_format()
    print("完成")
