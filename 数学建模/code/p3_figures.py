# -*- coding: utf-8 -*-
"""
问题 3 配图（数据来源：results/p3_fields.npz，由 p3_run.py 生成）：

  fig12_自适应步长.pdf   —— 宏观步长 H(t)、子步数 m(t)、误差估计 err(t) 的演化
  fig13_全烘干剖面.pdf   —— 烘干全过程浓度剖面 + 特征位置时程（中心/表面/0.5 cm）
  fig14_终止判定.pdf     —— max_r C(t) 的尾段放大、二分搜索、QS 外推与误差
  fig15_干燥前沿.pdf     —— r–t 平面上的等浓度线（干燥前沿推进）

配色沿用 p1/p2_figures.py 的设计令牌。
"""
import os
import sys
# 控制台按 UTF-8 输出（Windows GBK 终端下避免 UnicodeEncodeError 与乱码）
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

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


def load():
    p = os.path.join(RES, "p3_fields.npz")
    if not os.path.exists(p):
        raise SystemExit("缺少 %s，请先运行 p3_run.py" % p)
    return np.load(p, allow_pickle=True)


# ===========================================================================
def fig_adaptive(d):
    """宏观步长、子步数与误差估计的演化（创新点 1、2）。"""
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))

    t0 = d["log_t"] / 3600.0
    H = d["log_h"] / 60.0
    m = d["log_m"]
    err = d["log_err"]

    ax = axes[0]
    ax.plot(t0, H, color=CAT[0], linewidth=1.2)
    ax.axvline(4.0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
    ax.text(4.25, min(H) + 0.35 * (max(H) - min(H)), "预热期结束\nt = 4 h",
            color=MUTED, fontsize=7.5, va="bottom")
    style(ax, "时间 t / h", "宏观步长 H(t) / min", "（a）多时间尺度自适应步长")
    ax.set_ylim(0, max(H) * 1.15)
    ax.yaxis.set_major_locator(MultipleLocator(1.0))

    ax = axes[1]
    ax.plot(t0, m, color=CAT[1], linewidth=1.2)
    ax.axhline(1.0, color=AXIS, linewidth=0.9)
    style(ax, "时间 t / h", "子步数 m", "（b）步长加倍法的子步数")
    ax.set_ylim(0.5, max(m) * 1.15)

    ax = axes[2]
    flo = 1e-16
    ax.semilogy(t0, np.maximum(err, flo), color=CAT[2], linewidth=1.2)
    ax.axhline(1e-8, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.text(0.4, 2e-8, "目标容差 1e-8", color=CRIT, fontsize=7.5, va="bottom")
    over = err > 1e-8
    if np.any(over):
        ax.scatter(t0[over], np.maximum(err[over], flo), s=10, facecolors="none",
                   edgecolors=CRIT, linewidths=0.8, label="超限（按子步上限接受）")
        ax.legend(loc="lower left", handlelength=1.1, labelspacing=0.25)
    # 保留超出目标容差的初始瞬变段，避免被坐标上限裁掉
    ax.set_ylim(1e-16, max(1e-4, float(np.nanmax(err)) * 1.4))
    style(ax, "时间 t / h", "误差估计 err", "（c）步长加倍法误差估计")
    ax.text(0.97, 0.95, "N = %d；%d 步超限" % (len(d["r"]) - 1, over.sum()),
            transform=ax.transAxes, color=INK2, fontsize=7.5, ha="right", va="top")

    fig.tight_layout()
    save(fig, "fig12_自适应步长.pdf")


def fig_profiles(d):
    """全烘干过程的浓度剖面 + 特征位置时程。"""
    t = d["times"] / 3600.0
    r = d["r"] * 100.0
    C = d["C"]
    t_dry = float(d["t_dry"]) / 3600.0

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))

    # (a) 浓度剖面的时空演化
    ax = axes[0]
    marks = [0.0, 1.0, 3.0, 6.0, 12.0, 24.0, 36.0, 48.0, 56.0]
    for k, th in enumerate(marks):
        i = int(np.searchsorted(t, th))
        i = min(i, len(t) - 1)
        ax.plot(r, C[i], color=SEQ7[min(k * 7 // len(marks), 6)],
                linewidth=1.6 if th in (0.0, 56.0) else 1.0,
                label="t = %.0f h" % t[i])
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.text(0.10, 0.22, "判据 0.15", color=CRIT, fontsize=7.5, va="bottom")
    style(ax, "到中心距离 r / cm", "水分浓度 C / (kg/kg)", "（a）浓度剖面的演化")
    ax.legend(loc="upper right", ncol=1, handlelength=1.1, labelspacing=0.25)
    ax.set_ylim(0, 2.75)

    # (b) 特征位置的时程
    ax = axes[1]
    j05 = int(np.argmin(np.abs(r - 0.5)))
    j10 = int(np.argmin(np.abs(r - 1.0)))
    j20 = len(r) - 1
    for j, lab, c in ((0, "中心 r = 0", CAT[0]), (j05, "r = 0.5 cm", CAT[1]),
                      (j10, "r = 1.0 cm", CAT[2]), (j20, "表面 r = 2 cm", CAT[3])):
        ax.plot(t, C[:, j], color=c, linewidth=1.7, label=lab)
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.axvline(t_dry, color=INK, linewidth=1.0, linestyle=(0, (2, 2)))
    ax.text(t_dry - 1.0, 3.6, "t_num(N=640) = %.4f h" % t_dry, color=INK,
            fontsize=7.5, ha="right")
    style(ax, "时间 t / h", "水分浓度 C / (kg/kg)", "（b）特征位置的水分浓度")
    ax.set_xlim(0, 61)
    ax.legend(loc="lower left", handlelength=1.1, labelspacing=0.25)
    ax.set_yscale("log")
    ax.set_ylim(3e-2, 9)

    # (c) 浓度梯度的径向分布；量值按 cm 计算
    ax = axes[2]
    for k, th in enumerate([6.0, 24.0, 40.0, 50.0, 56.0]):
        i = min(int(np.searchsorted(t, th)), len(t) - 1)
        # dC/dr 以 cm 为横坐标，避免将 m^-1 误标为 cm^-1
        g = np.gradient(C[i], r)
        ax.plot(r, np.abs(g), color=SEQ7[k + 1], linewidth=1.4,
                label="t = %.0f h" % t[i])
    style(ax, "到中心距离 r / cm", r"$|\partial C/\partial r|$ / (kg/kg per cm)",
          "（c）浓度梯度的径向分布")
    ax.set_yscale("log")
    ax.legend(loc="lower left", handlelength=1.1, labelspacing=0.25)

    fig.tight_layout()
    save(fig, "fig13_全烘干剖面.pdf")


def fig_termination(d):
    """终止时间判定：尾段交点、二分区间和 QS 事后降阶检验。"""
    t = d["times"] / 3600.0
    C = d["C"]
    mx = C.max(axis=1)
    t_dry = float(d["t_dry"]) / 3600.0

    fig = plt.figure(figsize=(11.6, 3.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.05], wspace=0.46)

    # (a) 尾段 mx(t) 与判据的交点
    ax = fig.add_subplot(gs[0, 0])
    m = t >= 30.0
    ax.plot(t[m], mx[m], color=CAT[0], linewidth=1.8)
    ax.axhline(0.15, color=CRIT, linewidth=1.3, linestyle=(0, (4, 3)))
    ax.plot([t_dry], [0.15], marker="o", color=CRIT, markersize=5.5, zorder=5)
    ax.annotate("N = 640\n$t_{\\mathrm{num}}$ = %.4f h" % t_dry,
                xy=(t_dry, 0.15), xytext=(39.0, 0.166),
                color=INK, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=INK, linewidth=1.0))
    style(ax, "时间 t / h", r"$\max_r\,C$ / (kg/kg)", "（a）判据的尾段交点")
    ax.set_ylim(0.140, 0.20)
    ax.set_xlim(30, 58)

    # (b) 二分搜索的收敛过程
    ax = fig.add_subplot(gs[0, 1])
    hist = d["bisect_hist"] if "bisect_hist" in d else None
    if hist is not None and len(hist) > 1:
        ht = np.array([h[0] for h in hist], dtype=float)
        k = np.arange(1, len(ht) + 1)
        # 严格二分时相邻中点距离等于上一轮更新后区间宽度的一半。
        # 由保存的中点历史恢复实际区间，不使用浓度偏差冒充单调收敛量。
        widths_s = np.r_[2.0 * np.abs(np.diff(ht)), abs(ht[-1] - ht[-2])]
        if not np.allclose(widths_s[1:], widths_s[:-1] / 2.0):
            raise ValueError("二分中点历史不满足逐次折半，不能重建区间宽度")
        ax.semilogy(k, widths_s, marker="o", markersize=3.5,
                    color=CAT[2], linewidth=1.2)
        ax.axhline(1.0, color=CRIT, linewidth=1.0, linestyle=(0, (4, 3)))
        ax.text(1.15, 1.2, "终止容差 1 s", color=CRIT, fontsize=7.5)
        ax.set_xticks(k)
        ax.set_xticklabels(["%d" % v for v in k])
        style(ax, "已完成二分次数", "更新后区间宽度 / s", "（b）二分区间逐次缩小")
        ax.set_yscale("log")
    else:
        raise ValueError("缺少足够的二分中点历史")

    # (c) 冻结形状用 t0 时刻形状；漂移修正用 t >= t0 的完整后续轨道拟合
    ax = fig.add_subplot(gs[0, 2])
    xerr = d["qs_starts_h"] if "qs_starts_h" in d else None
    if xerr is not None:
        ef = np.abs(d["qs_err_frz"])
        ed = np.abs(d["qs_err_drf"])
        fin = np.isfinite(ef)
        ax.semilogy(xerr[fin], np.maximum(ef[fin], 1e-12), marker="s",
                    markersize=3.5, color=CAT[1], linewidth=1.2, label="QS 冻结形状")
        ax.semilogy(xerr, np.maximum(ed, 1e-12), marker="o",
                    markersize=3.5, color=CAT[2], linewidth=1.2, label="QS 漂移（后续轨道拟合）")
        if np.any(~fin):
            ax.plot(xerr[~fin], np.full(np.sum(~fin), 1700.0), linestyle="none",
                    marker="x", markersize=5, color=CAT[1], label="冻结形状无有限交点")
        ax.axhline(1.0, color=CRIT, linewidth=1.1, linestyle=(0, (4, 3)))
        ax.text(30.5, 1.3, "1 min", color=CRIT, fontsize=7.5)
        ax.legend(loc="upper left", bbox_to_anchor=(0, 0.90), handlelength=1.2,
                  labelspacing=0.3, fontsize=7.2)
        ax.set_ylim(3e-1, 3e3)
    style(ax, "起始时刻 $t_0$ / h", "|降阶结果 − 二分结果| / min", "（c）事后降阶检验")

    save(fig, "fig14_终止判定.pdf")


def fig_front(d):
    """r–t 平面上的等浓度线（干燥前沿）。"""
    t = d["times"] / 3600.0
    r = d["r"] * 100.0
    C = d["C"]
    C = np.where(C > 1e-12, C, 1e-12)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))
    ax = axes[0]
    levels = [0.05, 0.15, 0.3, 0.6, 1.0, 1.5, 2.0, 2.6]
    cs = ax.contourf(t, r, C.T, levels=levels, cmap="Blues", extend="both")
    cb = fig.colorbar(cs, ax=ax, pad=0.02, ticks=levels)
    cb.set_label("C / (kg/kg)", color=INK2)
    cb.outline.set_edgecolor(AXIS)
    ax.set_yticks(np.arange(0, 2.01, 0.25))
    ax.set_yticklabels(["%g" % v for v in np.arange(0, 2.01, 0.25)])
    ax.contour(t, r, C.T, levels=[0.15], colors=[CRIT], linewidths=1.8)
    ax.text(21.0, 1.62, "C = 0.15 等值线", color=CRIT, fontsize=7.8)
    style(ax, "时间 t / h", "到中心距离 r / cm", "（a）干燥前沿在 r–t 平面上的推进")
    t_num = float(d["t_dry"]) / 3600.0
    ax.set_xlim(0, t_num)

    # (b) 各半径的达标时刻
    ax = axes[1]
    reach = np.full(C.shape[1], np.nan)
    for j in range(C.shape[1]):
        idx = np.where(C[:, j] < 0.15)[0]
        if len(idx):
            i = idx[0]
            if i == 0:
                reach[j] = t[0]
            else:
                # 60 s 输出采样之间作线性过阈值插值，得到一致近似时刻
                c0, c1 = C[i - 1, j], C[i, j]
                reach[j] = t[i - 1] + (0.15 - c0) * (t[i] - t[i - 1]) / (c1 - c0)
    ax.plot(r, reach, color=CAT[0], linewidth=1.8)
    ax.axhline(float(d["t_dry"]) / 3600.0, color=CRIT, linewidth=1.3,
               linestyle=(0, (4, 3)))
    ax.text(1.95, float(d["t_dry"]) / 3600.0 - 2.2,
            "t_num(N=640) = %.4f h" % (float(d["t_dry"]) / 3600.0),
            color=CRIT, fontsize=7.8, ha="right", va="top")
    style(ax, "到中心距离 r / cm", "首次达标时刻 / h",
          "（b）各处首次达标时刻")
    ax.text(0.03, 0.05, "60 s 输出样本之间作线性过阈值插值",
            transform=ax.transAxes, color=INK2, fontsize=7.5)
    ax.set_ylim(0, 60)
    ax.set_xlim(0, 2)

    fig.tight_layout()
    save(fig, "fig15_干燥前沿.pdf")


# ===========================================================================
def main():
    d = load()
    print("读取 results/p3_fields.npz")
    fig_adaptive(d)
    fig_profiles(d)
    fig_termination(d)
    fig_front(d)


if __name__ == "__main__":
    main()
