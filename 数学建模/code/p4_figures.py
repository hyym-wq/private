# -*- coding: utf-8 -*-
"""
问题 4 配图（数据来源：results/p4_fields.npz，由 p4_run.py 生成）：

  fig16_半径收缩与对流.pdf —— 附件 2 半径 R(t)、中心差分导数 R'(t)（收缩集中在前 24 h）
  fig17_移动边界剖面.pdf   —— 物理坐标 r 剖面（移动边界）+ 归一化坐标 E 剖面（浓缩效应）
  fig18_特征位置时程.pdf   —— 中心/表面/固定距离时程 + 有收缩 vs 无收缩对比
  fig19_验证与对比.pdf     —— 空间收敛、对流项符号敏感性、烘干时长对比

配色沿用 p1/p2/p3_figures.py 的设计令牌。
"""
import os
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_model as m

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
    p = os.path.join(RES, "p4_fields.npz")
    if not os.path.exists(p):
        raise SystemExit("缺少 %s，请先运行 p4_run.py" % p)
    return np.load(p, allow_pickle=True)


def _profile_phys(d, i):
    """第 i 个输出时刻的物理坐标剖面 (r_cm, C)。"""
    E = d["E"]
    R = d["R_t"][i]
    return E * R * 100.0, d["C"][i]


def _fixed_position_series(d, rcm):
    """固定距离 rcm 处的水分浓度时程（材料收缩离开后置 NaN）。"""
    E = d["E"]; C = d["C"]; R_t = d["R_t"]; times = d["times"]
    r = rcm * 0.01
    out = np.full(len(times), np.nan)
    for i in range(len(times)):
        if r < R_t[i]:
            out[i] = np.interp(r / R_t[i], E, C[i])
    return out


# ===========================================================================
def fig_radius(d):
    """附件 2 半径与中心差分导数：收缩集中在前 ~24 h。"""
    tR = d["tR"] / 3600.0
    R = d["R_data"] * 100.0          # cm
    Rp = d["Rp_data"] * 3600.0 * 100.0   # cm/h（中心差分）

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.2))

    ax = axes[0]
    ax.plot(tR, R, color=CAT[0], linewidth=1.8)
    ax.axhline(1.198, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
    ax.text(1.0, 1.21, "收缩停止 ~1.198 cm", color=MUTED, fontsize=7.5, va="bottom")
    ax.axvline(24.0, color=AXIS, linewidth=0.9, linestyle=(0, (4, 3)))
    ax.text(24.5, 1.85, "t ≈ 24 h\n收缩基本停止", color=MUTED, fontsize=7.5, va="center")
    style(ax, "时间 t / h", "半径 R(t) / cm", "（a）附件 2 的收缩半径")
    ax.set_ylim(1.0, 2.15)

    ax = axes[1]
    ax.plot(tR, Rp, color=CAT[1], linewidth=1.8)
    ax.axhline(0.0, color=AXIS, linewidth=0.9)
    ax.fill_between(tR, Rp, 0, where=Rp < 0, color=CAT[1], alpha=0.10)
    style(ax, "时间 t / h", r"$R^\prime(t)$ / (cm/h)", "（b）中心差分数值微分")
    ax.set_ylim(min(Rp) * 1.15, max(Rp) * 0.15)

    fig.tight_layout()
    save(fig, "fig16_半径收缩与对流.pdf")


def fig_profiles(d):
    """物理坐标剖面（移动边界）与归一化坐标剖面（浓缩效应）。"""
    times = d["times"] / 3600.0
    E = d["E"]; C = d["C"]; R_t = d["R_t"]
    t_dry = float(d["times"][-1]) / 3600.0

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))

    # (a) 物理坐标：右端在 R(t) 处收缩
    ax = axes[0]
    marks = [0.0, 3.0, 6.0, 12.0, 24.0, 36.0, 48.0, 52.4]
    for k, th in enumerate(marks):
        i = min(int(np.searchsorted(times, th)), len(times) - 1)
        r, ci = _profile_phys(d, i)
        ax.plot(r, ci, color=SEQ7[min(k * 7 // len(marks), 6)],
                linewidth=1.6 if th in (0.0, 52.4) else 1.0,
                label="t = %.1f h" % times[i])
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.text(0.06, 0.26, "判据 0.15", color=CRIT, fontsize=7.5, va="bottom")
    style(ax, "到中心距离 r / cm", "水分浓度 C / (kg/kg)", "（a）物理坐标剖面（右端随收缩内移）")
    ax.set_ylim(0, 2.75)
    ax.set_xlim(0, 2.1)
    ax.legend(loc="upper right", ncol=1, handlelength=1.1, labelspacing=0.25, fontsize=7)

    # (b) 归一化坐标：边界恒在 E=1，剖面因对流被"抬升"
    ax = axes[1]
    for k, th in enumerate(marks):
        i = min(int(np.searchsorted(times, th)), len(times) - 1)
        ax.plot(E, C[i], color=SEQ7[min(k * 7 // len(marks), 6)],
                linewidth=1.6 if th in (0.0, 52.4) else 1.0,
                label="t = %.1f h" % times[i])
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    style(ax, "归一化坐标 E = r / R(t)", "水分浓度 C / (kg/kg)",
          "（b）归一化坐标剖面（固定域 [0,1]）")
    ax.set_ylim(0, 2.75)
    ax.set_xlim(0, 1.0)
    ax.legend(loc="upper right", ncol=1, handlelength=1.1, labelspacing=0.25, fontsize=7)

    fig.tight_layout()
    save(fig, "fig17_移动边界剖面.pdf")


def fig_timeseries(d, N_fixed=400):
    """特征位置时程 + 有/无收缩对比。"""
    times = d["times"] / 3600.0
    C = d["C"]
    t_dry = float(d["times"][-1]) / 3600.0

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))

    # (a) 中心 / 表面 / 固定距离（收缩离开后截断）
    ax = axes[0]
    ax.plot(times, C[:, 0], color=CAT[0], linewidth=1.7, label="中心 r = 0")
    ax.plot(times, C[:, -1], color=CAT[3], linewidth=1.7, label="表面 r = R(t)")
    for rcm, c in ((0.5, CAT[1]), (1.0, CAT[2])):
        s = _fixed_position_series(d, rcm)
        msk = ~np.isnan(s)
        ax.plot(times[msk], s[msk], color=c, linewidth=1.4,
                label="r = %.1f cm（收缩离开后截断）" % rcm)
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.axvline(t_dry, color=INK, linewidth=1.0, linestyle=(0, (2, 2)))
    ax.text(t_dry - 1.0, 2.5, "$t_{\\mathrm{dry}}$ = %.4f h" % t_dry,
            color=INK, fontsize=7.5, ha="right")
    style(ax, "时间 t / h", "水分浓度 C / (kg/kg)", "（a）特征位置的水分浓度")
    ax.set_yscale("log")
    ax.set_ylim(3e-2, 4)
    ax.set_xlim(0, 54)
    ax.legend(loc="lower left", handlelength=1.1, labelspacing=0.25, fontsize=7)

    # (b) 有收缩 vs 无收缩（中心浓度）
    ax = axes[1]
    ax.plot(times, C[:, 0], color=CAT[0], linewidth=1.7, label="有收缩（附件 2）")

    def Rf(t): return 0.02
    def Rpf(t): return 0.0
    T_inf, C_inf = m.make_ambient()
    res_f = m.solve(T_inf, C_inf, Rf, Rpf, N=N_fixed, dt=60.0, t_end=700000.0)
    tf = res_f["times"] / 3600.0
    ax.plot(tf, res_f["C"][:, 0], color=CAT[2], linewidth=1.7, label="无收缩（R 冻结 2 cm）")
    ax.axhline(0.15, color=CRIT, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.axvline(t_dry, color=CAT[0], linewidth=1.0, linestyle=(0, (2, 2)))
    ax.axvline(tf[-1], color=CAT[2], linewidth=1.0, linestyle=(0, (2, 2)))
    ax.text(t_dry + 0.5, 1.8, "%.1f h" % t_dry, color=CAT[0], fontsize=7.5)
    ax.text(tf[-1] - 2.0, 1.8, "%.1f h" % tf[-1], color=CAT[2], fontsize=7.5, ha="right")
    style(ax, "时间 t / h", "中心水分浓度 C(0,t) / (kg/kg)", "（b）收缩显著缩短烘干时长")
    ax.set_yscale("log")
    ax.set_ylim(0.1, 3)
    ax.set_xlim(0, 135)
    ax.legend(loc="lower left", handlelength=1.1, labelspacing=0.25, fontsize=7)

    fig.tight_layout()
    save(fig, "fig18_特征位置时程.pdf")


def fig_verify(d):
    """空间收敛、对流项符号敏感性、烘干时长对比（数据取自 problem4_verification.txt）。"""
    fig = plt.figure(figsize=(11.4, 3.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], wspace=0.42)

    # (a) 空间收敛
    ax = fig.add_subplot(gs[0, 0])
    Ns = np.array([100, 200, 400, 800])
    td = np.array([52.3667, 52.3833, 52.4000, 52.4000])
    ax.semilogx(Ns, td, marker="o", color=CAT[0], linewidth=1.6)
    ax.axhline(52.40, color=CRIT, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.text(110, 52.4008, "收敛值 52.40 h", color=CRIT, fontsize=7.2, va="bottom")
    style(ax, "网格节点数 N（对数）", r"$t_{\mathrm{dry}}$ / h", "（a）空间收敛")
    ax.set_ylim(52.35, 52.42)
    ax.set_yticks([52.36, 52.38, 52.40, 52.42])

    # (b) 对流项符号敏感性
    ax = fig.add_subplot(gs[0, 1])
    labels = ["−1（错误）", "0（无对流）", "+1（正确）"]
    vals = [48.650, 50.833, 52.400]
    cols = [MUTED, CAT[2], CAT[0]]
    ax.bar(labels, vals, color=cols, width=0.55)
    for x, v in enumerate(vals):
        ax.text(x, v + 0.25, "%.1f" % v, ha="center", fontsize=8, color=INK2)
    style(ax, "对流项 (E R′/R)·∂/∂E 的符号", r"$t_{\mathrm{dry}}$ / h",
          "（b）对流项符号的物理检验")
    ax.set_ylim(46, 54)

    # (c) 烘干时长对比
    ax = fig.add_subplot(gs[0, 2])
    cats = ["问题 4\n有收缩", "问题 4\n无收缩", "问题 3\n无收缩\n（附录 3）"]
    vals = [52.4, 129.1, 57.2]
    cols = [CAT[0], CAT[2], MUTED]
    ax.bar(cats, vals, color=cols, width=0.55)
    for x, v in enumerate(vals):
        ax.text(x, v + 2.5, "%.1f h" % v, ha="center", fontsize=8, color=INK2)
    style(ax, "", r"$t_{\mathrm{dry}}$ / h", "（c）收缩对烘干时长的缩短")
    ax.set_ylim(0, 145)
    ax.tick_params(axis="x", labelsize=7.5)

    fig.tight_layout()
    save(fig, "fig19_验证与对比.pdf")


# ===========================================================================
def main():
    d = load()
    print("读取 results/p4_fields.npz")
    fig_radius(d)
    fig_profiles(d)
    fig_timeseries(d)
    fig_verify(d)


if __name__ == "__main__":
    main()
