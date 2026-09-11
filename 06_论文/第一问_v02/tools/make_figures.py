"""Reproduce the Q1 paper figure from the final full-precision solution.

Run from any directory. Optional --source and --output-dir override the
repository-relative canonical input and the adjacent figures directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


def main() -> None:
    paper_dir = Path(__file__).resolve().parents[1]
    repo_dir = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path,
        default=repo_dir / "04_结果" / "第一问" / "q1_main_r16.npz",
    )
    parser.add_argument("--output-dir", type=Path, default=paper_dir / "figures")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    available = {font.name for font in font_manager.fontManager.ttflist}
    chinese_font = next(
        (name for name in ("SimSun", "Noto Serif CJK SC", "Microsoft YaHei", "SimHei") if name in available),
        None,
    )
    if chinese_font is None:
        raise RuntimeError("Install SimSun, Noto Serif CJK SC, Microsoft YaHei, or SimHei.")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [chinese_font, "DejaVu Sans"],
        "font.size": 9.5,
        "axes.labelsize": 9.5,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.unicode_minus": False,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
        "savefig.facecolor": "white",
    })

    times = [100, 300, 600, 900, 1200, 1500, 1800]
    with np.load(args.source) as data:
        r_cm = data["r"] * 100
        z = data["z"]
        temperature = data["snapshots_t"][:, 0, :]
        moisture = data["snapshots_c"][:, 0, :]
    if temperature.shape != (7, len(r_cm)) or moisture.shape != temperature.shape:
        raise ValueError("Expected seven saved time snapshots on the radial mesh.")
    if not np.isclose(z[0], 0) or not np.isclose(r_cm[-1], 2):
        raise ValueError("Expected the z=0 midplane and a 2 cm cylinder radius.")

    colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#222222"]
    styles = ["-", (0, (5, 2)), (0, (1.5, 1.5)), (0, (6, 2, 1.5, 2)),
              (0, (8, 2, 2, 2)), (0, (3, 1.5, 1, 1.5)), (0, (9, 2))]
    fig, axes = plt.subplots(1, 2, figsize=(16 / 2.54, 6.4 / 2.54))
    fig.subplots_adjust(left=0.087, right=0.977, bottom=0.29, top=0.87, wspace=0.30)
    for ax, values in zip(axes, (temperature, moisture)):
        for index, time in enumerate(times):
            ax.plot(r_cm, values[index], color=colors[index], linestyle=styles[index],
                    linewidth=1.15, label=f"{time} s")
        ax.set_xlim(0, 2)
        ax.set_xticks([0, 0.5, 1, 1.5, 2])
        ax.set_xlabel("径向位置 r / cm", labelpad=3)
        ax.tick_params(direction="in", length=3, width=0.65, pad=3)
        ax.grid(axis="y", color="#d7d7d7", linewidth=0.45, alpha=0.7)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylim(28, 38)
    axes[0].set_yticks([28, 30, 32, 34, 36, 38])
    axes[0].set_ylabel("温度 T / ℃", labelpad=3)
    axes[0].set_title("(a) 中截面温度分布", pad=6)
    axes[1].set_ylim(1.45, 2.6)
    axes[1].set_yticks([1.5, 1.75, 2, 2.25, 2.5])
    axes[1].set_ylabel("含水率 C / (kg/kg)", labelpad=3)
    axes[1].set_title("(b) 中截面含水率分布", pad=6)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=7, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, 0.035), handlelength=2.0,
               handletextpad=0.45, columnspacing=0.75, borderaxespad=0)
    for extension in ("pdf", "png"):
        output = args.output_dir / f"q1_radial_profiles.{extension}"
        fig.savefig(output, dpi=400)
        print(output)
    plt.close(fig)
    print(f"Source: {args.source.resolve()}; snapshots: {times}; z=0; font: {chinese_font}")


if __name__ == "__main__":
    main()
