# 第一问 LaTeX 论文稿

主文件：`0910_第一问_药材预热模型_v01_Codex.tex`。正文仅覆盖药材烘干问题的第一问，采用二维轴对称模型，以中截面 z=0 输出径向结果。

## 编译

安装 TeX Live（推荐2025或更新版本）或带有 ctex、Fandol、newtxmath 的其他 TeX 发行版。编译字体为随 TeX Live 提供的 Fandol 和 TeX Gyre Termes，无需 Windows 专有中文字体。

在本目录执行两遍：

```text
xelatex -interaction=nonstopmode -halt-on-error 0910_第一问_药材预热模型_v01_Codex.tex
xelatex -interaction=nonstopmode -halt-on-error 0910_第一问_药材预热模型_v01_Codex.tex
```

Windows PowerShell 也可运行 `./build.ps1`，中间文件会保存于 build，最终 PDF 复制至本目录。上传整个源码包到 Overleaf 后，选择 XeLaTeX 和上述主文件即可编译。图件已提供 PDF 矢量文件，不需要 Python 才能编译论文。

## 文件结构

- `sections/q1.tex`：正文、公式及模型论述。
- `tables/`：自动生成的四张表和正文数值宏。
- `figures/`：正式径向结果图，矢量 PDF 和400 dpi PNG。
- `tools/make_tables.py`：在原仓库中读取权威 JSON，生成表格与数值快照；不重新运行模型。
- `tools/make_figures.py`：读取原仓库的 `q1_main_r16.npz` 生成图；可用 `--source` 指定输入。
- `evidence/`：写作依据、逐项数值溯源和验证记录，不进入论文正文。

在原仓库中，使用安装 NumPy、Matplotlib 的 Python 执行绘图脚本即可重绘。模型原始数组和完整 result1.xlsx 仍在仓库 `04_结果/第一问/`，不为编译论文而复制大量中间计算文件。

本稿按用户要求直接使用已核验第一问内容进行写作。数值快照是本稿溯源依据，不伪称为用户签收或整场比赛的最终提交冻结。论文已明确经验含水率边界、未纳入机制及数值检验范围。
