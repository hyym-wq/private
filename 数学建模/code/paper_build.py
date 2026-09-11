# -*- coding: utf-8 -*-
"""
生成 2026 高教社杯 A 题「药材的烘干问题」完整论文（DOCX，CUMCM 保留模板）。

从 reference.docx 克隆模板，写入全部四个问题的模型与结果，
插入三线表、公式（含上下标）与真实结果图（PDF → PNG），输出到 paper/。
"""
import os
import shutil

import docx
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_TAB_ALIGNMENT, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import fitz  # PyMuPDF

BASE = r"C:\Users\lenovo\Desktop\数学建模"
FIG = os.path.join(BASE, "figures")
PNG = os.path.join(BASE, "figures", "_paper_png")
os.makedirs(PNG, exist_ok=True)
TEMPLATE = r"C:\Users\lenovo\.claude\skills\math-competition\writing\cumcm-template\assets\reference.docx"
OUT = os.path.join(BASE, "paper", "药材的烘干问题.docx")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

SONG = "宋体"
HEI = "黑体"
LATIN = "Times New Roman"
MONO = "Courier New"


# ---------------------------------------------------------------------------
# 基础格式工具
# ---------------------------------------------------------------------------
def fmt(run, east=SONG, latin=LATIN, size=12, bold=False):
    run.font.name = latin
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn('w:rFonts'))
    if rfonts is None:
        rfonts = OxmlElement('w:rFonts')
        rpr.append(rfonts)
    rfonts.set(qn('w:ascii'), latin)
    rfonts.set(qn('w:hAnsi'), latin)
    rfonts.set(qn('w:eastAsia'), east)
    run.font.size = Pt(size)
    run.font.bold = bold


def body(doc, text, indent=True):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing = 1.0
    pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    if indent:
        pf.first_line_indent = Pt(24)  # 首行缩进两个汉字
    add_math_runs(p, text)
    return p


def title(doc, text):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf.space_before = Pt(18)
    pf.space_after = Pt(18)
    r = p.add_run(text)
    fmt(r, east=HEI, size=16, bold=True)
    return p


def h1(doc, text):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf.space_before = Pt(12)
    pf.space_after = Pt(12)
    r = p.add_run(text)
    fmt(r, east=HEI, size=14, bold=True)
    return p


def h2(doc, text):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = WD_ALIGN_PARAGRAPH.LEFT
    pf.space_before = Pt(6)
    pf.space_after = Pt(6)
    r = p.add_run(text)
    fmt(r, east=HEI, size=12, bold=True)
    return p


def h3(doc, text):
    return h2(doc, text)


def add_math_runs(p, text, east=SONG, latin=LATIN, size=12, bold=False):
    """将 _{...} 与 ^{...} 分别解析为下标/上标 run，其余为普通 run。"""
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch in ('_', '^') and i + 1 < n and text[i + 1] == '{':
            j = text.index('}', i + 2)
            r = p.add_run(text[i + 2:j])
            fmt(r, east=east, latin=latin, size=size, bold=bold)
            if ch == '_':
                r.font.subscript = True
            else:
                r.font.superscript = True
            i = j + 1
        else:
            j = i
            while j < n and text[j] not in '_^{':
                j += 1
            if j == i:
                j += 1
            r = p.add_run(text[i:j])
            fmt(r, east=east, latin=latin, size=size, bold=bold)
            i = j


def eq(doc, text, num=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(6)
    pf.space_after = Pt(6)
    pf.tab_stops.add_tab_stop(Cm(8.0), WD_TAB_ALIGNMENT.CENTER)
    pf.tab_stops.add_tab_stop(Cm(16.0), WD_TAB_ALIGNMENT.RIGHT)
    r0 = p.add_run("\t")
    fmt(r0)
    add_math_runs(p, text)
    if num:
        r1 = p.add_run("\t(%s)" % num)
        fmt(r1)
    return p


def page_break(doc):
    p = doc.add_paragraph()
    p.add_run().add_break(WD_BREAK.PAGE)


# ---------------------------------------------------------------------------
# 三线表
# ---------------------------------------------------------------------------
def _table_borders_none(tbl):
    tblPr = tbl._tbl.tblPr
    borders = OxmlElement('w:tblBorders')
    for edge in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        el = OxmlElement('w:' + edge)
        el.set(qn('w:val'), 'nil')
        borders.append(el)
    tblPr.append(borders)


def _set_cell_border(cell, top=None, bottom=None):
    tcPr = cell._tc.get_or_add_tcPr()
    borders = tcPr.find(qn('w:tcBorders'))
    if borders is None:
        borders = OxmlElement('w:tcBorders')
        tcPr.append(borders)
    for edge, spec in (('top', top), ('bottom', bottom)):
        if spec is None:
            continue
        val, sz = spec
        el = borders.find(qn('w:' + edge))
        if el is None:
            el = OxmlElement('w:' + edge)
            borders.append(el)
        el.set(qn('w:val'), val)
        el.set(qn('w:sz'), str(sz))
        el.set(qn('w:color'), '000000')
        el.set(qn('w:space'), '0')


def add_table(doc, caption, header, rows, note=None, first_col_left=False):
    cp = doc.add_paragraph()
    cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cp.paragraph_format.space_before = Pt(6)
    cp.paragraph_format.space_after = Pt(3)
    r = cp.add_run(caption)
    fmt(r, east=HEI, size=10.5, bold=True)

    ncol = len(header)
    nrow = len(rows) + 1
    t = doc.add_table(rows=nrow, cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    _table_borders_none(t)

    for j, h in enumerate(header):
        c = t.cell(0, j)
        c.text = ""
        p = c.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        add_math_runs(p, str(h), east=HEI, size=10.5, bold=True)
        _set_cell_border(c, top=('single', 12), bottom=('single', 6))

    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            c = t.cell(i + 1, j)
            c.text = ""
            p = c.paragraphs[0]
            p.alignment = (WD_ALIGN_PARAGRAPH.LEFT if (first_col_left and j == 0)
                           else WD_ALIGN_PARAGRAPH.CENTER)
            add_math_runs(p, str(val), east=SONG, size=10.5)

    for j in range(ncol):
        _set_cell_border(t.cell(nrow - 1, j), bottom=('single', 12))

    if note:
        np_ = doc.add_paragraph()
        np_.paragraph_format.space_before = Pt(2)
        np_.paragraph_format.space_after = Pt(6)
        r = np_.add_run(note)
        fmt(r, east=SONG, size=9)
    return t


# ---------------------------------------------------------------------------
# 图
# ---------------------------------------------------------------------------
def _pdf2png(name):
    src = os.path.join(FIG, name)
    base = os.path.splitext(name)[0]
    dst = os.path.join(PNG, base + ".png")
    if not os.path.exists(dst):
        d = fitz.open(src)
        pix = d[0].get_pixmap(dpi=200)
        pix.save(dst)
        d.close()
    return dst


def fig(doc, name, caption, width=14.0):
    png = _pdf2png(name)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(3)
    p.add_run().add_picture(png, width=Cm(width))
    cp = doc.add_paragraph()
    cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cp.paragraph_format.space_after = Pt(6)
    r = cp.add_run(caption)
    fmt(r, east=HEI, size=10.5)


def code_block(doc, text):
    for line in text.split("\n"):
        p = doc.add_paragraph()
        pf = p.paragraph_format
        pf.line_spacing = 1.0
        pf.space_before = Pt(0)
        pf.space_after = Pt(0)
        r = p.add_run(line if line else " ")
        fmt(r, east=SONG, latin=MONO, size=9)


# ---------------------------------------------------------------------------
# 正文
# ---------------------------------------------------------------------------
def build():
    shutil.copyfile(TEMPLATE, OUT)
    doc = Document(OUT)

    # 删除模板全部正文内容（保留分节设置与页脚页码）
    for child in list(doc.element.body):
        if child.tag != qn('w:sectPr'):
            doc.element.body.remove(child)

    # ===================== 第 1 页：题目 + 摘要 + 关键词 =====================
    title(doc, "药材的烘干问题")

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(12)
    r = p.add_run("摘  要")
    fmt(r, east=HEI, size=14, bold=True)

    body(doc, "针对中药材热风烘干过程中药材内部温度与水分浓度的时空演化规律及烘干时长的定量确定问题，本文从能量与质量的守恒机制出发，建立径向一维热质耦合偏微分方程模型，采用守恒型有限体积离散与隐式时间推进求解，并以多重数值检验验证了模型的正确性与结论的稳定性。", indent=False)

    body(doc, "针对问题一，利用圆柱药材长径比远大于一的几何特征将三维问题化简为沿半径方向的一维轴对称问题，按附录二的常物性参数建立热传导与水分扩散耦合模型，其中水分扩散系数随水分浓度呈指数衰减，环境温度与水分浓度由附件一分段线性插值给定；采用守恒型有限体积法离散空间、隐式格式离散时间并配合牛顿迭代与追赶法求解。1800 s 时中心温度 33.5752°C、表面温度 36.7853°C，中心水分浓度 2.5500、表面 1.5120 kg/kg。数值解与级数解析解误差在万分之一量级，空间与时间均呈二阶收敛，全局能量与质量守恒偏差处于十的负十四次方量级。", indent=False)

    body(doc, "针对问题二，将物性推广为随水分浓度和温度变化的非线性关系（密度、比热、导热系数与扩散系数均按附录三随状态变化），建立整个烘干过程的非线性变参数热质耦合模型。给出 3 h 内每隔 0.5 h 沿半径方向的温度与水分浓度分布，3 h 时中心水分浓度 1.7662、表面 1.0078 kg/kg，中心温度升至 49.85°C。", indent=False)

    body(doc, "针对问题三，在问题二模型基础上确定满足\u201c各处水分浓度均低于 0.15 kg/kg\u201d的烘干时长。针对干燥末期表面扩散边界层极薄导致均匀网格收敛缓慢的问题，引入向表面等比加密的分级网格与基于步长加倍误差控制的自适应时间步长，并配合二分法精确定位。求得烘干时长 57.17 h。网格收敛、时间步无关性、准稳态外推独立校核与全局质量守恒等检验表明该时长在 0.3 min 量级内稳定。", indent=False)

    body(doc, "针对问题四，考虑水分流失引起的药材尺寸收缩，由附件二给出半径随时间的收缩规律，引入坐标变换将移动边界映射为固定计算域，精确处理坐标变换产生的对流项，收缩速率由中心差分数值微分得到，网格随边界等比例缩放以保证空间分辨率一致，物性参数采用附录四。求得考虑收缩后的烘干时长 52.4 h，较固定边界情形缩短约 2.46 倍，表明收缩显著加速干燥。对流项符号的物理检验、退化一致性与守恒性校验验证了移动边界处理的正确性。", indent=False)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(6)
    r = p.add_run("关键词：")
    fmt(r, east=HEI, size=12, bold=True)
    r2 = p.add_run("热质耦合；烘干时长；有限体积法；分级网格；移动边界；坐标变换")
    fmt(r2, east=SONG, size=12)

    page_break(doc)

    # ===================== 一、问题重述 =====================
    h1(doc, "一、问题重述")
    h2(doc, "1.1 问题背景")
    body(doc, "干燥是决定中药材成品品质的关键工序之一，热风烘干通过调控烘房温湿环境完成干燥，主要包括预热平衡与恒温干燥两个阶段。工艺参数选取不当容易导致干燥效率低、能耗高、成品品质不稳定等问题；传统的试验优化模式存在成本高、周期长等缺陷，亟需借助数理分析与数值仿真的方法得到干燥规律，为工艺参数优化提供依据。")
    h2(doc, "1.2 问题重述")
    body(doc, "本题研究对象为长 25 cm、半径 2 cm 的圆柱形中药材。烘干开始时药材温度为 28°C、水分浓度（干基含水率）为 2.55 kg/kg，烘房温度与水分浓度随时间的变化由附件 1 给出，烘干过程中药材半径随时间的变化由附件 2 给出。需建立数学模型依次解决以下四个递进问题：")
    body(doc, "（1）问题一：建立预热平衡阶段药材温度与水分浓度变化规律的数学模型（相关参数见附录 2），按表 1、表 2 格式给出 100、300、600、900、1200、1500、1800 s、到药材中心 0、0.5、1、1.5、2 cm 处的结果，并将完整结果保存到 result1.xlsx。")
    body(doc, "（2）问题二：建立整个烘干过程药材温度与水分浓度变化规律的数学模型（经验公式统一采用附录 3），按表 3、表 4 格式给出 3 h 内每隔 0.5 h 的结果，并保存到 result2.xlsx。")
    body(doc, "（3）问题三：按\u201c药材各处水分浓度低于 0.15 kg/kg\u201d的烘干要求确定烘干时长（单位 h），按表 5 格式给出结果，并保存到 result3.xlsx。")
    body(doc, "（4）问题四：考虑水分流失导致的尺寸变化（附件 2 给出半径变化），确定烘干时长（相关经验公式见附录 4），按表 6 格式给出结果，并保存到 result4.xlsx。")

    # ===================== 二、问题分析 =====================
    h1(doc, "二、问题分析")
    body(doc, "四个问题本质上是同一物理过程——圆柱药材在热风环境下的热质耦合传递——在不同物性参数、不同时间尺度、不同边界条件下的递进细化。问题一处理预热平衡阶段（30 min）的短时行为，物性按常物性（附录 2）处理，此时温度场尚未完全均匀、水分仅在表面附近扩散；问题二将时间尺度扩展到整个烘干过程（2~3 天），物性随水分浓度与温度非线性变化（附录 3）；问题三在问题二基础上进一步要求精确确定烘干时长；问题四引入水分流失导致的尺寸收缩（移动边界，附件 2 + 附录 4），是前序模型在变形几何方向上的最终扩展。因此整个求解采用统一的\u201c径向一维热质耦合偏微分方程 + 守恒型有限体积离散 + 隐式时间推进\u201d框架，各问题按需替换物性关系、边界条件与网格策略。")
    h2(doc, "2.1 问题一的分析")
    body(doc, "预热平衡阶段时间短（30 min），由附录 2 可知密度、比热、导热系数均为常数，仅水分扩散系数随水分浓度变化。此时温度场与水分场通过物性弱耦合，但两场仍需同时求解。关键判断：热扩散系数远大于初始湿扩散系数（约 34 倍），故温度在数分钟内趋于均匀，而水分在 30 min 内仅表面薄层明显下降。边界条件上，表面为第三类对流边界（Robin），中心为轴对称零通量。环境温度与水分浓度由附件 1 分段线性插值得到。")
    h2(doc, "2.2 问题二的分析")
    body(doc, "整个烘干过程时间跨度 2~3 天，覆盖预热平衡与恒温干燥两个阶段。附录 3 中密度、比热、导热系数均随水分浓度变化，扩散系数同时随水分浓度与温度变化，热传导方程因此成为非线性方程，热与湿两场通过变物性强耦合。需在问题一格式基础上引入牛顿型线性化与迭代。恒温干燥阶段环境参数维持预热末值。该模型给出整个烘干过程的完整时空分布，为问题三的时长判定奠定基础。")
    h2(doc, "2.3 问题三的分析")
    body(doc, "问题三的关键是\u201c各处水分浓度均低于 0.15 kg/kg\u201d这一终止判据的精确实现。由问题二可知，干燥末期表面出现极薄的扩散边界层（厚度约 3e-5 m），均匀网格无法分辨该薄层，导致烘干时长对网格敏感、空间收敛缓慢；同时初期环境阶跃引起的瞬变要求小时间步，而漫长的恒温干燥段可用大步长。因此需引入向表面等比加密的分级网格与自适应时间步长，并对终止时刻用二分法精确定位，以保证烘干时长的精度与稳定性。")
    h2(doc, "2.4 问题四的分析")
    body(doc, "实际烘干中水分流失使药材收缩，半径随时间减小（附件 2），扩散路径缩短且浓缩效应使水分向内部对流输送，二者共同改变干燥进程，固定边界计算域不再适用，需引入坐标变换将移动边界映射为固定域。关键难点在于坐标变换产生的对流项符号与离散的正确性，以及收缩速率的获得（由附件 2 中心差分数值微分）。物性采用附录 4（与附录 3 形式相同、系数不同）。该模型是问题二、三在移动边界下的推广，退化（收缩速率为零）时应精确回到固定边界结果。")

    # ===================== 三、基本假设 =====================
    h1(doc, "三、基本假设")
    body(doc, "（1）药材为半径 2 cm、长 25 cm 的圆柱，长径比 12.5，忽略两端效应，按无限长圆柱处理，问题化为沿半径方向的一维轴对称问题。")
    body(doc, "（2）药材内部为各向同性、均匀的连续介质，初始温度与水分浓度沿半径均匀分布（T_{0} = 28°C，C_{0} = 2.55 kg/kg）。")
    body(doc, "（3）水分以有效扩散的形式由内部向表面迁移，毛细流动、汽化与收缩等机理集总于题目给出的经验扩散系数中，不显式刻画相变潜热对温度场的影响。")
    body(doc, "（4）表面边界条件为牛顿冷却与对流传质形式（第三类），对流换热系数与对流传质系数为常数（h = 25 W/(m²·K)，hm = 8×10⁻⁷ m/s）。")
    body(doc, "（5）水分浓度采用干基含水率（kg/kg）定义，与题设一致。")
    body(doc, "（6）预热平衡阶段环境温度与水分浓度按附件 1 分段线性变化，恒温干燥阶段维持预热末值。")
    body(doc, "（7）问题四中，药材收缩为各向同性均匀收缩，收缩速率仅由附件 2 的半径数据确定，沿长度方向一致（仍为一维轴对称）。")

    # ===================== 四、符号说明 =====================
    h1(doc, "四、符号说明")
    add_table(doc, "符号说明", ["符号", "含义", "单位"], [
        ["T", "药材温度", "°C"],
        ["C", "药材水分浓度（干基含水率）", "kg/kg"],
        ["t", "时间", "s"],
        ["r", "到药材中心的径向距离", "m"],
        ["R", "药材半径", "m"],
        ["ρ", "密度", "kg/m³"],
        ["c_{p}", "比热容", "J/(kg·K)"],
        ["k", "导热系数", "W/(m·K)"],
        ["D", "水分扩散系数", "m²/s"],
        ["h", "对流换热系数", "W/(m²·K)"],
        ["h_{m}", "对流传质系数", "m/s"],
        ["T_{∞}、C_{∞}", "烘房环境温度、环境水分浓度", "°C、kg/kg"],
        ["α", "热扩散系数", "m²/s"],
        ["E", "归一化坐标 E = r/R(t)", "—"],
        ["N", "径向网格节点数", "—"],
        ["Δt、Δr", "时间步长、空间步长", "s、m"],
        ["t_{dry}", "烘干时长", "h"],
        ["θ", "时间格式权重（θ = 1/2）", "—"],
    ], first_col_left=True)

    # ===================== 五、模型建立与求解 =====================
    h1(doc, "五、模型建立与求解")

    # ---------- 5.1 问题一 ----------
    h2(doc, "5.1 问题一的模型建立与求解")
    h3(doc, "5.1.1 模型建立")
    body(doc, "（1）控制方程。基于假设（1）（2），预热平衡阶段药材温度与水分浓度满足径向一维轴对称热质耦合方程（0 < r < R，t > 0）：")
    eq(doc, "ρc_{p} ∂T/∂t = (1/r) ∂/∂r (k r ∂T/∂r)", 1)
    eq(doc, "∂C/∂t = (1/r) ∂/∂r (D r ∂C/∂r)", 2)
    body(doc, "（2）定解条件。初值、对称边界与表面第三类对流边界分别为：")
    eq(doc, "T(r,0)=T_{0}=28°C,  C(r,0)=C_{0}=2.55 kg/kg", 3)
    eq(doc, "r=0:  ∂T/∂r=0,  ∂C/∂r=0", 4)
    eq(doc, "r=R:  -k ∂T/∂r = h(T-T_{∞}(t)),  -D ∂C/∂r = h_{m}(C-C_{∞}(t))", 5)
    body(doc, "（3）物性参数（附录 2）。密度 ρ = 820 kg/m³、比热容 c_{p} = 2600 J/(kg·K)、导热系数 k = 0.36 W/(m·K)、对流换热系数 h = 25 W/(m²·K)、对流传质系数 h_{m} = 8×10⁻⁷ m/s，水分扩散系数随水分浓度呈指数衰减：")
    eq(doc, "D(C) = 7×10^{-9} e^{-0.89/C}", 6)
    body(doc, "环境温度 T_{∞}(t) 与环境水分浓度 C_{∞}(t) 由附件 1 分段线性插值给定。")
    body(doc, "（4）无量纲特征分析。热扩散系数 α = k/(ρc_{p}) = 1.69×10⁻⁷ m²/s，初始湿扩散系数 D(C_{0}) = 7×10⁻⁹ e^{−0.89/2.55} = 4.94×10⁻⁹ m²/s，二者之比约 34；热毕渥数 Bi_{T} = hR/k = 1.39，质毕渥数 Bi_{m} = h_{m} R/D ≈ 3.24。这表明温度远快于水分趋于均匀、且表面传质阻力不可忽略，与后续计算结果一致。")
    h3(doc, "5.1.2 模型求解")
    body(doc, "采用守恒型有限体积法离散空间。对控制体 [r_{j−1/2}, r_{j+1/2}] 积分式(1)(2)，得半离散格式：")
    eq(doc, "ρc_{p} V_{j} dT_{j}/dt = [k r ∂T/∂r]_{r_{j+1/2}} - [k r ∂T/∂r]_{r_{j-1/2}}", 7)
    eq(doc, "V_{j} dC_{j}/dt = [D r ∂C/∂r]_{r_{j+1/2}} - [D r ∂C/∂r]_{r_{j-1/2}}", 8)
    eq(doc, "V_{j} = (r_{j+1/2}^{2} - r_{j-1/2}^{2})/2", 9)
    body(doc, "其中 V_{j} 为控制体体积权重，满足恒等式 Σ_{j} V_{j} = R²/2，保证离散格式的严格守恒。界面通量用中心差分，中心节点用洛必达法则处理，表面节点用半控制体积精确平衡（使 Robin 边界进入守恒式）。时间方向采用 θ 加权隐式格式：")
    eq(doc, "(u_{j}^{n+1} - u_{j}^{n})/Δt = θ F^{n+1} + (1-θ) F^{n},  θ = 1/2", 10)
    body(doc, "θ = 1/2 即 Crank–Nicolson 格式，无条件稳定且时间二阶收敛。因 D(C) 非线性，式(2)在每个时间步对扩散系数做滞后线性化并配合牛顿迭代修正（含 ∂D/∂C 敏度项），2~3 次迭代收敛；三对角线性方程组用托马斯追赶法直接求解。")
    body(doc, "计算步骤为：读入附件 1 并分段线性插值得到任意时刻的环境温湿度；初始化 T_{j}^{0} = 28、C_{j}^{0} = 2.55；逐时间步先解热方程得 T^{n+1}，再解质方程（含线性化）得 C^{n+1}，迭代至收敛；输出指定时刻剖面并写出 result1.xlsx。")
    h3(doc, "5.1.3 结果与分析")
    body(doc, "预热平衡阶段（前 30 min）的温度与水分浓度分布见表 1、表 2，剖面演化见图 1。")
    add_table(doc, "表 1  30 分钟内药材的温度（°C）", ["时间/s", "0 cm", "0.5 cm", "1 cm", "1.5 cm", "2 cm"], [
        ["100", "28.0001", "28.0004", "28.0041", "28.0326", "28.1786"],
        ["300", "28.0411", "28.0638", "28.1514", "28.3676", "28.8469"],
        ["600", "28.4537", "28.5363", "28.8039", "29.3153", "30.1640"],
        ["900", "29.3245", "29.4584", "29.8753", "30.6156", "31.7297"],
        ["1200", "30.5428", "30.7098", "31.2222", "32.1121", "33.4269"],
        ["1500", "31.9957", "32.1867", "32.7658", "33.7460", "35.1204"],
        ["1800", "33.5752", "33.7719", "34.3640", "35.3619", "36.7853"],
    ])
    add_table(doc, "表 2  30 分钟内药材的水分浓度（kg/kg，干基）", ["时间/s", "0 cm", "0.5 cm", "1 cm", "1.5 cm", "2 cm"], [
        ["100", "2.5500", "2.5500", "2.5500", "2.5500", "2.2880"],
        ["300", "2.5500", "2.5500", "2.5500", "2.5486", "2.0661"],
        ["600", "2.5500", "2.5500", "2.5500", "2.5342", "1.8846"],
        ["900", "2.5500", "2.5500", "2.5495", "2.5041", "1.7596"],
        ["1200", "2.5500", "2.5500", "2.5478", "2.4647", "1.6619"],
        ["1500", "2.5500", "2.5499", "2.5440", "2.4212", "1.5812"],
        ["1800", "2.5500", "2.5496", "2.5376", "2.3763", "1.5120"],
    ])
    fig(doc, "fig1_预热平衡阶段剖面.pdf", "图 1  预热平衡阶段药材温度与水分浓度剖面")
    body(doc, "由表 1、表 2 可知：1800 s 时中心温度 33.5752°C、表面温度 36.7853°C，温度场在径向已明显趋于均匀；而水分浓度中心 2.5500、表面 1.5120 kg/kg，仅在表面约 1 cm 的薄层内显著下降，中心几乎未变。这一\u201c温度先于水分趋于均匀\u201d的现象正是热扩散系数约为湿扩散系数 34 倍、且表面传质阻力较大（Bi_{m} ≈ 3.24）的体现。")
    h3(doc, "5.1.4 检验分析")
    body(doc, "为验证模型的正确性，进行以下检验：")
    body(doc, "（1）解析解对比：在常环境温度、纯导热情形下，数值解与分离变量级数解析解的最大误差在 1×10⁻⁴ 量级，验证了 Robin 边界与中心奇点处理的正确性（图 2）。")
    body(doc, "（2）收敛阶：固定时间步，空间方向以 N = 20、40、80、160 加密，温度与水分浓度的无穷范数误差均呈二阶收敛（阶约 2.0）；固定空间步，时间方向以 Δt 加倍，误差亦呈二阶收敛，符合 Crank–Nicolson 格式的理论阶（图 3）。")
    body(doc, "（3）守恒性：全局能量与质量平衡的相对偏差分别为 2.9×10⁻¹⁴ 与 1.9×10⁻¹⁴，在机器精度水平上精确守恒。")
    body(doc, "（4）稳定性：Crank–Nicolson 格式在 Δt 放大到 300 s 时仍与 Δt = 1 s 的结果一致（有界收敛），而显式格式在傅里叶数超过 0.5 后立即非物理放大（图 4），印证了隐式格式的无条件稳定性。")
    fig(doc, "fig4_解析解验证.pdf", "图 2  数值解与级数解析解对比（常环境纯导热）")
    fig(doc, "fig5_收敛阶验证.pdf", "图 3  空间与时间方向的收敛阶验证")
    fig(doc, "fig6_离散格式对比.pdf", "图 4  Crank–Nicolson 与显式格式的稳定性对比")
    h3(doc, "5.1.5 小结")
    body(doc, "问题一建立了预热平衡阶段的径向一维热质耦合模型，以守恒型有限体积法配合隐式 Crank–Nicolson 格式与牛顿迭代求解，给出了前 30 min 的温度与水分浓度分布。结果定量反映了\u201c温度远快于水分趋于均匀、水分仅在表面薄层下降\u201d的物理规律，并通过解析解、收敛阶、守恒性与稳定性四类检验验证了模型的正确性。")

    # ---------- 5.2 问题二 ----------
    h2(doc, "5.2 问题二的模型建立与求解")
    h3(doc, "5.2.1 模型建立")
    body(doc, "问题二将时间尺度扩展到整个烘干过程（2~3 天），物性由附录 3 给出，随水分浓度（及温度）变化：")
    eq(doc, "ρ(C)=650+128C,  c_{p}(C)=1450+2736 C/(C+1),  k(C)=0.21+0.38 C/(C+1)", 11)
    eq(doc, "D(C,T)=2.4×10^{-3} e^{-0.45/C} e^{-3850/T_{K}}", 12)
    body(doc, "其中 T_K 为开尔文温度（T_K = T + 273.15）。控制方程仍为式(1)(2)，但密度、比热、导热系数随水分浓度变化使热传导方程也成为非线性方程，热与湿两场通过变物性强耦合。环境条件在预热平衡阶段按附件 1 分段线性变化，进入恒温干燥阶段后维持附件 1 的末值不变。")
    h3(doc, "5.2.2 模型求解")
    body(doc, "在问题一离散格式的基础上，将物性系数在每个时间步冻结于当前状态并对非线性项做牛顿型线性化（扩散系数在状态点处泰勒展开一阶项，含 ∂D/∂C 与 ∂D/∂T 敏度），温度与水分方程弱耦合交替求解，每步 Picard 迭代 2~3 次收敛，三对角方程组用托马斯算法求解。图 5 给出了附录 3 变参数物性随水分浓度的变化。")
    fig(doc, "fig9_变参数物性.pdf", "图 5  附录 3 变参数物性随水分浓度的变化")
    h3(doc, "5.2.3 结果与分析")
    body(doc, "3 h 内的温度与水分浓度分布见表 3、表 4，整个烘干过程的剖面演化见图 6。")
    add_table(doc, "表 3  3 小时内药材的温度（°C）", ["时间/h", "0 cm", "0.5 cm", "1 cm", "1.5 cm", "2 cm"], [
        ["0.5", "32.1894", "32.3822", "32.9660", "33.9608", "35.4137"],
        ["1.0", "40.3818", "40.5538", "41.0609", "41.8790", "42.9985"],
        ["1.5", "45.8470", "45.9351", "46.1934", "46.6056", "47.1412"],
        ["2.0", "48.4504", "48.4883", "48.5987", "48.7739", "49.0038"],
        ["2.5", "49.4671", "49.4793", "49.5138", "49.5655", "49.6604"],
        ["3.0", "49.8496", "49.8554", "49.8747", "49.9102", "49.9664"],
    ])
    add_table(doc, "表 4  3 小时内药材的水分浓度（kg/kg，干基）", ["时间/h", "0 cm", "0.5 cm", "1 cm", "1.5 cm", "2 cm"], [
        ["0.5", "2.5499", "2.5487", "2.5249", "2.3267", "1.6484"],
        ["1.0", "2.5249", "2.4941", "2.3577", "2.0234", "1.4704"],
        ["1.5", "2.3854", "2.3252", "2.1344", "1.8020", "1.3470"],
        ["2.0", "2.1707", "2.1083", "1.9236", "1.6258", "1.2307"],
        ["2.5", "1.9566", "1.9006", "1.7359", "1.4718", "1.1163"],
        ["3.0", "1.7662", "1.7165", "1.5701", "1.3331", "1.0078"],
    ])
    fig(doc, "fig10_全烘干过程剖面.pdf", "图 6  整个烘干过程的温度与水分浓度剖面")
    body(doc, "由表 3、表 4 可知：预热阶段（约前 4 h）药材整体升温，至 3 h 时中心温度升至 49.85°C、表面 49.97°C，径向温差已小于 0.2°C；水分浓度则持续下降，3 h 时中心 1.7662、表面 1.0078 kg/kg，表面下降速度明显快于内部。模拟至 3 天时，中心温度 50.17°C、中心水分浓度 0.1364 kg/kg、表面水分浓度 0.0516 kg/kg，药材已接近干燥终点。")
    h3(doc, "5.2.4 检验分析")
    body(doc, "牛顿型线性化较零阶滞后格式显著加快收敛（每步 2~3 次迭代即达 10⁻¹⁰ 容差）；全局质量守恒在离散意义下精确成立。问题二的格式是问题一的直接推广，问题一的收敛阶、守恒性与稳定性结论同样适用。")
    h3(doc, "5.2.5 小结")
    body(doc, "问题二建立了整个烘干过程的非线性变参数热质耦合模型，通过牛顿型线性化与 Picard 迭代稳定求解强非线性方程组，给出了 3 h 内及 3 天全程的温度与水分浓度演化，为问题三的烘干时长判定提供了完整的时空场。")

    # ---------- 5.3 问题三 ----------
    h2(doc, "5.3 问题三的模型建立与求解")
    h3(doc, "5.3.1 模型建立")
    body(doc, "问题三沿用问题二的模型（附录 3 变参数物性），终止判据为：当药材各处水分浓度均低于 0.15 kg/kg 时达到烘干要求，即 max_{r} C(r,t) < 0.15 kg/kg。由问题二的末期剖面可知，干燥末期表面出现极薄的扩散边界层，中心成为最迟达标位置，故该判据等价于中心水分浓度低于 0.15 kg/kg（后续检验予以证实）。")
    h3(doc, "5.3.2 模型求解")
    body(doc, "（1）分级网格。干燥末期表面扩散边界层厚度约 3×10⁻⁵ m，均匀网格需 Δr 约 1×10⁻⁵ m 才能分辨，导致烘干时长空间收敛缓慢（表观阶仅约 1.0）。为此引入向表面等比加密的分级网格：")
    eq(doc, "r_{j} = R[1-(1-ξ_{j})^{q}],  ξ_{j}=j/N,  q=2", 13)
    body(doc, "q = 2 时用约 60 个节点即可在表面薄层内布置 3 个单元，节点数较均匀网格降低约一个量级而精度相当（同一精度下 CPU 时间减少约 5 倍）。")
    body(doc, "（2）自适应时间步长。采用步长加倍法做误差控制：以 60 s 为基准宏观步长，逐步长内用变子步长积分并估计误差，误差超限则减半、偏低则加倍。预热期瞬变段自动细分至约 2 s 子步，恒温干燥段放大到 300 s，兼顾精度与效率（图 8）。")
    body(doc, "（3）终止时刻二分定位。以 60 s 为输出节拍判断判据，在跨越判据的相邻时刻间二分，精确定位 max_{r} C = 0.15 的终止时刻。")
    h3(doc, "5.3.3 结果与分析")
    body(doc, "烘干时长 t_{dry} = 205771.8 s = 57.16 h（按分钟向上取整的保守口径为 57.17 h，即 57 h 10 min）。每隔 6 h 的水分浓度见表 5，剖面与干燥前沿见图 7、图 9。")
    add_table(doc, "表 5  药材烘干过程的水分浓度（kg/kg，干基）", ["时间/h", "0 cm", "0.5 cm", "1 cm", "1.5 cm", "2 cm"], [
        ["6", "1.0156", "0.9868", "0.9006", "0.7546", "0.5335"],
        ["12", "0.4548", "0.4418", "0.4019", "0.3289", "0.1639"],
        ["18", "0.2979", "0.2906", "0.2679", "0.2247", "0.0842"],
        ["24", "0.2371", "0.2321", "0.2161", "0.1851", "0.0663"],
        ["30", "0.2053", "0.2013", "0.1887", "0.1637", "0.0597"],
        ["36", "0.1854", "0.1821", "0.1714", "0.1500", "0.0565"],
        ["42", "0.1716", "0.1687", "0.1593", "0.1404", "0.0547"],
        ["48", "0.1614", "0.1588", "0.1503", "0.1331", "0.0536"],
        ["54", "0.1535", "0.1511", "0.1433", "0.1274", "0.0528"],
        ["57.1667（烘干结束）", "0.1500", "0.1477", "0.1402", "0.1249", "0.0525"],
    ])
    fig(doc, "fig13_全烘干剖面.pdf", "图 7  烘干全过程的水分浓度剖面与特征位置时程")
    fig(doc, "fig15_干燥前沿.pdf", "图 8  干燥前沿随时间的推进")
    fig(doc, "fig12_自适应步长.pdf", "图 9  自适应时间步长与误差控制")
    body(doc, "由图 8 可知，水分浓度等值线（干燥前沿）由表面向内部逐步推进，中心是最后达标的位置；表 5 中 57.17 h 时中心水分浓度恰降至 0.1500，各处均已低于 0.15，满足烘干要求。")
    h3(doc, "5.3.4 检验分析")
    body(doc, "（1）网格收敛：均匀网格 N 由 20 增至 640 时 t_{dry} 由 56.18 h 单调升至 57.16 h，Richardson 外推极限 57.17 h；分级网格 q = 2 仅用 60 节点即达 57.17 h，与均匀网格极限一致。")
    body(doc, "（2）时间步无关性：改变输出节拍（30~120 s）、子步上限（1~10）与误差容差（1×10⁻⁶ ~ 1×10⁻⁸），t_{dry} 变化均在 ±0.01 min 内。")
    body(doc, "（3）准稳态外推独立校核：自 36 h 起，计入形状漂移的准稳态外推与二分结果一致到 0.34 min 以内，两种误差来源不同的方法互相印证。")
    body(doc, "（4）判据核查：逐节点统计最迟达标时刻，恒为中心节点（r = 0），故\u201c各处低于 0.15\u201d与\u201c中心低于 0.15\u201d等价。")
    body(doc, "（5）质量守恒：单步离散质量恒等式的残差约 1×10⁻⁹ 量级，边界通量积分的两种求积格式相对偏差约 4×10⁻⁵，均在输出精度（5×10⁻⁵）以内。")
    h3(doc, "5.3.5 小结")
    body(doc, "问题三在问题二模型基础上，以分级网格与自适应时间步长解决干燥末期薄边界层分辨与长时间步推进的矛盾，用二分法精确定位终止时刻，得到烘干时长 57.17 h，并通过网格收敛、时间步无关性、准稳态外推、判据核查与质量守恒多重检验确认了结论的稳定性。")

    # ---------- 5.4 问题四 ----------
    h2(doc, "5.4 问题四的模型建立与求解")
    h3(doc, "5.4.1 模型建立")
    body(doc, "（1）移动边界。实际烘干中水分流失使药材收缩，附件 2 给出半径 R(t) 由 2.0 cm 收缩至 1.198 cm（收缩集中在前约 24 h）。固定计算域 [0, R] 不再适用，引入归一化坐标将移动边界映射到固定域：")
    eq(doc, "E = r/R(t)", 14)
    body(doc, "（2）坐标变换（Landau 变换）。由链式法则，固定物理坐标 r 处的时间导数与固定归一化坐标 E 处的时间导数满足：")
    eq(doc, "∂C/∂t|_{r} = ∂C/∂τ|_{E} - (E R′/R) ∂C/∂E", 15)
    body(doc, "代入式(2)得变换后的质扩散方程：")
    eq(doc, "∂C/∂τ = (1/R^{2})(1/E) ∂/∂E (E D ∂C/∂E) + (E R′/R) ∂C/∂E", 16)
    body(doc, "其中对流项 +(E R′/R)·∂C/∂E 由边界收缩引起：收缩使同一 E 处的物理位置内移，水分被\u201c向内输送\u201d，故符号为正（由链式法则严格导出，见 5.4.4 检验）。R′ = dR/dt 为收缩速率。温度方程作同样变换，仅系数为 ρ(C)c_{p}(C)、k(C)。")
    body(doc, "（3）收缩速率。R(t) 由附件 2 分段线性插值，R′(t) 由中心差分数值微分得到：")
    eq(doc, "R′(t_{k}) = [R(t_{k+1}) - R(t_{k-1})]/(2Δt_{R})", 17)
    body(doc, "（4）物性参数（附录 4）。与附录 3 形式相同、系数不同：")
    eq(doc, "ρ(C)=760+90C,  c_{p}(C)=1850+2150 C/(C+1),  k(C)=0.12+0.20 C/(C+1)", 18)
    eq(doc, "D(C,T)=4.2×10^{-4} e^{-0.30/C} e^{-3850/T_{K}}", 19)
    body(doc, "网格在 E 域固定，物理网格随 R(t) 等比例缩放，从而保证整个烘干过程中空间分辨率一致。")
    h3(doc, "5.4.2 模型求解")
    body(doc, "在 E 域 [0,1] 上用守恒型有限体积离散，扩散项与对流项分别中心差分（对流项在表面用单侧差分），守恒权重满足恒等式 Σ_{j} V_{j} = 1/2；时间方向用 Crank–Nicolson 格式，非线性用牛顿线性化加 Picard 迭代，三对角方程组用托马斯算法求解。针对收缩初始段收缩速率较大引起的 Crank–Nicolson 起步振荡，前 1 h 用 Δt = 2 s 子步消除表面浓度的非物理过冲。")
    h3(doc, "5.4.3 结果与分析")
    body(doc, "考虑收缩后的烘干时长 t_{dry} = 188610 s = 52.39 h（按分钟向上取整的保守口径为 52.40 h）。每隔 6 h 的水分浓度见表 6，半径收缩、移动边界剖面与时程见图 10~12。")
    add_table(doc, "表 6  药材烘干过程的水分浓度（kg/kg，干基）", ["时间/h", "0 cm", "0.5 cm", "1.0 cm", "1.5 cm", "药材表面"], [
        ["6.0000", "1.9524", "1.7828", "1.2571", "——", "0.5684"],
        ["12.0000", "0.8705", "0.7739", "0.4848", "——", "0.2145"],
        ["18.0000", "0.4572", "0.4112", "0.2641", "——", "0.1016"],
        ["24.0000", "0.3058", "0.2793", "0.1895", "——", "0.0710"],
        ["30.0000", "0.2370", "0.2189", "0.1549", "——", "0.0608"],
        ["36.0000", "0.1992", "0.1854", "0.1352", "——", "0.0565"],
        ["42.0000", "0.1754", "0.1641", "0.1223", "——", "0.0543"],
        ["48.0000", "0.1591", "0.1494", "0.1132", "——", "0.0531"],
        ["52.3917（烘干结束）", "0.1500", "0.1412", "0.1081", "——", "0.0525"],
    ], note="注：\u201c——\u201d表示该固定距离处材料已因收缩离开（r ≥ R(t)，无药材内部水分）。")
    fig(doc, "fig16_半径收缩与对流.pdf", "图 10  附件 2 收缩半径 R(t) 及其中心差分数值微分")
    fig(doc, "fig17_移动边界剖面.pdf", "图 11  移动边界下的物理坐标剖面与归一化坐标剖面")
    fig(doc, "fig18_特征位置时程.pdf", "图 12  特征位置水分浓度时程与收缩效应对比")
    body(doc, "由图 10 可知，收缩集中在前约 24 h，此后半径稳定在约 1.198 cm；图 11 中物理坐标剖面的右端随收缩内移，归一化坐标下边界恒在 E = 1、剖面因对流项被整体抬升（浓缩效应）。收缩对烘干时长的影响显著：有收缩时 t_{dry} = 52.4 h，而将半径冻结为 2 cm（无收缩）时 t_{dry} = 129.1 h，缩短约 2.46 倍（对应 R² 缩小 2.79 倍、扩散路径缩短与 1/R² 放大共同作用）。")
    h3(doc, "5.4.4 检验分析")
    body(doc, "（1）收敛性：空间方向 N = 100、200、400、800 时 t_{dry} 依次为 52.37、52.38、52.40、52.40 h，N ≥ 400 已收敛；时间方向 Δt = 120、60、30 s 均给出 52.39~52.40 h，结论与步长无关（图 13）。")
    body(doc, "（2）对流项符号物理检验：对流项取 +1（正确）、0（无对流）、−1（错误）时 t_{dry} 分别为 52.40、50.83、48.65 h。正号使 t_{dry} 增大，符合\u201c收缩向内输送水分、浓缩减缓干燥\u201d的物理；负号非物理地加速干燥，印证链式法则导出的正号正确，且对流项贡献约占 t_{dry} 的 3%，不可忽略。")
    body(doc, "（3）退化一致性：当 R′ = 0（无收缩）时，对流算子逐位为零，模型精确退化到固定边界情形。")
    body(doc, "（4）判据等价性与质量守恒：max_{r} C 恒取在中心（r = 0），故判据可用中心值；固定域下前 2 h 离散质量恒等式与解析通量的最大偏差约 1.2×10⁻⁵，守恒在离散意义下成立。")
    fig(doc, "fig19_验证与对比.pdf", "图 13  问题四的空间收敛、对流项符号检验与烘干时长对比")
    h3(doc, "5.4.5 小结")
    body(doc, "问题四通过 Landau 坐标变换将移动边界问题映射为固定域问题，精确处理了收缩引起的对流项，并以中心差分获得收缩速率，采用附录 4 物性求解。结果表明收缩使烘干时长由 129.1 h 缩短至 52.4 h（缩短 2.46 倍），并通过收敛性、对流项符号、退化一致性与守恒性检验验证了移动边界处理的正确性。")

    # ===================== 六、模型评价与推广 =====================
    h1(doc, "六、模型评价与推广")
    h2(doc, "6.1 模型评价")
    body(doc, "模型的优点在于：（1）统一的守恒型有限体积与隐式时间框架，能量与质量在离散意义下精确守恒（相对偏差约 1×10⁻¹⁴），全程无条件稳定；（2）各问题按需递进、复用充分——问题二引入非线性变物性，问题三引入分级网格与自适应步长，问题四引入移动边界坐标变换，层次清晰；（3）验证手段多元，包括解析解对比、空间/时间收敛阶、守恒性、退化一致性、对流项符号物理检验以及准稳态外推与二分两种独立口径的互相印证。")
    body(doc, "模型的不足在于：（1）一维轴对称假设忽略了长度方向不均匀与两端效应，对非均匀环境或短粗药材需三维推广；（2）有效扩散模型将毛细流动、汽化、收缩耦合等机理集总于经验公式，难以刻画机理细节，且附录经验公式本身带有适用范围；（3）收缩律假设各向同性且沿长度均匀，实际收缩可能各向异性并产生裂纹，未建模；（4）环境参数按分段线性插值，未考虑烘房温湿度的反馈与空间不均匀。")
    h2(doc, "6.2 模型推广")
    body(doc, "（1）将一维轴对称模型推广到二维非轴对称或三维（如矩形截面、多层结构药材），以处理端部效应与非均匀环境；（2）引入孔隙度、毛细压力等多孔介质模型，与经验有效扩散公式互相印证，刻画更细的输运机理；（3）结合能耗与品质指标（如色泽、有效成分降解动力学），将烘干时长与温度制度联合优化，反求最优工艺参数；（4）移动边界坐标变换方法可推广到冻干、溶胀、反应边界等其他自由边界问题。")

    # ===================== 参考文献 =====================
    h1(doc, "参考文献")
    refs = [
        "杨世铭, 陶文铨. 传热学[M]. 4版. 北京: 高等教育出版社, 2006.",
        "陶文铨. 数值传热学[M]. 2版. 西安: 西安交通大学出版社, 2001.",
        "Crank J, Nicolson P. A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type[J]. Proceedings of the Cambridge Philosophical Society, 1947, 43(1): 50-67.",
        "李庆扬, 王能超, 易大义. 数值分析[M]. 5版. 北京: 清华大学出版社, 2008.",
        "Thomas L H. Elliptic problems in linear difference equations over a network[R]. New York: Watson Scientific Computing Laboratory, Columbia University, 1949.",
        "Crank J. The Mathematics of Diffusion[M]. 2nd ed. Oxford: Clarendon Press, 1975.",
        "Luikov A V. Heat and Mass Transfer in Capillary-Porous Bodies[M]. Oxford: Pergamon Press, 1966.",
        "Whitaker S. Simultaneous heat, mass, and momentum transfer in porous media: a theory of drying[J]. Advances in Heat Transfer, 1977, 13: 119-203.",
        "潘永康, 王喜忠, 刘相东. 现代干燥技术[M]. 2版. 北京: 化学工业出版社, 2007.",
        "Mujumdar A S. Handbook of Industrial Drying[M]. 4th ed. Boca Raton: CRC Press, 2014.",
        "Dincer I, Zamfirescu C. Drying Phenomena: Theory and Applications[M]. Chichester: John Wiley & Sons, 2015.",
        "Landau H G. Heat conduction in a melting solid[J]. Quarterly of Applied Mathematics, 1950, 8(1): 81-94.",
        "Crank J. Free and Moving Boundary Problems[M]. Oxford: Clarendon Press, 1984.",
        "LeVeque R J. Finite Volume Methods for Hyperbolic Problems[M]. Cambridge: Cambridge University Press, 2002.",
        "刘相东, 于才渊, 周德仁. 常用工业干燥设备及应用[M]. 北京: 化学工业出版社, 2005.",
        "王补宣. 工程传热传质学[M]. 北京: 科学出版社, 1998.",
    ]
    for i, ref in enumerate(refs, 1):
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.0
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run("[%d] %s" % (i, ref))
        fmt(r, size=10.5)

    # ===================== 附录 =====================
    page_break(doc)
    h1(doc, "附  录")
    h2(doc, "附录一  支撑材料目录")
    add_table(doc, "支撑材料目录", ["类型", "文件", "说明"], [
        ["代码", "code/p1_model.py、p1_run.py、p1_verify.py、p1_figures.py", "问题一模型、求解、验证与配图"],
        ["代码", "code/p2_model.py、p2_run.py", "问题二变参数热质耦合模型与求解"],
        ["代码", "code/p3_model.py、p3_run.py、p3_verify.py、p3_grid.py", "问题三分级网格、自适应步长与烘干时长判定"],
        ["代码", "code/p4_model.py、p4_run.py、p4_verify.py、p4_figures.py", "问题四 Landau 变换移动边界模型与求解"],
        ["结果", "results/result1.xlsx ~ result4.xlsx", "四个问题的完整结果（题设模板格式）"],
        ["结果", "results/problem1_tables.txt 等", "各表、各检验的数值记录"],
        ["图件", "figures/fig1 ~ fig19", "论文所引用的全部数据图"],
    ], first_col_left=True)
    h2(doc, "附录二  问题一核心代码（节选）")
    code_block(doc, '''# 三对角方程组求解（Thomas 追赶法）
def thomas(a, b, c, d):
    n = b.size
    psi, phi = np.empty(n), np.empty(n)
    psi[0] = c[0] / b[0]; phi[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * psi[i - 1]
        psi[i] = c[i] / m
        phi[i] = (d[i] - a[i] * phi[i - 1]) / m
    x = np.empty(n); x[-1] = phi[-1]
    for i in range(n - 2, -1, -1):
        x[i] = phi[i] - psi[i] * x[i + 1]
    return x

# 热传导算子（中心 r=0 用 L'Hopital；内点守恒型有限体积；表面半控制体积精确平衡 Robin）
def build_op_heat(dr, N, R=R_CYL, bc=BC_DEFAULT):
    r = np.arange(N + 1) * dr
    r_half = (np.arange(N) + 0.5) * dr
    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)
    b[0] = -4.0 * ALPHA / dr ** 2; c[0] = 4.0 * ALPHA / dr ** 2   # r=0
    for j in range(1, N):
        A = ALPHA * r_half[j - 1] / (r[j] * dr ** 2)
        Cc = ALPHA * r_half[j] / (r[j] * dr ** 2)
        a[j] = A; c[j] = Cc; b[j] = -(A + Cc)
    if bc == "halfcell":
        VN = R * dr / 2.0 - dr ** 2 / 8.0
        a[N] = ALPHA * r_half[N - 1] / (VN * dr)
        b[N] = -a[N] - R * H_HEAT / (RHO * CP * VN)
        c[N] = 0.0; gT = R * H_HEAT / (RHO * CP * VN)
    return a, b, c, gT''')
    h2(doc, "附录三  问题二/三变参数物性（节选）")
    code_block(doc, '''# 附录 3 变参数经验公式（问题 2/3 统一采用）
def rho(C):    return 650.0 + 128.0 * C
def cp(C):     return 1450.0 + 2736.0 * C / (C + 1.0)
def kcond(C):  return 0.21 + 0.38 * C / (C + 1.0)
def diff_coef(C, T):
    return 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850.0 / (T + 273.15))''')
    h2(doc, "附录四  问题四移动边界核心代码（节选）")
    code_block(doc, '''# 对流算子（系数 E R'/R）：内部中心差分、表面单侧差分、中心为 0
def op_conv(g, R, Rp):
    N = g["N"]; dE = g["dE"]; E = g["E"]
    vE = E * (Rp / R)
    a = np.zeros(N + 1); b = np.zeros(N + 1); c = np.zeros(N + 1)
    j = np.arange(1, N)
    a[j] = -vE[j] / (2.0 * dE); c[j] = +vE[j] / (2.0 * dE)
    a[N] = -vE[N] / dE; b[N] = +vE[N] / dE
    return a, b, c

# 传质算子：扩散（非线性 D(C,T)）+ 对流 + 表面流出
def op_mass_diff(g, C, T, R, Rp):
    N = g["N"]
    rl, rr, V, dl, dr = g["rl"], g["rr"], g["V"], g["dl"], g["dr"]
    D = diff_coef(C, T); Dh = 0.5 * (D[:-1] + D[1:])
    invR2 = 1.0 / (R * R); gC = H_MASS / (R * V[N])
    out = np.zeros(N + 1)
    out[0] = invR2 * Dh[0] * rr[0] * (C[1] - C[0]) / (dr[0] * V[0])
    j = np.arange(1, N)
    out[j] = invR2 * (Dh[j-1]*rl[j]*(C[j-1]-C[j])/dl[j]
                      + Dh[j]*rr[j]*(C[j+1]-C[j])/dr[j]) / V[j]
    out[N] = invR2 * Dh[N-1] * rl[N] * (C[N-1]-C[N]) / (dl[N]*V[N]) - gC * C[N]
    ac, bc, cc = op_conv(g, R, Rp)
    out += _apply((ac, bc, cc), C)
    return out''')

    doc.save(OUT)
    print("已生成：", OUT)


if __name__ == "__main__":
    build()
