# 第1问：Euler与CN数值比较

保持一维径向固定圆柱物理模型不变。Euler和CN共用空间通量离散、物性和边界条件。CN含两个后向Euler半步启动，水分方程用Picard迭代。主结果为81节点、0.25秒；同网格算法比较为41节点、0.25秒。

## 目录
- 本目录：求解脚本、LaTeX结果表生成脚本、Artifact Tool Excel导出脚本及检查脚本。
- `02_数据/问题一`：真实附件1、附件2、用户最新指定的result1模板及每秒边界数据。附件2仅检查，不参与第1问。
- `04_结果/问题一`：CN完整result1.xlsx、指定点CSV、未舍入NPZ以及输入和数值检查JSON。
- `05_图片/问题一`：六幅中文图。方法差异不是对真实解的误差。
- `06_论文/问题一`：一个第1问主节的LaTeX源文件、两份自动生成表格及编译PDF，无摘要和问题陈述。

## 复现
Python需要numpy、matplotlib、openpyxl。无需SciPy；追赶法由代码实现。中文图使用Microsoft YaHei，其他系统需替换为可用中文字体，不能屏蔽缺字警告。

从项目根目录运行：

```powershell
python 03_代码/问题一/问题一_数值比较.py
python 03_代码/问题一/生成论文表格.py
python 03_代码/问题一/验证结果.py
```

求解脚本会生成数值、图和`work/表格数据.json`，不直接写xlsx。Excel导出由`导出结果.mjs`使用`@oai/artifact-tool`执行，需要安装此依赖的Codex运行时。当前环境在本目录建立了指向运行时node_modules的本地junction，该依赖目录不纳入版本控制。其他环境需自行配置模块解析，然后运行：

```powershell
node 03_代码/问题一/导出结果.mjs
```

论文需要XeLaTeX和ctex包。从`06_论文/问题一`运行两遍：

```powershell
xelatex -interaction=nonstopmode -halt-on-error 问题一模型建立.tex
xelatex -interaction=nonstopmode -halt-on-error 问题一模型建立.tex
```

时间离散检查采用1、0.5、0.25秒，空间检查采用21、41、81节点，总计六组求解。四位小数是题目格式，不代表所有点的第四位已收敛。尤其初始表面排湿的网格误差仍较明显；81节点结果是较细数值参考，不是解析真值。无实验数据校准，不改变物理模型。

原始赛题及物性由用户提供，result1模板最后指定来源为`D:/正式比赛/附件/附件3/result1.xlsx`，与早先模板哈希一致。论文采用ctexart而非未提供的Mrite format.cls。统计计时是本机一次运行结果，不是通用性能基准。
