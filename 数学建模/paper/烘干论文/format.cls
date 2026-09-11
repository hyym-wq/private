% format.cls —— CUMCM 竞赛论文格式类
% 用于全国大学生数学建模竞赛论文排版
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{format}[2026/09/11 CUMCM paper class v1.0]

% ============ 基于 ctexart 载入 ============
\LoadClass[12pt,a4paper]{ctexart}

% ============ 基本宏包 ============
\RequirePackage{amsmath,amssymb,amsthm}
\RequirePackage[margin=2.5cm]{geometry}
\RequirePackage{booktabs}
\RequirePackage{graphicx}
\RequirePackage{array}
\RequirePackage{caption}
\RequirePackage{float}
\RequirePackage{enumitem}
\RequirePackage{multirow}
\RequirePackage{xcolor}
\RequirePackage{hyperref}

% ============ 格式设置 ============
\captionsetup{font=small,labelfont=bf}
\setlength{\parindent}{2em}
\allowdisplaybreaks
\graphicspath{{../../figures/}{../figures/}{./figures/}{figures/}}
\hypersetup{
  colorlinks=true,
  linkcolor=black,
  citecolor=black,
  urlcolor=blue
}
\setlist{leftmargin=2em,itemsep=2pt,topsep=4pt}

% ============ 自定义数学命令 ============
\newcommand{\dd}{\mathrm{d}}
\newcommand{\e}{\mathrm{e}}
\newcommand{\E}{\mathbb{E}}
\newcommand{\Var}{\operatorname{Var}}
\newcommand{\Cov}{\operatorname{Cov}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\argmin}{\operatorname*{arg\,min}}
\newcommand{\argmax}{\operatorname*{arg\,max}}
\newcommand{\sgn}{\operatorname{sgn}}
\newcommand{\ind}{\mathbf{1}}
\newcommand{\T}{\mathsf{T}}

\endinput
