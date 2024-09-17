%% This is file `nudtproposal.cls',
% 出现任何问题请联系（zhangxinhang19@foxmail.com）
% 该模板的著作权属于（zhangxinhang19@foxmail.com）
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{gmcmthesis}[2024/09/10 v2.4 update by Zhangxinhang]

% --- Declare Options --- %

% --- terminates all options processing --- %

% --- Loadclass Infomation --- %
\LoadClass[a4paper,12pt]{article}

% --- Required Packages --- %
\RequirePackage[fontset=none]{ctex}
\RequirePackage{fontspec}
\AddToHook{package/xeCJK/after}{\defaultCJKfontfeatures{}}
\RequirePackage{xeCJK}
% 页面布局
\RequirePackage[includeheadfoot]{geometry}
\RequirePackage{pdfpages}
% 标题、目录
\RequirePackage[center,pagestyles]{titlesec}
\RequirePackage{titletoc}
% 图片
\RequirePackage{graphicx}
\RequirePackage{subfig}
\graphicspath{{image/svg/}{image/}{images/}{figure/}{figures/}{resources/}}
% 表格
\RequirePackage{array}
\RequirePackage{tabularx}
%% 长表格
\RequirePackage{longtable}
%% booktabs 提供了\toprule 等命令.
\RequirePackage{booktabs}
\RequirePackage{makecell}
%% multirow 支持在表格中跨行
\RequirePackage{multirow}
%% 调整间隔, 让表格更好看些
\RequirePackage{bigstrut}
%% 在跨行表格中输入定界符
\RequirePackage{bigdelim}
% 数学宏包
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{bm}
% 设置颜色
\RequirePackage{color,xcolor}
\RequirePackage{shadowtext}
% 保护脆落命令
\RequirePackage{cprotect}
% 超链接
\RequirePackage{hyperref}
\hypersetup{
    % xetex,
    % unicode=false, % hyperref 和 xetex 同时使用时不能开启 Unicode 选项.
    pdfstartview=FitH,
    CJKbookmarks=true,
    bookmarksnumbered=true,
    bookmarksopen=true,
    colorlinks, %注释掉此项则交叉引用为彩色边框(将colorlinks和pdfborder同时注释掉)
    pdfborder=001,   %注释掉此项则交叉引用为彩色边框
    allcolors=black,
    breaklinks=true
}
% 智能引用
\RequirePackage{cleveref}
%% crefformat settings
\crefformat{figure}{#2图~#1#3}
\crefrangeformat{figure}{图~#3#1#4\;到{}\;#5#2#6}
\crefmultiformat{figure}{图~#2#1#3}{和~#2#1#3}{，#2#1#3}{和~#2#1#3}
\crefformat{table}{#2表#1#3}
\crefrangeformat{table}{表#3#1#4\;到{}\;#5#2#6}
\crefmultiformat{table}{表~#2#1#3}{和~#2#1#3}{，#2#1#3}{和~#2#1#3}
\crefformat{equation}{#2式~(#1#3)}
\crefrangeformat{equation}{式~(#3#1#4)\;\~{}\;(#5#2#6)}
\crefmultiformat{equation}{式~(#2#1#3)}{ 和~(#2#1#3)}{，(#2#1#3)}{ 和~(#2#1#3)}
\crefformat{algocf}{#2算法~#1#3}
\crefrangeformat{algocf}{算法~#3#1#4\;到{}\;#5#2#6}
\crefmultiformat{algocf}{算法~#2#1#3}{和~#2#1#3}{，#2#1#3}{和~#2#1#3}
\RequirePackage{datetime}
% 首行缩进
\RequirePackage{indentfirst}
% 设置浮动体的标题
\RequirePackage{caption}
\DeclareCaptionLabelFormat{gmcm}{{\zihao{-4}\song#1~\rmfamily #2}}
\DeclareCaptionLabelSeparator{gmcm}{\hspace{1em}}
\DeclareCaptionFont{gmcm}{\zihao{-4}}
\captionsetup{labelformat=gmcm,labelsep=gmcm,font=gmcm}
\captionsetup[table]{position=top,belowskip=0bp \@plus 2bp \@minus 2bp,aboveskip=6bp \@plus 2bp \@minus 2bp}%
\captionsetup[figure]{position=bottom,belowskip=-3bp \@plus 2bp \@minus 2bp,aboveskip=6bp \@plus 2bp \@minus 2bp}%
\captionsetup[subfloat]
{labelformat=simple,font=gmcm,captionskip=6bp,nearskip=6bp,farskip=0bp,topadjust=0bp}
\renewcommand{\thesubfigure}{(\alph{subfigure})}
\renewcommand{\thesubtable}{(\alph{subtable})}
\renewcommand{\footnotesize}{\zihao{5}}
% 定制列表环境
\RequirePackage{enumitem}
% 下划线
\RequirePackage{xeCJKfntef}
% 尺寸计算
\RequirePackage{calc}
\RequirePackage{url}

% ---  Fonts --- %
% ttf字体 (Windows)
\newcommand{\installttf}{
    %% Windows Thesis Fonts
    \setmainfont{Times New Roman}
    \setsansfont{Arial}
    \setmonofont{Courier New}
    %% Using Office Family Fonts
    \setCJKmainfont[BoldFont={STZhongsong}]{SimSun}
    \setCJKsansfont{SimHei} % Hei
    \setCJKmonofont{FangSong} % Fangsong
    %% alias
    \setCJKfamilyfont{song}[AutoFakeBold={2.17}]{SimSun}
    \setCJKfamilyfont{hei}[AutoFakeBold={2.17}]{SimHei}
    \setCJKfamilyfont{fs}[AutoFakeBold={2.17}]{仿宋_GB2312.ttf} % fang-song
    \setCJKfamilyfont{kai}{KaiTi} % Kai
    \setCJKfamilyfont{hwxinwei}[AutoFakeBold={2.17}]{STXINWEI.TTF} % huawenxinwei
    \setCJKfamilyfont{li}[AutoFakeBold={2.17}]{SIMLI.TTF} % lishu
}
\installttf
\newcommand{\cusong}{\bfseries}
\newcommand{\song}{\CJKfamily{song}}     % 宋体
\newcommand{\fs}{\CJKfamily{fs}}         % 仿宋体
\newcommand{\kai}{\CJKfamily{kai}}       % 楷体
\newcommand{\hei}{\CJKfamily{hei}}       % 黑体
\newcommand{\hwxw}{\CJKfamily{hwxinwei}}   % 华文新魏
\newcommand{\li}{\CJKfamily{li}}   % 隶书
\def\songti{\song}
\def\fangsong{\fs}
\def\kaishu{\kai}
\def\heiti{\hei}
\def\huawenxinwei{\hwxw}
\def\lishu{\li}
% 重定义图表标题
\renewcommand\figurename{图}
\renewcommand\tablename{表}

% --- Page Layout --- %
\geometry{top=30.0mm,bottom=17.5mm,left=22.5mm,right=22.5mm,headsep=8mm}
\renewcommand*{\baselinestretch}{1.38}
% % 修改tabular 环境, 设置表格中的行间距为正文行间距.
% \let\gmcm@oldtabular\tabular
% \let\gmcm@endoldtabular\endtabular
% \renewenvironment{tabular}%
% 	{\bgroup%
% 	\renewcommand{\arraystretch}{1.38}%
% 	\gmcm@oldtabular}%
% 	{\gmcm@endoldtabular\egroup}
% 每行缩进两个汉字
\setlength\parindent{2em}

% --- 标题样式 --- % 
\titleformat{\section}{\centering\bfseries\heiti \zihao{4}}{\zhnumber{\thesection}、}{0 pt}{}
\titlespacing*{\section} {0pt}{0ex}{.2ex} % 设置标题前后间距
\titleformat{\subsection}{\zihao{-4} \bfseries}{\thesubsection}{5 pt}{}
\titlespacing*{\subsection} {0pt}{0ex}{.2ex}
\titleformat{\subsubsection}{\zihao{-4}}{\thesubsubsection}{5 pt}{}
\titlespacing*{\subsubsection} {0pt}{0ex}{.2ex}

% --- 参考文献 --- %
\RequirePackage[numbers]{natbib}
\setlength{\bibsep}{0pt plus 0.3ex}
\RequirePackage{etoolbox}
\AtBeginEnvironment{thebibliography}{%
    \phantomsection
    \addcontentsline{toc}{section}{\refname}
}
\newcommand{\upcite}[1]{\textsuperscript{\cite{#1}}} % 上标形式引用


% --- 摘要 --- % 
\newtoks\tenure
\newtoks\title
\renewenvironment{abstract}{%
  \begin{center}
    {\zihao{-2}\bfseries\huawenxinwei \shadowtext{中国研究生创新实践系列大赛}}\\[10pt]
    {\zihao{2}\bfseries\huawenxinwei \shadowtext{\textcolor{black}{“}}\shadowtext{华为杯”第\the\tenure 届中国研究生}}\\[12pt]
    {\zihao{2}\bfseries\huawenxinwei \shadowtext{数学建模竞赛}}
  \end{center}
  \vspace*{16pt}
  \begin{center}
    \begin{tabular}{cc}
      \zihao{-2}\lishu 题\quad 目： &  \zihao{3}\heiti\the\title \\
      \cline{2-2}    
    \end{tabular}%
  \end{center}
  \if@twocolumn
    \section*{\abstractname}%
  \else
    \begin{center}%
    {\zihao{-2}\lishu{摘\quad 要:}\vspace{\z@}}%
    \end{center}%
    \quotation
  \fi}
  {\if@twocolumn\else\endquotation\newpage\null\fi}
\renewenvironment{quotation}
  {\list{}{\listparindent 2em%
   \itemindent \listparindent
   \rightmargin\z@
   \leftmargin\z@
   \parsep \z@ \@plus\p@}%
  \item\relax}
  {\endlist}

% --- 关键词 --- % 
\newcommand\keywords[1]{%
  \newcommand{\gmcm@tokens@keywords}{#1}
  \par
  \vskip1ex
  {\noindent\zihao{-2}\lishu 关键词：}{\zihao{-4}\songti\gmcm@tokens@keywords}
}

% --- 定义附录 --- % 
\RequirePackage[title]{appendix}
\AtBeginEnvironment{appendices}{
  \titleformat{\section}{\centering\bfseries\heiti \zihao{4}}{附录~\thesection}{5 pt}{}
\titlespacing*{\section} {0pt}{0ex}{.2ex} % 设置标题前后间距
}
% --- 数学环境, 定理等设置 --- %
\newcommand*{\gmcm@cap@definition}{定义}
\newcommand*{\gmcm@cap@theorem}{定理}
\newcommand*{\gmcm@cap@lemma}{引理}
\newcommand*{\gmcm@cap@corollary}{推论}
\newcommand*{\gmcm@cap@assumption}{假设}
\newcommand*{\gmcm@cap@conjecture}{猜想}
\newcommand*{\gmcm@cap@axiom}{公理}
\newcommand*{\gmcm@cap@principle}{定律}
\newcommand*{\gmcm@cap@problem}{问题}
\newcommand*{\gmcm@cap@example}{例}
\newcommand*{\gmcm@cap@proof}{证明}
\newcommand*{\gmcm@cap@solution}{解}
\newtheorem{definition}{\gmcm@cap@definition}
\newtheorem{theorem}{\gmcm@cap@theorem}
\newtheorem{lemma}{\gmcm@cap@lemma}
\newtheorem{corollary}{\gmcm@cap@corollary}
\newtheorem{assumption}{\gmcm@cap@assumption}
\newtheorem{conjecture}{\gmcm@cap@conjecture}
\newtheorem{axiom}{\gmcm@cap@axiom}
\newtheorem{principle}{\gmcm@cap@principle}
\newtheorem{problem}{\gmcm@cap@problem}
\newtheorem{example}{\gmcm@cap@example}
\newtheorem{proof}{\gmcm@cap@proof}
\newtheorem{solution}{\gmcm@cap@solution}

% --- 浮动环境设置 --- %
%% 默认情况下, \LaTeX{} 要求每页的文字至少占据 20%,否则该页就只单独放置一个浮动环境, 而这通常不是我们想要的, 我们将这个要求降低到 5%.
\renewcommand*{\textfraction}{0.05}
%% 有时如果多个浮动环境连续放在一起, \LaTeX{}会将它们分在几个不同页，即使它们可在同一页放得下. 我们可以通过修改 |\topfraction| 和 |\bottomfraction| 分别设置顶端和底端的浮% 动环境的最大比例.
\renewcommand*{\topfraction}{0.9}
\renewcommand*{\bottomfraction}{0.8}
%% 有时\LaTeX{}会把一个浮动环境单独放在一页, 我们要求这个环境至少要占据 85% 才能单独放在一页. 注意:  |\floatpagefraction| 的数值必须小于 |\topfraction|.
\renewcommand*{\floatpagefraction}{0.85}

% --- tokens definition --- %
\newtoks\schoolname
\newtoks\teamid
\newtoks\membera
\newtoks\memberb
\newtoks\memberc
%  --- Make Titlepage --- %
\renewcommand\maketitle{%
  \pagenumbering{Roman}
  \thispagestyle{empty}
  \begin{center}
    \includegraphics{figures/logo.pdf}
  \end{center}

  \vspace*{25pt}

  \begin{center}
    {\zihao{-2}\bfseries\huawenxinwei \shadowtext{中国研究生创新实践系列大赛}}\\[10pt]
    {\zihao{2}\bfseries\huawenxinwei \shadowtext{\textcolor{black}{“}}\shadowtext{华为杯”第\the\tenure 届中国研究生}}\\[12pt]
    {\zihao{2}\bfseries\huawenxinwei \shadowtext{数学建模竞赛}}
  \end{center}

  \vspace*{56pt}
  \begin{center}
    \renewcommand\arraystretch{1.6}
    \begin{tabular}{c<{\raggedleft}p{11.48cm}}
      \multicolumn{1}{c}{\zihao{-2}\bfseries\songti 学\hfill{}校} & \multicolumn{1}{c}{\zihao{-2}\bfseries\songti \the\schoolname} \\
      \hline
          &  \\
      \multicolumn{1}{c}{\zihao{-2}\bfseries\songti 参赛队号} & \multicolumn{1}{c}{\zihao{-2}\bfseries\songti \the\teamid} \\
      \hline
          & \zihao{-2}\bfseries\songti 1. \the\membera \\
      \cline{2-2}    \multicolumn{1}{l}{\zihao{-2}\bfseries\songti 队员姓名} & \zihao{-2}\bfseries\songti 2. \the\memberb \\
      \cline{2-2}        & \zihao{-2}\bfseries\songti 3. \the\memberc \\[0.1cm]
      \hline
    \end{tabular}%
  \end{center}
  \clearpage
  \pagenumbering{arabic}
  \zihao{-4}
  \renewcommand{\normalsize}{\zihao{-4}} %
  \renewcommand{\baselinestretch}{1.35}
}
\AtBeginEnvironment{table}{
  \renewcommand\arraystretch{1.2}
}
\AtBeginEnvironment{longtable}{
  \renewcommand\arraystretch{1.2}
}
% --- pdf元数据 --- %
\hypersetup{
  pdftitle={偷偷把自己的信息写入pdf源数据},
  pdfauthor={张鑫航},
  pdfsubject={研究生数模竞赛},
  pdfkeywords={数模，LaTeX模板}
}
\endinput