% If you are writing a PhD dissertation, leave the document class unchanged.
% If you are writing a Master's thesis, change the `dissertation' in document class to `thesis'
\documentclass[dissertation, bookmarks]{grizz_thesis}

% To only compile a subset of the document uncomment the next line
%\includeonly{front_matter, 01Intro}

% If you want to add an extra package, add it here
%\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage[shortlabels]{enumitem}
%\usepackage{enumerate}
%\setlist{nosep}
\usepackage{tikz}
\usetikzlibrary{arrows}
\usepackage{amsmath,amsfonts,amssymb,amsthm,epsfig,epstopdf,titling,url,array,comment}
\usepackage{wasysym}
\usepackage{scalerel}
\usepackage{graphicx}
\usepackage{float}
\usepackage{centernot}
%\usepackage[mathscr]{eucal}
%\DeclareMathAlphabet{\mathcal}{OMS}{cmsy}{m}{n}
%\usepackage{bbold}

% Define custom macros here
\newcommand{\indep}{\rotatebox[origin=c]{90}{$\models$}}
\newcommand{\ind}{\nolinebreak{\perp \!\!\! \perp}}
\newcommand{\Hom}{\mathrm{Hom}}
\newcommand{\rank}{\mathrm{rank}}
\newcommand{\Mat}{\mathrm{Mat}}
\newcommand{\tr}{\mathrm{tr}}
\newcommand{\Scan}{\mathrm{Scan}}
\newcommand{\G}{\mathcal{G}}
\newcommand{\Sub}{\mathrm{\hat{S}ub}}
\newcommand{\Mixt}{\mathrm{Mixt}}
\newcommand{\minor}{\mathrm{minor}}
\newcommand{\sgn}{\mathrm{sgn}}
\newcommand{\Sym}{\mathrm{Sym}}
\newcommand{\LCM}{\mathrm{LCM}}
\newcommand{\LT}{\mathrm{LT}}
\newcommand{\jh}{\mathrm{\hat{j}}}
\newcommand{\kh}{\mathrm{\hat{k}}}
\newcommand{\ih}{\mathrm{\hat{i}}}
\newcommand{\ev}{\mathbf{e}}


\theoremstyle{plain}
\newtheorem{thm}{Theorem}[section]
\newtheorem{lem}[thm]{Lemma}
\newtheorem{prop}[thm]{Proposition}
\newtheorem{cor}[thm]{Corollary}

\theoremstyle{definition}
\newtheorem{defn}[thm]{Definition}%[section] %use to number by section
\newtheorem{conj}[thm]{Conjecture}%[section]
\newtheorem{exmp}[thm]{Example}%[section]
\newtheorem{rem}[thm]{Remark}




\begin{document}

\def\thesistitle{GR\"OBNER BASIS FOR THE DOUBLE DETERMINANTAL IDEAL}%Type your title here in all caps
\def\studentname{Josua Illian}%Type your name here
\def\degreename{Doctor of Philosophy in Applied Mathematical Sciences}%Type your degree here
\def\CommitteeChair{Li Li}%Type your committee chair name here
\def\AdvisorOne{Aycil Cesmelioglu}
\def\AdvisorTwo{Charles Ching-an Cheng}%Type your committee member name here
\def\AdvisorThree{Eddie Cheng}%Type your committee member name here
\def\AdvisorFour{Daniel Steffy}%Type your committee member name here



\extraFormating  % Formatting code -- DO NOT EDIT
\makeTheFrontMatter % Formatting code -- DO NOT EDIT
%%% Start main body
%------------------------------------------
% Formatting code -- DO NOT EDIT
%------------------------------------------
\pagenumbering{gobble}
\pagenumbering{arabic}
\newpage
\titleformat
	{\chapter}
	[display]
	{\normalfont\filcenter\singlespacing}
	{\vspace{0.15in}\MakeUppercase{\chaptertitlename~\thechapter}}
	{1pc}
	{#1}
	[\vspace*{2pc}]
\doublespacing
%------------------------------------------


%% Equation spacing - adjust if necessary
\setlength{\abovedisplayskip}{0.5pc}
\setlength{\belowdisplayskip}{0.5pc}
\setlength{\abovedisplayshortskip}{0pc}
\setlength{\belowdisplayshortskip}{0.5pc}
%%------------------------------------------




% Include each separate chapter source file here
\include{01Chapter}
\include{02Chapter}
\include{New03Chapter}
\include{04Chapter}%Uncomment to add another chapter
\include{05Chapter}%Uncomment to add another chapter
%\include{06Chapter}%Uncomment to add another chapter
%\include{07Chapter}%Uncomment to add another chapter
%\include{08Chapter}%Uncomment to add another chapter
%\include{09Chapter}%Uncomment to add another chapter

% The \StartGrizzAppendices command switches the formatting for the 
% remainder of the document -- Leave this alone
\StartGrizzAppendices

% Include each separate appendix source file here
%\include{appendixA}
%\include{appendixB}%Uncomment to add another appendix
%\include{appendixC}%Uncomment to add another appendix
%\include{appendixD}%Uncomment to add another appendix

%This makes the bibliography
\makeTheReferences


\end{document}
