% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TgView-part1.R
\name{GetResultPieActionEffect}
\alias{GetResultPieActionEffect}
\title{Get Result for Composition of Action Effect}
\usage{
GetResultPieActionEffect(
  object,
  limits.exprs.change = kexprs.change,
  limits.ext.action.effect = kpred.ext.action.effect,
  color.ext.action.effect = kpred.color.ext.effect,
  legend.title = "Action Effect",
  legend.title.size = element_text(size = 14),
  legend.key.size = unit(7, "mm"),
  legend.text.size = element_text(size = 12),
  legend.box.margin = margin(0, 0, 0, 6),
  show.note = TRUE,
  note.label.size = 3,
  note.lineheight = 1,
  note.postion.xy = c(0.5, 0.5),
  note.hjust = 0.5,
  note.vjust = 0.5
)
}
\arguments{
\item{object}{A \code{InterCell} object, which is created by \code{\link{CreateInterCellObject}}.}

\item{limits.exprs.change}{Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'. It selects 
the part of result to be shown.}

\item{limits.ext.action.effect}{Options are listed in \code{kpred.ext.action.effect}. It specifies 
the range of shown extended action effect.}

\item{color.ext.action.effect}{The colors to be used. It is one-to-one matching with \code{limits.ext.action.effect}.}

\item{legend.title}{The content of legend title.}

\item{legend.title.size}{The size of legend title.}

\item{legend.key.size}{The size of keys in legend.}

\item{legend.text.size}{The size of legend text.}

\item{legend.box.margin}{The margin of legend box, and should be given in \code{margin()}. See \code{?margin} for further help.}

\item{show.note}{If set TRUE, the note will be shown in the graph.}

\item{note.label.size}{Numeric. The size of note label.}

\item{note.lineheight}{Numeric. The lineheight of note.}

\item{note.postion.xy}{The position of note, should be 2 numbers that give the coordinate. Allowed range is (0~1, 0~1).}

\item{note.hjust}{The horizontal alignment of note. Allowed range is 0~1.}

\item{note.vjust}{The vertical alignment of note. Allowed range is 0~1.}
}
\value{
A list. Use \code{Tool.ShowGraph()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
\itemize{
  \item plot: NULL.
  \item grid.plot: the graph.
  \item table: a list of \code{data.frame}.
}
}
\description{
This function focuses on one interaction pair and its detailed information(expression changes, 
gene-gene action effect), and get result from it.
}
