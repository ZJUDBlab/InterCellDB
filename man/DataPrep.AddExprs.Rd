% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InterCellDBclass.R
\name{DataPrep.AddExprs}
\alias{DataPrep.AddExprs}
\title{Add Expression Data in InterCell Object}
\usage{
DataPrep.AddExprs(fgenes, exprs.data, force.overwrite = FALSE)
}
\arguments{
\item{fgenes}{feature gene table, stored in slot \code{fgenes} of \code{InterCell} object.}

\item{exprs.data}{normalized count matrix}

\item{force.overwrite}{If column 'Exprs' already exists, set this parameter \code{TRUE} to overwrite that.}
}
\description{
This function adds expression data to InterCell Object. It is mostly internally used.
}
