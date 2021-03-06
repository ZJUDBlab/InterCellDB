% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{Tool.SplitToGenDataFrame}
\alias{Tool.SplitToGenDataFrame}
\title{Split Character To Data.frame}
\usage{
Tool.SplitToGenDataFrame(to.splits.string, to.split.by, res.colnames)
}
\arguments{
\item{to.splits.string}{Character. The string to be splited.}

\item{to.split.by}{Character(1). The string will be splited by this parameter, and it will be directly 
passed to parameter \code{split} of function \code{strsplit}.}

\item{res.colnames}{Character(2). It gives 2 column names in the result table.
The splited string gets to have 2 parts of equal length, and 
they will be reconstructed to form a table so that 2 column names are needed.}
}
\description{
This function generates n-column data frame by spliting character by some specific letters or phrases.
}
