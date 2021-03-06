% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InterCellDBclass.R
\name{MergeGeneOI}
\alias{MergeGeneOI}
\title{Merge Gene of Interest}
\usage{
MergeGeneOI(gene.oi.1, gene.oi.2, merge.option = "intersect")
}
\arguments{
\item{gene.oi.1}{Genes of interest, generated by \code{FetchGeneOI}.}

\item{gene.oi.2}{Genes of interest, generated by \code{FetchGeneOI}.}

\item{merge.option}{The supported options are 'intersect' and 'union'.}
}
\value{
Character, the selected genes. Or, a list, with gene property attached.
}
\description{
This function is to merge genes of interest fetched by \code{FetchGeneOI}.
}
