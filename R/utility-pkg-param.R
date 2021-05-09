
#' Description of parameter species
#'
#' @description
#' Use to provide the same parameter description for fgenes.remapped.all
#'
#' @param species \pkg{InterCellDB} supports species in either 'human' or 'mouse'.
#' Other species are not supported to be used for analysis. 
#'
InsideParam.species <- function(
	species
) {
	NULL
}



#'
#' Description of parameter fgenes.remapped.all
#'
#' @description
#' Use to provide the same parameter description for fgenes.remapped.all
#'
#' @param fgenes.remapped.all Data.frame. Informations about feature genes(markers), prefer using return value of 
#' \code{\link{DataPrep.RemapClustersMarkers}} and similar functions start with \code{DataPrep.} here.
#'
#'
Inside.DummyFgenes <- function(
	fgenes.remapped.all
) {
	NULL
}

#'
#' Description of parameter genes.ref.db
#'
#' @description
#' Use to provide the same parameter description for genes.ref.db
#'
#' @param genes.ref.db List. It is taxonomy specific, which contains 3 parts:
#' 1. gene.ncbi.db        	NCBI-gene database for specific species
#' 2. gene.synonyms.db    	GeneID-GeneName-Synonym pairs extracted from (1)
#' 3. gene.dup.synonyms.db  Some genes may have same synonyms, and this can cause mistakes 
#' when match target synonym back to its authorized genename.
#' 
#'
Inside.DummyGenesRefDB <- function(
  genes.ref.db
) {
  NULL
}

#'
#' Description of parameter go.ref.db
#'
#' @description
#' Use to provide the same parameter description for go.ref.db
#'
#' @param go.ref.db Data.frame. Reference GO database provided in this package, the list is given as follows:
#' \itemize{
#'   \item For human: \code{go.human.ref.db}
#'   \item For mouse: \code{go.mouse.ref.db}
#' }
#'
#'
Inside.DummyGORefDB <- function(
	go.ref.db
) {
	NULL
}


#'
#' Description of parameter actions.ref.db
#'
#' @description
#' Use to provide the same parameter description for actions.ref.db
#'
#' @param actions.ref.db Data.frame. Reference pairs actions database provided in this package, the list is given as follows:
#' \itemize{
#'   \item For human: \code{actions.human.ref.db}
#'   \item For mouse: \code{actions.mouse.ref.db}
#' }
#'
#'
Inside.DummyActionsRefDB <- function(
	actions.ref.db
) {
	NULL
}