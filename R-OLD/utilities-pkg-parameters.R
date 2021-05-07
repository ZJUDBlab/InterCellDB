

#' @details
#' This package provides databases and tools for scRNA-seq data analysis.
#' @keywords internal
"_PACKAGE"



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
#' Description of parameter interact.pairs.acted
#'
#' @description
#' Use to provide the same parameter description for interact.pairs.acted
#'
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}, 
#' which is formatted data about target interaction pairs.
#'
#'
Inside.DummyInteractPairsActed <- function(
  interact.pairs.acted
) {
  NULL
}



#'
#' Description of parameter VEinfos
#'
#' @description
#' Use to provide the same parameter description for VEinfos
#'
#' @param VEinfos List. It contains informations about vertices and edges, and is exactly return value of
#' \code{GenerateVEinfos()} or \code{TrimVEinfos()}.
#'
#'
Inside.DummyVEinfos <- function(
	VEinfos
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
#' Description of parameter pairs.ref
#'
#' @description
#' Use to provide the same parameter description for pairs.ref
#'
#' @param pairs.ref Data.frame. Reference interaction database provided in this package, the list is given as follows: 
#' \itemize{
#'	 \item For human: \code{pairs.human.experiments.db}, \code{pairs.human.knowledge.db}, \code{pairs.human.prediction.db}
#'	 \item For mouse: \code{pairs.mouse.experiments.db}, \code{pairs.mouse.knowledge.db}, \code{pairs.mouse.prediction.db}
#' }
#'
#'
Inside.DummyPairsRef <- function(
  pairs.ref
) {
  NULL
}



#'
#' Description of parameter anno.location.ref
#'
#' @description
#' Use to provide the same parameter description for anno.location.ref
#'
#' @param anno.location.ref Data.frame. Annotation of subcellular locations for genes provided in this package, the list is given as follows:
#' \itemize{
#'   \item For human: \code{anno.location.human.ref.db}
#'   \item For mouse: \code{anno.location.mouse.ref.db}
#' }
#'
#'
Inside.DummyAnnoLocationRefDB <- function(
	anno.location.ref
) {
	NULL
}



#'
#' Description of parameter anno.type.ref
#'
#' @description
#' Use to provide the same parameter description for anno.type.ref
#'
#' @param anno.type.ref Data.frame. Annotation of types(molecular functions) for genes provided in this package, the list is given as follows:
#' \itemize{
#'   \item For human: \code{anno.type.human.ref.db}
#'   \item For mouse: \code{anno.type.mouse.ref.db}
#' }
#'
#'
Inside.DummyAnnoTypeRefDB <- function(
	anno.type.ref
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
#' Description of parameter direction.A.to.B
#'
#' @description
#' Use to provide the same parameter description for direction.A.to.B
#'
#' @param direction.A.to.B Logic. this parameter gets 3 options: TRUE, FALSE, NULL. In InterCellDB-defined gene pairs, gene pairs have 2 genes, with 1 named "inter.GeneName.A", 
#' and the other named "inter.GeneName.B". As for interactions between specific cluster group like G1~G2, the gene A is expressed by the former, which is G1, and 
#' gene B is expressed by the latter, which is G2. So, the 3 options mean: With TRUE, all directional gene pairs 
#' that goes from the former cluster to the latter, like A->B or A--|B, will be extracted, but if it is set FALSE, then all directional gene pairs 
#' goes the opposite way, which is something like A<--B or A|--B, will be extracted. If NULL is set, then all directional gene pairs will be reserved.
#'
#'
Inside.DummyDirectionAtoB <- function(
	direction.A.to.B
) {
	NULL
}

