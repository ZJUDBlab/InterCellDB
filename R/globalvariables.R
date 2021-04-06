

#' Cluster split character (global).
#' 
#' @description
#' This is the character used to split clusternames in name of interaction pairs.
#'
#' @details
#' When representing a interaction pair, the names of both clusters involved are put together,
#' and kClustersSplit is the character used to split them.
#' Please make sure that kClustersSplit does not appear in any name of clusters, or you can 
#' modify this by hand.
#'
#'
#'
#' @export
#'
kClustersSplit <- "~"



#' Genes split character in one gene pair (global).
#' 
#' @description
#' This is the character used to split gene names in gene pairs, e.g EGF>EGFR.
#'
#' @details
#' When representing a gene pair, the names of both genes participated are put together,
#' and kGenesSplit is the character used to split them.
#' Please make sure that kGenesSplit does not appear in any name of genes, or you can 
#' modify this by hand.
#'
#'
#'
#' @export
#'
kGenesSplit <- "-~-"



#' Predefined expression change selections
#'
#' @description
#' There are 4 status that one gene pairs can have upon expression change, and this package 
#' defines them in these 4 short expressions. 
#'
#'
#'
#' @export
#'
kpred.exprs.change <- c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn")



#' Predefined mode
#'
#' @description
#' This is extracted from databases depicting mode of actions for functional usages.
#' expression is transcriptional regulation.
#'
#'
#' @export
#'
kpred.mode <- c("activation", "inhibition", "binding", "catalysis", "reaction", "expression", "ptmod", "other")  # "phenotype", "other" not exist anymore



#' Predefined action effect
#'
#' @description
#' This is extracted from databases depicting action effect of actions for functional usages.
#'
#'
#'
#' @export
#'
kpred.action.effect <- c("positive", "negative", "unspecified", "undirected")



# [inside usage]
# For utility, use global value to give the definition
kaction.id.mapped <- c(
	"A---B", # #1, undirected, others are directed
	"A-->B", # #2
	"A<--B", # #3
	"A--|B", # #4
	"A|--B", # #5
	"A--oB", # #6
	"Ao--B"  # #7
	# "undefined for (> 6) and all other(< 0)"
)


