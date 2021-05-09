
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pre-defined default global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Predefined data format in given fgenes
#'
#' @description
#' Given data on differentially expressed genes and their belonging clusters should be in the right format.
#' Some columns must be maintained, which will be further used in downstream analysis.
#'
#'
kmusthave.colnames <- c("cluster", "gene", "LogFC", "PVal")  # NOT CHANGE

#' Predefined mode & corresponding colour usage for plotting
#'
#' @description
#' This is extracted from databases depicting mode of actions for functional usages.
#' expression is transcriptional regulation.
#'
kpred.mode <- c("activation", "inhibition", "binding", "catalysis", "reaction", "expression", "ptmod", "other")
kpred.color.mode <- c("#FB8072", "#B3DE69", "#80B1D3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDB462", "#FCCDE5")

#' Predefined action effect
#'
#' @description
#' This is extracted from databases depicting action effect of actions for functional usages.
#'
kpred.action.effect <- c("positive", "negative", "unspecified", "undirected")
kpred.color.effect <- c("#FB8072", "#B3DE69", "#80B1D3", "#8DD3C7")

#' Predefined extended action effect
#'
#' @description
#' This depicts action effect and direction of action together.
#'
kpred.ext.action.effect <- c(
	"A---B", # #1, undirected, others are directed
	"A-->B", # #2
	"A<--B", # #3
	"A--|B", # #4
	"A|--B", # #5
	"A--oB", # #6
	"Ao--B"  # #7
	# "undefined for (> 7) and all other(< 0)"
)
kpred.color.ext.effect <- c("#8DD3C7", "#FB8072", "#FF5740", "#B3DE69", "#81EF48", "#80B1D3", "#6B6AEA")


#' Formula Used in Network Analysis
#'
#' @description
#' This is the default formula for calculating overall strength for interactions.
#'
#' @param pairs.ref data.frame. Data recording all gene pairs as well as their attributes in one interaction.
#' @param colname.eval character. The column names to be used to evaluate strength.
#'
FullView.formula.Evaluation.default <- function(
	pairs.ref,
	colname.eval = c("inter.LogFC.A", "inter.LogFC.B")
) {
	eval.res <- 0.0
	if (nrow(pairs.ref) > 0) {
		eval.res <- sum(abs(pairs.ref[, colname.eval[1]] * pairs.ref[, colname.eval[2]]))
	}
	eval.res
}



#' Formula on Power for One Gene Pair
#'
#' @description
#' This is the default formula on calculating power for every gene pair, and the 
#' sum value of LogFC for gene partners will be regarded as power.
#'
#' @param data.f vector. The LogFC values for one list of gene partners. LogFC values are 
#' specified when creating \code{\link{InterCell}} object. 
#' @param data.b vector. The LogFC values for another list of gene partners, which is 
#' one-by-one matched to those in parameter \code{data.f}.
#'
TgView.formula.onLogFC.default <- function(
	data.f, 
	data.b
) {
	if (length(data.f) != length(data.b)) {
		stop("Unexpected non-identical length data input!")
	}
	return(data.f + data.b)
}



#' Formula on Confidence for One Gene Pair
#'
#' @description
#' This is the default formula on calculating confidence for every gene pair, and the 
#' multiplying value of PVal for gene partners will be regarded as confidence. 
#'
#' @param data.f vector. The PVal values for one list of gene partners. PVal values are 
#' specified when creating \code{\link{InterCell}} object. 
#' @param data.b vector. The PVal values for another list of gene partners, which is 
#' one-by-one matched to those in parameter \code{data.f}.
#'
TgView.formula.onPVal.default <- function(
	data.f, 
	data.b,
	pval.log.max = 1000
) {
	# 'pval.log.max' is 1000 in default setting.  
	# The use of 1000 as default maximum, because e-999 are usual lowest limit when all log(values) are infinite.
	default.max.replace <- pval.log.max
	if (length(data.f) != length(data.b)) {
		stop("Unexpected non-identical length data input!")
	}
	tmp.f <- abs(log10(data.f))
	tmp.b <- abs(log10(data.b))
	inds.tmp.f <- which(is.finite(tmp.f))
	inds.tmp.b <- which(is.finite(tmp.b))
	if (length(inds.tmp.f) == 0 && length(inds.tmp.b) == 0) {
		tmp.f <- rep(default.max.replace, times = length(tmp.f))
		tmp.b <- rep(default.max.replace, times = length(tmp.b))
	} else {
		max.f <- suppressWarnings(max(tmp.f[inds.tmp.f]))  # if length(inds.tmp.f) == 0, max returns -Inf
		max.f	<- ifelse(is.infinite(max.f), default.max.replace, max.f)
		max.b <- suppressWarnings(max(tmp.b[inds.tmp.b]))  # if length(inds.tmp.b) == 0, max returns -Inf
		max.b	<- ifelse(is.infinite(max.b), default.max.replace, max.b)
		tmp.f[which(is.infinite(tmp.f))] <- 10 * max.f
		tmp.b[which(is.infinite(tmp.b))] <- 10 * max.b
	}
	return(tmp.f * tmp.b)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# InterCell Object definitions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The Class for Database Storage in InterCellDB
#'
#' The InterCellDBPack class provides a standard storage format of all databases for running analysis. 
#' This object is mostly created along with \code{\link{InterCell-class}}. 
#' Users need to specify the species firstly. Then, user-defined subset of database could be 
#' collected by function \code{\link{SelectDBSubset}}.
#'
#' @slot species InterCellDB only allows 'human' and 'mouse' species.
#' @slot genes.db The gene reference database, which is generated from NCBI gene site. It records all the authorized gene names.
#' @slot pairs.db The database records all gene pairs, which are transferred from protein interaction database - 'STRING'.
#' @slot actions.db The database records all action properties (including 'action.mode', 'action.effect') for gene pairs in slot:`pairs.db`.
#' @slot anno.location.db The database collects all subcellular location from COMPARTMENTS for genes in slot:`pairs.db`.
#' @slot anno.type.db The database collects all molecular function from Uniprot for genes in slot:`pairs.db`.
#' @slot go.ref.db The database collects all GO terms for genes in slot:`pairs.db`.
#' @slot accessory.db A list of accessory database, for example, the mapping list of raw keywords from Uniport. 
#' @slot misc The preserved area for miscellaneous things. Currently no usage.
#'
#' @name InterCellDBPack-class
#' @rdname InterCellDBPack-class
#' @exportClass InterCellDBPack
#'
InterCellDBPack <- setClass(
	Class = "InterCellDBPack",
	slots = c(
		species = "character",
		genes.db = "list", 
		pairs.db = "data.frame",
		actions.db = "data.frame", 
		anno.location.db = "data.frame",
		anno.type.db = "data.frame",
		go.ref.db = "data.frame",
		accessory.db = "list",
		misc = "list"
	)
)



#' The InterCell Class
#'
#' The InterCell object is used to in all the analysis \pkg{InterCellDB} provided.
#' It incoporates the necessary databases, and stores every important intermediate result.
#'
#' @slot fgenes InterCellDB needs the differential expressed genes and their belonging clusters as input. 
#' The given data need to have columns named 'gene', 'cluster', 'LogFC', 'PVal'. Those each means, 'gene': the 
#' authorized gene names, 'cluster': cell cluster, 'LogFC': the log transformed fold change of genes, 'PVal': the 
#' P val for regarding one gene as differentially expressed one. 
#' @slot database Database stored in \code{InterCellDBPack} object. 
#' @slot formulae The default formulae to be used in analysis. 
#' @slot inter.fullview This is used to store the result of network analysis. 
#' @slot tg.action.pairs This is used to store the intermediate result of action properties and gene pairs.  
#' @slot tg.veinfo This is used to store the intermediate result of one intercellular communication (for one interacting 2-cell group).
#' @slot tg.action.comp This is used to store the result of analyzing composition of action mode and action effect in one interacellular communication.
#' @slot tg.spgenes This is used to store the selected gene pairs and all their participating interacitons in one intercellular communication.
#' @slot tool.vars This stores several variables to be used embedded in program. 
#' @slot misc The not important intermediate result will be put in this parameter as well as some 
#' pre-defined settings. 
#'
#' @name InterCell-class
#' @rdname InterCell-class
#' @exportClass InterCell
#'
InterCell <- setClass(
	Class = "InterCell",
	slots = c(
		fgenes = "data.frame", 
		database = "InterCellDBPack", 
		formulae = "list", 
		inter.fullview = "list", 
		tg.action.pairs = "list", 
		tg.veinfo = "list", 
		tg.action.comp = "list",
		tg.spgenes = "list", 
		#pred.action = "list", pred.action This stores the pre-defined action modes and action effects, and terms beyond definitions in this parameter will be regarded as invalid ones.
		tool.vars = "list", 
		misc = "list" 
	)
)



# validation function for \code{InterCellDBPack-class} 
validInterCellDBPackObject <- function(object) {
	# check species
	allowed.species <- c("human", "mouse")
	if (length(object@species) != 1 || sum(allowed.species %in% object@species) == 0) {
		return(paste("InterCellDB only supports species in human and mouse."))
	}
	TRUE
}
setValidity("InterCellDBPack", validInterCellDBPackObject)

# validation function for \code{InterCell-class}
# only the most un-modifiable slots should be checked, 
# and all checked ones should set its initial value in the first time, or get errors
validInterCellObject <- function(object) {
	# check pred.action
	# if (length(object@pred.action$action.mode) == 0 || !all(object@pred.action$action.mode %in% kpred.mode)) {
	# 	return(paste0("Using undefined action mode: ", paste0(setdiff(object@pred.action$action.mode, kpred.mode), collapse = ", ")))
	# }
	# if (length(object@pred.action$action.mode) != length(object@pred.action$color.mode)) {
	# 	return(paste0("Using different length of 'action.mode' and its used color 'color.mode'!"))
	# }
	# if (length(object@pred.action$action.effect) == 0 || !all(object@pred.action$action.effect %in% kpred.action.effect)) {
	# 	return(paste0("Using undefined action effect: ", paste0(setdiff(object@pred.action$action.effect, kpred.action.effect), collapse = ", ")))
	# }
	# if (length(object@pred.action$action.effect) != length(object@pred.action$color.effect)) {
	# 	return(paste0("Using different length of 'action.effect' and its used color 'color.effect'!"))
	# }

	# check k***Split
	if (length(object@tool.vars$gene.pair.split) != 1) {
		return(paste0("Invalid symbol for splitting genes in gene pairs!"))
	}
	if (length(object@tool.vars$cluster.split) != 1) {
		return(paste0("Invalid symbol for splitting clusters in cluster groups!"))
	}

	# check musthave colnames in given fgenes
	if (!all(object@misc$musthave.colnames %in% kmusthave.colnames)) {
		return(paste0("Detect manual modification on object internal variables. ", 
			"Program will be failed in unexpected situations. "))
	}

	validInterCellDBPackObject(object@database)
	TRUE
}
setValidity("InterCell", validInterCellObject)



# initialize \code{InterCellDBPack-class}
setMethod(
	f = "initialize",
	signature = c("InterCellDBPack"),
	definition = function(.Object, ...) {
		.Object <- callNextMethod(.Object, ...)
		.Object@accessory.db <- list(Uniprot.key.map = Uniprot.key.map.list)
		.Object@misc <- list(TAKEN = "nothing yet")
		validObject(.Object)
		return(.Object)
	}
)

# initialize \code{InterCell-class}
setMethod(
	f = "initialize",
	signature = c("InterCell"),
	definition = function(.Object, ...) {
		.Object <- callNextMethod(.Object, ...)
		.Object@formulae <- list(FULLVIEW = FullView.formula.Evaluation.default,
			TG.LOGFC = TgView.formula.onLogFC.default,
			TG.PVAL = TgView.formula.onPVal.default)
		#.Object@pred.action <- list(action.mode = kpred.mode, action.effect = kpred.action.effect)
		#.Object@tool.vars <- list(gene.pair.split = "-~-", cluster.split = "~")
		validObject(.Object)
		return(.Object)
	}
)



setMethod(
	f = "show",
	signature = "InterCellDBPack",
	definition = function(object) {
		cat("Using '", object@species, "' database.\n", sep = "")
	}
)

setMethod(
	f = "show",
	signature = "InterCell",
	definition = function(object) {
		cat("A InterCell object, with ", nrow(object@fgenes), " differentially expressed genes spanning ", 
			length(unique(object@fgenes$cluster)), " clusters.\n", sep = "")
		show(object@database)
		if (!is.null(object@tg.veinfo) && length(object@tg.veinfo) > 0) {
			this.involved.clusters <- getOrigClusterNameTgVEInfo(object)
			cat("Intercellular analysis performed between ", this.involved.clusters$cluster.name.A, " and ", this.involved.clusters$cluster.name.B)
		} else {
			cat("Intercellular analysis not processed yet.")
		}
	}
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Creation Method
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create InterCellDBPack Object
#'
#' @description
#' This function is used to create \code{InterCellDBPack-class}.
#'
#' @inheritParams InsideParam.species
#'
#' @return A InterCellDBPack object.
#'
#' @export
#'
CreateDBPackObject <- function(
	species
) {
	# Check species, only support human and mouse yet
	species <- CheckSpeciesValidity(species)
	# set initial database
	ret.obj <- NULL
	if (species == "human") {
		ret.obj <- new(
			Class = "InterCellDBPack",
			species = species, 
			genes.db = genes.human.ref.db,
			pairs.db = pairs.human.db,
			actions.db = actions.human.ref.db,
			anno.location.db = anno.location.human.ref.db,
			anno.type.db = anno.type.human.ref.db,
			go.ref.db = go.human.ref.db
		)
	}
	if (species == "mouse") {
		ret.obj <- new(
			Class = "InterCellDBPack",
			species = species, 
			genes.db = genes.mouse.ref.db,
			pairs.db = pairs.mouse.db,
			actions.db = actions.mouse.ref.db,
			anno.location.db = anno.location.mouse.ref.db,
			anno.type.db = anno.type.mouse.ref.db,
			go.ref.db = go.mouse.ref.db
		)
	}

	# get database slimmed to useful data in analyzing only
	genes.cols.reserve <- c("tax_id", "GeneID", "Symbol", "Synonyms", "Symbol_from_nomenclature_authority")
	ret.obj@genes.db$gene.ncbi.db <- ret.obj@genes.db$gene.ncbi.db[, genes.cols.reserve]
	anno.locs.cols.reserve <- c("GeneID", "Gene.name", "GO.ID.target", "GO.Term.target", "Source", "Evidence", "score")
	ret.obj@anno.location.db <- ret.obj@anno.location.db[, anno.locs.cols.reserve]
	anno.type.cols.reserve <- c("GeneID", "Gene.name", "Keyword.ID", "Keyword.Name")
	ret.obj@anno.type.db <- ret.obj@anno.type.db[, anno.type.cols.reserve]

	return(ret.obj)
}


#' Analyze interaction network in full view
#'
#' @description
#' This function analyzes count and power of interaction pairs among all given clusters.
#'
#' @param DEG.table A table on differentially expressed genes and their belonging cell clusters.
#' @inheritParams InsideParam.species
#' @param cluster.split The letters used to split 2 cell clusters in one interaction. It can also be modified later by 
#' using \code{setClusterSplit} if not decided yet.
#' @param gene.pair.split The letters used to split 2 gene partners in one gene pair. It can also be modified later by 
#' using \code{setGenePairSplit} if not decided yet.
#'
#' @details
#' The parameter \code{DEG.table} is recommended to be generated by \pkg{Seurat}. Other packages are also applicable, if
#' they work with scRNA-seq, cell clustering and calculation on cluster-specific differentially expressed genes. The input 
#' format of \code{DEG.table} should be one data.frame with 4 required columns that are named 'cluster', 'gene', 'LogFC', 'PVal'.
#' \itemize{
#'   \item{cluster}: the cell cluster.
#'   \item{gene}: differentially expressed genes, which are grouped by their belonging clusters.
#'   \item{LogFC}: the fold change of genes.
#'   \item{PVal}: the P value for gene being calculated as differentially expressed gene.
#' }
#'
#' To represent one interaction, like interaction between 'Myeloid cell' and 'T cell', 
#' it will be looked like 'Myeloid cell~T cell' if \code{cluster.split = "~"}.
#' To represent oen gene pair, like gene pair of IL6 and IL6R,
#' it will be looked like 'IL6-~-IL6R' if \code{gene.pair.split = "-~-"}.
#' 
#' To avoid program failure, the letters appearing in gene names are not recommended for \code{gene.pair.split},
#' and the letters appearing in cluster names are not recommended for \code{cluster.split}. 
#' The program will test for those situation but users should keep this in mind.
#'
#' @return A InterCell object. 
#'
#' @export
#'
CreateInterCellObject <- function(
	DEG.table,
	species,
	cluster.split = "~",
	gene.pair.split = "-~-"
) {
	# pre-check part
	if (!all(kmusthave.colnames %in% colnames(DEG.table))) {
		stop("Required columns are not given!\n",
			paste("Column named ", paste0(setdiff(kmusthave.colnames, colnames(DEG.table)), collapse = ", "), 
			" are not included in given feature genes table.",
			"Please use colnames(<var>)[<index>] <- '<name>' to set proper columns corresponding to those.",
			"Contents represented by <> should be replaced by user definitions!")
		)
	}

	# set object
	DBPack.obj <- CreateDBPackObject(species)
	# re-align input table with gene reference database in InterCellDB
	DEG.align.res <- DataPrep.RemapClustersMarkers(DEG.table, species)
	IT.InterCell.Obj <- new(
		Class = "InterCell",
		fgenes = DEG.align.res$result,
		database = DBPack.obj,
		#pred.action = list(action.mode = kpred.mode, action.effect = kpred.action.effect,
		#	color.mode = kpred.color.mode, color.effect = kpred.color.effect),
		tool.vars = list(gene.pair.split = gene.pair.split, cluster.split = cluster.split), 
		misc = list(musthave.colnames = kmusthave.colnames)
	)
	# give the result of alignment in the misc
	IT.InterCell.Obj@misc$input.align.result <- DEG.align.res
	# set object
	IT.InterCell.Obj <- setGenePairSplit(IT.InterCell.Obj, gene.pair.split)
	IT.InterCell.Obj <- setClusterSplit(IT.InterCell.Obj, cluster.split)

	return(IT.InterCell.Obj)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Accessor Function
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Using FullView Result
#' 
#' The \code{setFullViewResult} function is to \bold{set} the result of network analysis.
#'
#' @rdname FullViewResult-InterCell
#' @export
#'
setGeneric(name = "setFullViewResult", def = function(object, ...) {
	standardGeneric("setFullViewResult")
	}
)
#' @param new.inter.fullview A new result of network analysis.
#'
#' @examples
#' \dontrun{
#'   setFullViewResult(object, some.new.fullview.result)
#' }
#'
#' @rdname FullViewResult-InterCell
#'
setMethod(
	f = "setFullViewResult",
	signature = "InterCell",
	definition = function(object, new.inter.fullview) {
		if (!is.null(object@inter.fullview) && length(object@inter.fullview) != 0) {
			warning("Overwrite existed result on network analysis. ")
		}
		object@inter.fullview <- new.inter.fullview
		return(object)
	}
)

#' Using FullView Result
#'
#' The \code{getFullViewResult} function is to \bold{get} the result of network analysis.
#'
#' @rdname FullViewResult-InterCell
#' @export
#'
setGeneric(name = "getFullViewResult", def = function(object, ...) {
	standardGeneric("getFullViewResult")
	}
)

#' @examples
#' \dontrun{
#'   getFullViewResult(object)
#' }
#' @rdname FullViewResult-InterCell
#'
setMethod(
	f = "getFullViewResult",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@inter.fullview) || length(object@inter.fullview) == 0) {
			stop("Network analysis is not performed yet. ")
		}
		return(object@inter.fullview)
	}
)

#' Using TargetView Action Pairs
#' 
#' The \code{setTgActionPairs} function is to \bold{set} the result of intercellular analysis 
#' on action properties for one interaction between 2 cell clusters.
#'
#' @rdname TgActionPairs-InterCell
#' @export
#'
setGeneric(name = "setTgActionPairs", def = function(object, ...) {
	standardGeneric("setTgActionPairs")
	}
)

#' @param new.action.pairs A new result of intercellular analysis on action properties
#'
#' @examples
#' \dontrun{
#'   setTgActionPairs(object, some.new.tg.action.pairs)
#' }
#'
#' @rdname FullViewResult-InterCell
#'
setMethod(
	f = "setTgActionPairs",
	signature = "InterCell",
	definition = function(object, new.action.pairs) {
		object@tg.action.pairs <- new.action.pairs
		return(object)
	}
)

#' Using TargetView Action Pairs
#' 
#' The \code{getTgActionPairs} function is to \bold{get} the result of intercellular analysis 
#' on action properties for one interaction between 2 cell clusters.
#'
#' @rdname TgActionPairs-InterCell
#' @export
#'
setGeneric(name = "getTgActionPairs", def = function(object, ...) {
	standardGeneric("getTgActionPairs")
	}
)

#' @examples
#' \dontrun{
#'   getTgActionPairs(object)
#' }
#' @rdname TgActionPairs-InterCell
#'
setMethod(
	f = "getTgActionPairs",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@tg.action.pairs) || length(object@tg.action.pairs) == 0) {
			stop("TG action pairs not run [TODO] name alignment")
		}
		return(object@tg.action.pairs)
	}
)

#' Using TargetView VEinfo
#' 
#' The \code{setTgVEInfo} function is to \bold{set} the result of intercellular analysis
#' on detailed gene pairs (forming the interaction network) for one interaction between 2 cell clusters.
#'
#' @rdname TgVEInfo-InterCell
#' @export
#'
setGeneric(name = "setTgVEInfo", def = function(object, ...) {
	standardGeneric("setTgVEInfo")
	}
)

#' @param new.veinfo A new result of intercellular analysis on detailed gene pairs.
#'
#' @examples
#' \dontrun{
#'   setTgVEInfo(object, some.new.veinfo)
#' }
#'
#' @rdname TgVEInfo-InterCell
#'
setMethod(
	f = "setTgVEInfo",
	signature = "InterCell",
	definition = function(object, new.veinfo) {
		object@tg.veinfo <- new.veinfo
		return(object)
	}
)

#' Using TargetView VEinfo
#'
#' The \code{getTgVEInfo} function is to \bold{get} the result of intercellular analysis
#' on detailed gene pairs (forming the interaction network) for one interaction between 2 cell clusters.
#'
#' @rdname TgVEInfo-InterCell
#' @export
#'
setGeneric(name = "getTgVEInfo", def = function(object, ...) {
	standardGeneric("getTgVEInfo")
	}
)

#' @examples
#' \dontrun{
#'   getTgVEInfo(object)
#' }
#' @rdname TgVEInfo-InterCell
#'
setMethod(
	f = "getTgVEInfo",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@tg.veinfo) || length(object@tg.veinfo) == 0) {
			stop("TgVEinfo not run [TODO] name alignment")
		}
		return(object@tg.veinfo)
	}
)


#' Using TargetView Action Composition
#' 
#' The \code{setTgActionComp} function is to \bold{set} the result of intercellular analysis
#' on composition of action mode and action effect for one interaction between 2 cell clusters.
#'
#' @rdname TgActionComp-InterCell
#' @export
#'
setGeneric(name = "setTgActionComp", def = function(object, ...) {
	standardGeneric("setTgActionComp")
	}
)

#' @param new.action.comp A new result of intercellular analysis on composition of action mode and action effect.
#'
#' @examples
#' \dontrun{
#'   setTgActionComp(object, some.new.action.comp)
#' }
#'
#' @rdname TgActionComp-InterCell
#'
setMethod(
	f = "setTgActionComp",
	signature = "InterCell",
	definition = function(object, new.action.comp) {
		object@tg.action.comp <- new.action.comp
		return(object)
	}
)

#' Using TargetView Action Composition
#'
#' The \code{getTgActionComp} function is to \bold{get} the result of intercellular analysis
#' on composition of action mode and action effect for one interaction between 2 cell clusters.
#'
#' @rdname TgActionComp-InterCell
#' @export
#'
setGeneric(name = "getTgActionComp", def = function(object, ...) {
	standardGeneric("getTgActionComp")
	}
)

#' @examples
#' \dontrun{
#'   getTgActionComp(object)
#' }
#' @rdname TgActionComp-InterCell
#'
setMethod(
	f = "getTgActionComp",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@tg.action.comp) || length(object@tg.action.comp) == 0) {
			stop("TgActionComp not run [TODO] name alignment")
		}
		return(object@tg.action.comp)
	}
)


#' Using TargetView Pair Specificity
#' 
#' The \code{setTgSpGenes} function is to \bold{set} the result of intercellular analysis
#' on gene pair specificity for one interaction between 2 cell clusters.
#'
#' @rdname TgSpGenes-InterCell
#' @export
#'
setGeneric(name = "setTgSpGenes", def = function(object, ...) {
	standardGeneric("setTgSpGenes")
	}
)

#' @param new.spgenes A new result of intercellular analysis on gene pair specificity.
#'
#' @examples
#' \dontrun{
#'   setTgSpGenes(object, some.new.spgenes)
#' }
#'
#' @rdname TgSpGenes-InterCell
#'
setMethod(
	f = "setTgSpGenes",
	signature = "InterCell",
	definition = function(object, new.spgenes) {
		object@tg.spgenes <- new.spgenes
		return(object)
	}
)

#' Using TargetView Pair Specificity
#'
#' The \code{getTgSpGenes} function is to \bold{get} the result of intercellular analysis
#' on gene pair specificity for one interaction between 2 cell clusters.
#'
#' @rdname TgSpGenes-InterCell
#' @export
#'
setGeneric(name = "getTgSpGenes", def = function(object, ...) {
	standardGeneric("getTgSpGenes")
	}
)

#' @examples
#' \dontrun{
#'   getTgSpGenes(object)
#' }
#' @rdname TgSpGenes-InterCell
#'
setMethod(
	f = "getTgSpGenes",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@tg.spgenes) || length(object@tg.spgenes) == 0) {
			stop("TgSpGenes not run [TODO] name alignment")
		}
		return(object@tg.spgenes)
	}
)


#' Using Gene Pair Split
#' 
#' The \code{setGenePairSplit} function is to \bold{set} the gene pair split, 
#' which is either one charater or string, and used to split 2 gene partners in one gene pair. 
#'
#' @rdname GenePairSplit-InterCell
#' @export
#'
setGeneric(name = "setGenePairSplit", def = function(object, ...) {
		standardGeneric("setGenePairSplit")
	}
)

#' @param new.gene.pair.split A new gene pair split, which is either one character or string.
#'
#' @examples
#' \dontrun{
#'   setGenePairSplit(object, some.new.gene.pair.split)
#' }
#'
#' @rdname GenePairSplit-InterCell
#'
setMethod(
	f = "setGenePairSplit", 
	signature = "InterCell",
	definition = function(object, new.gene.pair.split) {
		if (length(new.gene.pair.split) != 1) {
			stop("Only one character or string is allowed, like '~' or '-~-'.")
		}
		# check if symbol are in fgenes
		use.genes <- unique(object@fgenes$gene)
		check.split.res <- strsplit(use.genes, split = new.gene.pair.split, fixed = TRUE)
		if (length(unlist(check.split.res)) > length(use.genes)) {
			stop("Used gene pair split symbol exists in gene names given in processing data (fgenes), and will be discarded.")
		}
		object@tool.vars$gene.pair.split <- new.gene.pair.split
		return(object)
	}
)


#' Using Gene Pair Split
#' 
#' The \code{getGenePairSplit} function is to \bold{get} the gene pair split, 
#' which is either one charater or string, and used to split 2 gene partners in one gene pair. 
#'
#' @rdname GenePairSplit-InterCell
#' @export
#'
setGeneric(name = "getGenePairSplit", def = function(object) {
		standardGeneric("getGenePairSplit")
	}
)

#' @examples
#' \dontrun{
#'   getGenePairSplit(object)
#' }
#' @rdname GenePairSplit-InterCell
#'
setMethod(
	f = "getGenePairSplit", 
	signature = "InterCell",
	definition = function(object) {
		return(object@tool.vars$gene.pair.split)
	}
)

#' Using Cluster Group Split
#' 
#' The \code{setClusterSplit} function is to \bold{set} the cluster group split, 
#' which is either one charater or string, and used to split 2 cell clusters in one interaction.
#'
#' @rdname ClusterSplit-InterCell
#' @export
#'
setGeneric(name = "setClusterSplit", def = function(object, ...) {
		standardGeneric("setClusterSplit")
	}
)

#' @param new.cluster.split A new cluster group split, which is either one character or string.
#'
#' @examples
#' \dontrun{
#'   setClusterSplit(object, some.new.cluster.split)
#' }
#'
#' @rdname ClusterSplit-InterCell
#'
setMethod(
	f = "setClusterSplit", 
	signature = "InterCell",
	definition = function(object, new.cluster.split) {
		if (length(new.cluster.split) != 1) {
			stop("Only one character or string is allowed, like '~' or '~~'.")
		}
		# check if symbol are in cluster names
		use.clusters <- unique(object@fgenes$cluster)
		check.split.res <- strsplit(use.clusters, split = new.cluster.split, fixed = TRUE)
		if (length(unlist(check.split.res)) > length(use.clusters)) {
			stop("Used cluster group split symbol exists in cluster name given in processing data (fgenes), and will be discarded.")
		}
		object@tool.vars$cluster.split <- new.cluster.split
		return(object)
	}
)

#' Using Cluster Group Split
#'
#' The \code{getClusterSplit} function is to \bold{set} the cluster group split, 
#' which is either one charater or string, and used to split 2 cell clusters in one interaction.
#'
#' @rdname ClusterSplit-InterCell
#' @export
#'
setGeneric(name = "getClusterSplit", def = function(object) {
		standardGeneric("getClusterSplit")
	}
)

#' @examples
#' \dontrun{
#'   getClusterSplit(object)
#' }
#' @rdname ClusterSplit-InterCell
#'
setMethod(
	f = "getClusterSplit", 
	signature = "InterCell",
	definition = function(object) {
		return(object@tool.vars$cluster.split)
	}
)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function that returns object
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# This is used to modify the InterCellDBPack object.
# like use high-confidence Experimentally validated database
#
# TO NOTE, exp > know > pred. Which means get <know> involved, <exp> score must be 0.
#
SelectDBSubset.default <- function(
	object,
	combined.score.range = c(1, 1000),
	use.exp = TRUE,
	exp.score.range = c(1, 1000),
	use.know = TRUE,
	know.score.range = c(1, 1000),
	use.pred = TRUE,
	pred.score.range = c(1, 1000),
	sel.physical = FALSE,
	sel.action.mode = "ALL",  # "ALL" consider "Other", so pairs not in actions.db will be reserved as well
	sel.action.effect = "ALL",
	sel.action.merge.option = "intersect",  # or "union"
	slim.along.with.pairs = TRUE  # genes.db will be exception
) {
	# check score range
	if (length(combined.score.range) != 2) {
		stop("Given `combined.score.range` not in right format, should be like c(1, 1000).")
	}
	if (length(exp.score.range) != 2) {
		stop("Given `exp.score.range` not in right format, should be like c(1, 1000).")
	}
	if (length(know.score.range) != 2) {
		stop("Given `know.score.range` not in right format, should be like c(1, 1000).")
	}
	if (length(pred.score.range) != 2) {
		stop("Given `pred.score.range` not in right format, should be like c(1, 1000).")
	}
	# check action 
	if (sel.action.mode[1] != "ALL") {
		not.valid.action.mode <- setdiff(sel.action.mode, kpred.mode)
		if (length(not.valid.action.mode) > 0) {
			warning("Given undefined action mode: ", paste0(not.valid.action.mode, collapse = ", ", ". "))
		}
		sel.action.mode <- intersect(sel.action.mode, kpred.mode)
		if (length(sel.action.mode) == 0) {
			stop("No valid action mode is selected!")
		}
	}
	if (sel.action.effect[1] != "ALL") {
		not.valid.action.effect <- setdiff(sel.action.effect, kpred.action.effect)
		if (length(not.valid.action.effect) > 0) {
			warning("Given undefined aciton effect: ", paste0(not.valid.action.effect, collapse = ", "), ". ")
		}
		sel.action.effect <- intersect(sel.action.effect, kpred.action.effect)
		if (length(sel.action.effect) == 0) {
			stop("No valid action effect is selected!")
		}
	}

	this.pairs.db <- object@pairs.db
	this.actions.db <- object@actions.db
	# top-range subset
	this.pairs.db <- this.pairs.db[intersect(which(this.pairs.db$inter.Combined.Score >= combined.score.range[1]),
		which(this.pairs.db$inter.Combined.Score <= combined.score.range[2])), ]
	if (nrow(this.pairs.db) == 0) {
		stop("Select too small range of score, and get no valid gene pairs!")
	}

	# select from scores
	retDB.list <- list()
	inds.no.exp.score <- which(this.pairs.db$inter.Experiments.Score == 0)
	inds.no.know.score <- which(this.pairs.db$inter.Database.Score == 0)
	# get exp part
	if (use.exp == TRUE) {
		retDB.list <- c(retDB.list, list(this.pairs.db[intersect(which(this.pairs.db$inter.Experiments.Score >= exp.score.range[1]),
			which(this.pairs.db$inter.Experiments.Score <= exp.score.range[2])), ]))
	}
	if (use.know == TRUE) {
		tmp.know.sub.db <- this.pairs.db[inds.no.exp.score, ]
		retDB.list <- c(retDB.list, list(tmp.know.sub.db[intersect(which(tmp.know.sub.db$inter.Database.Score >= know.score.range[1]),
			which(tmp.know.sub.db$inter.Database.Score <= know.score.range[2])), ]))
	}
	if (use.pred == TRUE) {
		tmp.pred.sub.db <- this.pairs.db[intersect(inds.no.exp.score, inds.no.know.score), ]
		retDB.list <- c(retDB.list, list(tmp.pred.sub.db[intersect(which(tmp.pred.sub.db$inter.Predicted.Score >= pred.score.range[1]),
			which(tmp.pred.sub.db$inter.Predicted.Score <= pred.score.range[2])), ]))
	}
	# collect result from score selection
	this.pairs.db <- Reduce(rbind, retDB.list)


	## select from actions
	use.action.pairs.list <- list()
	use.cut.symbol <- "->-"
	# in mode
	if (sel.physical == TRUE) {  # currently it's the same to set sel.action.mode => 'binding'
		if (sel.action.mode[1] == "ALL") {
			sel.action.mode <- "binding"
		} else {
			sel.action.mode <- unique(c(sel.action.mode, "binding"))
		}
	}
	if (sel.action.mode[1] != "ALL") {
		tmp.sel.mode.pairs.df <- FastAlignPairs(this.actions.db[which(this.actions.db$mode %in% sel.action.mode), 
			c("inter.GeneID.A", "inter.GeneID.B", "inter.GeneName.A", "inter.GeneName.B")], 4)
		tmp.sel.mode.pairs <- paste(tmp.sel.mode.pairs.df[, "inter.GeneID.A"], tmp.sel.mode.pairs.df[, "inter.GeneID.B"], sep = use.cut.symbol)
		use.action.pairs.list <- c(use.action.pairs.list, list(action.mode = tmp.sel.mode.pairs))
	}
	# in effect
	if (sel.action.effect[1] != "ALL") {
		# action
		tmp.col.action <- tapply(seq_len(nrow(this.actions.db)), this.actions.db$action, function(x) { x })
		# set list name right to be "non", "activation", "inhibition"
		ind.non.action <- setdiff(seq_along(tmp.col.action), match(c("activation", "inhibition"), names(tmp.col.action)))
		if (length(ind.non.action) != 1) {
			stop("Database Failed by incorporating undefined action mode.")
		}
		names(tmp.col.action)[ind.non.action] <- "non"
		# direction
		tmp.col.direction <- tapply(seq_len(nrow(this.actions.db)), this.actions.db$is_directional, function(x) { x })
		#tmp.col.a.act <- tapply(seq_len(nrow(this.actions.db)), this.actions.db$a_is_acting, function(x) { x })
		
		tmp.ret.inds <- integer()
		if ("positive" %in% sel.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["activation"]], tmp.col.direction[["t"]])))
		}
		if ("negative" %in% sel.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["inhibition"]], tmp.col.direction[["t"]])))	
		}
		if ("unspecified" %in% sel.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["non"]], tmp.col.direction[["t"]])))
		}
		if ("undirected" %in% sel.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["non"]], tmp.col.direction[["f"]])))	
		}
		tmp.sel.effect.pairs.df <- FastAlignPairs(this.actions.db[tmp.ret.inds, 
			c("inter.GeneID.A", "inter.GeneID.B", "inter.GeneName.A", "inter.GeneName.B")], 4)
		tmp.sel.effect.pairs <- paste(tmp.sel.effect.pairs.df[, "inter.GeneID.A"], tmp.sel.effect.pairs.df[, "inter.GeneID.B"], sep = use.cut.symbol)
		use.action.pairs.list <- c(use.action.pairs.list, list(action.effect = tmp.sel.effect.pairs))
	}
	# merge result
	use.action.pairs <- character()
	if (length(use.action.pairs.list) == 2 && all(c("action.mode", "action.effect") %in% names(use.action.pairs.list))) {
		use.action.pairs <- switch(sel.action.merge.option,
			"intersect" = {
				intersect(use.action.pairs.list[["action.mode"]], use.action.pairs.list[["action.effect"]])
			},
			"union" = {
				union(use.action.pairs.list[["action.mode"]], use.action.pairs.list[["action.effect"]])
			},
			stop("Undefined result merging option!")
		)
	} else {
		use.action.pairs <- as.character(unlist(use.action.pairs.list))
	}
	use.action.pairs <- unique(use.action.pairs)


	# collect result from actions selection
	if (length(use.action.pairs) > 0) {
		this.pairs.db <- FastAlignPairs(this.pairs.db, 4)
		tmp.cur.pairs <- paste(this.pairs.db[, "inter.GeneID.A"], this.pairs.db[, "inter.GeneID.B"], sep = use.cut.symbol)
		this.pairs.db <- this.pairs.db[which(tmp.cur.pairs %in% use.action.pairs), ]
	}

	## get pairs settled, and slim other accessory database
	object@pairs.db <- this.pairs.db
	if (slim.along.with.pairs == TRUE) {
		use.cut.symbol.slim <- "->-"
		#involved.genes <- unique(c(object@pairs.db$inter.GeneName.A, object@pairs.db$inter.GeneName.B))
		tmp.sel.pairs.conv <- paste(object@pairs.db[, "inter.GeneID.A"], object@pairs.db[, "inter.GeneID.B"], sep = use.cut.symbol.slim)
		tmp.sel.pairs.rev <- paste(object@pairs.db[, "inter.GeneID.B"], object@pairs.db[, "inter.GeneID.A"], sep = use.cut.symbol.slim)
		# slim actions
		tmp.sel.actions.p <- paste(object@actions.db[, "inter.GeneID.A"], object@actions.db[, "inter.GeneID.B"], sep = use.cut.symbol.slim)
		object@actions.db <- object@actions.db[union(which(tmp.sel.actions.p %in% tmp.sel.pairs.conv), which(tmp.sel.actions.p %in% tmp.sel.pairs.rev)), ]
	}

	return(object)
}


setGeneric(name = "SelectDBSubset", def = function(object, ...) {
		standardGeneric("SelectDBSubset")
	}
)
setMethod(
	f = "SelectDBSubset", 
	signature = "InterCellDBPack",
	definition = SelectDBSubset.default
)
setMethod(
	f = "SelectDBSubset", 
	signature = "InterCell",
	definition = function(object, ...) {
		object@database <- SelectDBSubset.default(object@database, ...)
		return(object)
	}
)





# %%%%%%%%%%%%%%%%%%
# Other functions
# %%%%%%%%%%%%%%%%%%

#' Remap gene symbols to their authorized genenames
#'
#' @description
#' Use species-specific reference database to remap genes to their authorized genenames
#'
#' @param markers.all Data.frame. Feature genes that are generated from \code{Seurat::FindAllMarkers()} or 
#' similar functions in other packages.
#' @param species [TODO]
#' @param if.used.inside Logic. If used inside, some process will not run.
#'
#'
#'
#' @export
#'
DataPrep.RemapClustersMarkers  <- function(
	markers.all,
	species,
	if.used.inside = FALSE
) {
	# force character
	markers.all$gene <- as.character(markers.all$gene)
	markers.all$cluster <- as.character(markers.all$cluster)

	## set each database
	# Check species, only support human and mouse yet
	genes.ref.db <- NULL
	species <- CheckSpeciesValidity(species)
	if (species == "human") {
		genes.ref.db <- genes.human.ref.db
	} else {
		genes.ref.db <- genes.mouse.ref.db
	}
	
	entrez.db <- genes.ref.db$gene.ncbi.db
	map.synonyms.db <- genes.ref.db$gene.synonyms.db

	# split genes to already authorized ones and un-authorized ones
	inds.raw.match <- which(markers.all$gene %in% entrez.db$Symbol_from_nomenclature_authority)
	markers.raw.match   <- markers.all[inds.raw.match, ]
	markers.raw.unmatch <- markers.all[setdiff(seq_len(nrow(markers.all)), inds.raw.match), ]
	proc.unmatch.genes <- unique(markers.raw.unmatch$gene)
	
	## un-authorized genes have 2 ways to go:
	# 1. cannot mapping from synonyms, keep still
	# 2. mapping from synonyms
	#
	ret.d0.cannot.match <- ret.d1.map.to.diff <- ret.d2.diff.map.one <- ret.d3.map.to.exist <- character()
	unmatch.map.results <- data.frame(unmatched = character(), match.res = character(), stringsAsFactors = FALSE)
	if (nrow(markers.raw.unmatch) > 0) {
		# 1> extract unmatch-able genes
		ret.d0.cannot.match <- setdiff(proc.unmatch.genes, map.synonyms.db$Synonym.each)
		
		# short warning
		if (length(ret.d0.cannot.match) > 0 && if.used.inside == FALSE) {
			print(paste0("Get ", length(ret.d0.cannot.match), 
				" genes cannot be mapped from synonyms.",
				" See return value `$unmapping.genes`."))
		}

		# 2> mapping the rest
		genes.poss.map <- setdiff(proc.unmatch.genes, ret.d0.cannot.match)
		inds.map.possible <- which(map.synonyms.db$Synonym.each %in% genes.poss.map)
		
		# check un-identical synonyms. One synonyms could map to different authorized genes
		check.non.identical <- tapply(map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.possible], map.synonyms.db$Synonym.each[inds.map.possible], function(x) {x}, simplify = FALSE)
		inds.dups <- sapply(check.non.identical, function(x) { ifelse(length(x) > 1, TRUE, FALSE) })
		# short warning
		ret.d1.warning <- names(check.non.identical)[inds.dups]
		if (length(ret.d1.warning) > 0 && if.used.inside == FALSE) {
			print(paste0("Get ", length(ret.d1.warning), 
				" synonyms that cannot be mapped to unified authorized genes.",
				" See return value `$one.synonym.map.to.some.genes`."))  
		}
		# detailed hint
		ret.d1.map.to.diff <- vapply(which(inds.dups == TRUE), all.non.identical = check.non.identical, 
			function(x, all.non.identical) {
				this.gene <- names(all.non.identical)[x]
				against.authorized.genes <- paste0("(", paste0(all.non.identical[[x]], collapse = ", "), ")")
				paste(this.gene, "-->", against.authorized.genes)
			},
			FUN.VALUE = character(1), USE.NAMES = FALSE
		)

		# check if mapping result the same as other mapped genes
		check.map.dups <- tapply(map.synonyms.db$Synonym.each[inds.map.possible], 
			map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.possible], 
			function(x) {x}, simplify = FALSE)
		inds.map.dups <- sapply(check.map.dups, function(x) { ifelse(length(x) > 1, TRUE, FALSE) })
		# short warning
		ret.d2.warning <- unique(as.character(unlist(lapply(which(inds.map.dups == TRUE), all.map.dups = check.map.dups, 
			function(x, all.map.dups) {
				all.map.dups[[x]]
			}
		))))
		if (length(ret.d2.warning) > 0 && if.used.inside == FALSE) {
			print(paste0("Get ", length(ret.d2.warning), 
				" synonyms that get overlap mapping genes with at least one other.",
				" See return value `$some.synonyms.map.to.one.gene`."))  
		}
		# detailed hint
		ret.d2.diff.map.one <- vapply(which(inds.map.dups == TRUE), all.map.dups = check.map.dups,
			function(x, all.map.dups) {
				this.tg <- names(all.map.dups)[x]
				from.synonyms <- paste0("(", paste0(all.map.dups[[x]], collapse = ", "), ")")
				paste(from.synonyms, "-->", this.tg)
			},
			FUN.VALUE = character(1), USE.NAMES = FALSE
		)

		# check if mapping to exist genes
		check.map.to.exist <- intersect(markers.raw.match$gene, names(check.map.dups))
		# short warning
		ret.d3.warning <- as.character(unlist(lapply(check.map.to.exist, all.map.dups = check.map.dups, 
			function(x, all.map.dups) {
				all.map.dups[[which(names(all.map.dups) == x)]]
			}
		)))
		if (length(ret.d3.warning) > 0 && if.used.inside == FALSE) {
			print(paste0("Get ", length(ret.d3.warning), 
				" synonyms that get mapping to existing authorized genes.",
				" See return value `$synonyms.map.to.exist.gene`."))
		}
		# detailed hint
		ret.d3.map.to.exist <- vapply(check.map.to.exist, all.map.dups = check.map.dups, 
			function(x, all.map.dups) {
				this.to.map.genes <- all.map.dups[[which(names(all.map.dups) == x)]]
				orig.to.map <- paste0("(", paste0(this.to.map.genes, collapse = ", "), ")")
				paste(orig.to.map, "-->", x)
			},
			FUN.VALUE = character(1), USE.NAMES = FALSE
		)
		
		## get mapping result
		# unmatch-able ones
		unmatch.map.result.0 <- data.frame(unmatched = ret.d0.cannot.match, match.res = ret.d0.cannot.match, stringsAsFactors = FALSE)
		# match-able ones, in default: the first matched gene name will be used 
		inds.map.matches <- match(genes.poss.map, map.synonyms.db$Synonym.each)
		unmatch.map.result.1 <- data.frame(unmatched = map.synonyms.db$Synonym.each[inds.map.matches], match.res = map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.matches], stringsAsFactors = FALSE)
		# collect all
		unmatch.map.results <- rbind(unmatch.map.result.0, unmatch.map.result.1)
	}
	
	# collect markers result after mapping
	markers.raw.unmatch.dummy <- left_join(markers.raw.unmatch[, "gene", drop = FALSE],
		unmatch.map.results, by = c("gene" = "unmatched"))
	markers.raw.unmatch$gene <- markers.raw.unmatch.dummy$match.res
	markers.all <- rbind(markers.raw.unmatch, markers.raw.match)

	# after remapping, genes get to be duplicate with existing ones, re-check if mapping result has duplicate genes
	fcheck.result <- lapply(unique(markers.all$cluster), ref.markers = markers.all,
		function(x, ref.markers) {
			this.c.markers <- ref.markers[which(ref.markers$cluster == x), ]
			this.c.len <- tapply(seq_along(this.c.markers$gene), this.c.markers$gene, length)
			this.f.dup.genes <- names(this.c.len)[which(this.c.len > 1)]
			# in default: remove the latter one in given data
			this.c.markers <- DoPartUnique(this.c.markers, match(c("gene", "cluster"), colnames(this.c.markers)))
			list(dup.genes = this.f.dup.genes, markers = this.c.markers)
		})
	names(fcheck.result) <- unique(markers.all$cluster)

	# ret 
	ret.markers.all <- bind_rows(lapply(fcheck.result, function(x) { x$markers }))
	# detailed hint
	ret.dx.fcheck.dup <- lapply(seq_along(fcheck.result), all.fcheck.res = fcheck.result,
		function(x, all.fcheck.res) {
			this.cluster <- names(all.fcheck.res)[x]
			this.dup.genes <- all.fcheck.res[[x]][["dup.genes"]]
			ret.detailed <- NA
			if (length(this.dup.genes) > 0) {
				ret.detailed <- paste(this.cluster, "~", paste0(this.dup.genes, collapse = ", "))
			}
			list(detailed = ret.detailed, raw = unique(this.dup.genes))
		}
	)
	un.fin.check <- as.character(sapply(ret.dx.fcheck.dup, function(x) { x$detailed }))
	un.fin.check <- un.fin.check[which(!is.na(un.fin.check))]
	if (length(un.fin.check) > 0) {
		warning("There remain genes group by clusters to be checked manually. ",
			paste0(un.fin.check, collapse = "; "),
			". The program automatically remove duplicate ones in cluster scale.")
	}
	ret.d4.final.dup.genes <- lapply(ret.dx.fcheck.dup, function(x) { x$raw })
	names(ret.d4.final.dup.genes) <- names(fcheck.result)
	ret.d4.len <- vapply(ret.d4.final.dup.genes, FUN = length, FUN.VALUE = integer(1))
	ret.d4.final.dup.genes <- ret.d4.final.dup.genes[which(ret.d4.len > 0)]

	return(list(result = ret.markers.all, 
		unmapping.genes = ret.d0.cannot.match, 
		one.synonym.map.to.some.genes = ret.d1.map.to.diff, 
		some.synonyms.map.to.one.gene = ret.d2.diff.map.one, 
		synonyms.map.to.exist.gene = ret.d3.map.to.exist,
		after.map.dup.genes = ret.d4.final.dup.genes))
}



ListAllClusters <- function(
	object
) {
	return(unique(object@fgenes$cluster))
}



ReplaceClusterName <- function(
	object,
	cluster.names.current,
	cluster.names.replace
) {
	markers.all <- object@fgenes
	colnames.cluster <- "cluster"
	# process
	if ((colnames.cluster %in% colnames(markers.all)) == FALSE) {
		stop("Selected column name defining clusters is not in the given data!")
	}
	if (length(cluster.names.current) != length(cluster.names.replace)) {
		stop("The replaced names are of different length of the current used ones.")
	}
	reserve.oldnames.col <- reserve.oldnames.col.proto <- paste(colnames.cluster, "oldv", sep = ".")
	for (try.i in 1:100) {
		reserve.oldnames.col <- paste(reserve.oldnames.col.proto, as.character(try.i), sep = ".")
		if ((reserve.oldnames.col %in% colnames(markers.all)) == FALSE) {
			break
		}
		if (try.i == 100) {
			stop("Cannot allocate proper colnames for old cluster names, program failed! Please check given data!")
		}
	}
	markers.all[, reserve.oldnames.col] <- markers.all[, colnames.cluster]
	tmp.fac <- factor(markers.all[, colnames.cluster])
	lvl.tmp.fac <- levels(tmp.fac)
	inds.match <- match(cluster.names.current, lvl.tmp.fac)
	if (length(which(is.na(inds.match))) != 0) {
		stop(paste0("Please give right current used cluster names! ",
			"Wrong given ones are: ", paste0(cluster.names.current[which(is.na(inds.match))], collapse = ", "),
			"."))    
	}
	lvl.tmp.fac[inds.match] <- cluster.names.replace
	levels(tmp.fac) <- lvl.tmp.fac
	markers.all[, colnames.cluster] <- as.character(tmp.fac)
	# return
	object@fgenes <- markers.all
	return(object)
}


# fetch genes of interest
FetchGeneOI.default <- function(
	object,  # InterCellDBPack
	select.location = NULL,
	select.location.score = c(1:5),
	select.type = NULL,
	select.merge.type = NULL,  # merged type only have 16 options. Prioritize upon select.type
	select.go.terms = NULL,  # ID and Term are supported
	go.use.relative = TRUE, 
	go.relative.option = "offspring"
) {
	#this.genes.db <- object@genes.db
	#this.pairs.db <- object@pairs.db
	#this.locs.db <- object@anno.location.db
	#this.type.db <- object@anno.type.db
	#this.go.db <- object@go.ref.db

	## input parameter process
	# check location
	if (!is.null(select.location)) {
		avb.opt.location <- unique(object@anno.location.db$GO.Term.target)
		not.valid.location <- setdiff(select.location, avb.opt.location)
		if (length(not.valid.location) > 0) {
			warning("Given undefined location: ", paste0(not.valid.location, collapse = ", ", ". "))
		}
		select.location <- intersect(select.location, avb.opt.location)
		if (length(select.location) == 0) {
			select.location <- NULL
		}
	}
	if ((length(select.location.score) > 1 && !is.integer(select.location.score)) || 
		(length(select.location.score) == 1 && !is.numeric(select.location.score))) {
		stop("Location score ranges from 1 to 5, and only those 5 integers are supported!")
	}
	if (!is.null(select.location.score)) {
		avb.location.score <- c(1:5)
		not.valid.location.score <- setdiff(select.location.score, avb.location.score)
		if (length(not.valid.location.score) > 0) {
			warning("Given invalid location score: ", paste0(not.valid.location.score, collapse = ", ", ". "))
		}
		select.location.score <- intersect(select.location.score, avb.location.score)
	}

	# check type
	if (!is.null(select.type)) {
		avb.opt.type <- unique(object@anno.type.db$Keyword.Name)
		not.valid.type <- setdiff(select.type, avb.opt.type)
		if (length(not.valid.type) > 0) {
			warning("Given undefined type: ", paste0(not.valid.type, collapse = ", ", ". "))
		}
		select.type <- intersect(select.type, avb.opt.type)
		if (length(select.type) == 0) {
			select.type <- NULL
		}
	}

	# check merged type (if given, then override type)
	if (!is.null(select.merge.type)) {
		if (!is.null(select.type)) {
			warning("Select genes using both parameter 'select.type' & 'select.merge.type'. Only the options in 'select.merge.type' are used!")
		}
		avb.opt.mg.type <- unique()
	}

	# check GO terms
	# - put in the tool function


	## fetch genes
	ret.gene.oi <- character()
	inside.set.gene.oi <- function(gene.oi.res, new.get.genes) {
		if (length(gene.oi.res) == 0) {
			gene.oi.res <- new.get.genes
		} else {
			gene.oi.res <- intersect(gene.oi.res, new.get.genes)
		}
		return(gene.oi.res)
	}
	# use location
	if (!is.null(select.location)) {
		gene.oi.from.locs <- object@anno.location.db[intersect(which(object@anno.location.db$GO.Term.target %in% select.location), 
			which(object@anno.location.db$score %in% select.location.score)), "Gene.name"]
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.locs)
	}
	if (!is.null(select.type)) {
		gene.oi.from.type <- object@anno.type.db[which(object@anno.type.db$Keyword.Name %in% select.type), "Gene.name"]
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.type)
	}
	if (!is.null(select.go.terms)) {
		gene.oi.from.go <- Tool.FindGenesFromGO(select.go.terms, object@genes.db, object@go.ref.db, 
			go.use.relative = go.use.relative, go.relative.option = go.relative.option)
		gene.oi.from.go <- as.character(unlist(gene.oi.from.go))
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.go)
	}
	# get unique result
	ret.gene.oi <- unique(ret.gene.oi)

	print(paste("Fetch", length(ret.gene.oi), "genes of interest."))
	return(ret.gene.oi)
}

setGeneric(name = "FetchGeneOI", def = function(object, ...) {
		standardGeneric("FetchGeneOI")
	}
)
setMethod(
	f = "FetchGeneOI", 
	signature = "InterCellDBPack",
	definition = FetchGeneOI.default
)
setMethod(
	f = "FetchGeneOI", 
	signature = "InterCell",
	definition = function(object, ...) {
		return(FetchGeneOI.default(object@database, ...))
	}
)



#' Format User Given Gene Pairs
#' 
#' This function makes 2-column table (each row represents one gene pair) to be 
#' formatted to be used in this package.
#' The gene names will be remapped to the authorized gene names recorded in gene 
#' reference database embedded in this package.
#'
#' @param gene.pairs.table [TODO] in data.frame 2 columns
#' @param species  [TODO] which species to be add database
#'
#' @return A list with \code{$result} storing the formatted gene pairs.
#'
#' @export
#'
FormatCustomGenePairs <- function(
	gene.pairs.table,
	species
) {
	# check input
	if (class(gene.pairs.table) != "data.frame") {
		stop("Given gene pairs must be stored in table (R `data.frame` structure).")
	}
	if (ncol(gene.pairs.table) != 2) {
		warning("Given table has more than 2 columns, only the first 2 columns will be used.")
	}

	# check species
	species <- CheckSpeciesValidity(species)

	# process
	result.list <- list()
	apx.list <- list()
	for (i in 1:2) {  # only use the first 2 columns
		tmp.genes <- gene.pairs.table[, i]
		# create dummy df to meet the requirement of musthave columns
		dummy.fgenes <- data.frame(gene = tmp.genes, 
			cluster = seq_along(tmp.genes),  # give one gene one cluster to make duplicate ones preserved, but program may get quite slow
			LogFC = 1, PVal = 0, 
			num.id = seq_along(tmp.genes), stringsAsFactors = FALSE)
		dummy.remap.res <- suppressWarnings(DataPrep.RemapClustersMarkers(dummy.fgenes, species, if.used.inside = TRUE))
		dummy.fgenes <- dummy.remap.res$result
		dummy.fgenes <- dummy.fgenes[order(dummy.fgenes$num.id, decreasing = FALSE), ]
		result.list <- c(result.list, list(dummy.fgenes$gene))
		apx.list <- c(apx.list, dummy.remap.res[setdiff(names(dummy.remap.res), "result")])
	}

	result.IT <- data.frame(gene.A = result.list[[1]],
		gene.B = result.list[[2]],
		stringsAsFactors = FALSE)
	colnames(result.IT) <- colnames(gene.pairs.table)[1:2]

	# return
	list(result = result.IT, match.status = apx.list)
}



