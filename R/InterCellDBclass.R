
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

#' Predefined Action Mode & Color Usage for It
#'
#' \code{kpred.action.mode} is extracted from databases depicting mode of actions for functional usages.
#' expression is transcriptional regulation.
#'
#' @rdname pred-action-mode
#' @export
#'
kpred.action.mode <- c("activation", "inhibition", "binding", "catalysis", "reaction", "expression", "ptmod", "other")

#' Predefined Action Mode & Color Usage for It 
#'
#' \code{kpred.color.mode} is aligned with \code{kpred.action.mode} one-to-one. In default setting,
#' the color used for \bold{action mode} in plotting will use colors listed in this variable.
#'
#' @rdname pred-action-mode
#' @export
#'
kpred.color.mode <- c("#D70051", "#00913A", "#1296D4", "#956134", "#F46D42", "#0A0AFF", "#762A83", "#B5B5B6")
#c("#FB8072", "#B3DE69", "#80B1D3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDB462", "#FCCDE5")

#' Predefined Action Effect & Color Usage for It
#'
#' \code{kpred.action.effect} is extracted from databases depicting action effect of actions for functional usages.
#'
#' @rdname pred-action-effect
#' @export
#'
kpred.action.effect <- c("positive", "negative", "unspecified", "undirected")

#' Predefined Action Effect & Color Usage for It
#'
#' \code{kpred.color.effect} is aligned with \code{kpred.action.effect} one-to-one. In default setting,
#' the color used for \bold{action effect} in plotting will use colors listed in this variable.
#'
#' @rdname pred-action-effect
#' @export
#'
kpred.color.effect <- c("#FB8072", "#B3DE69", "#80B1D3", "#8DD3C7")

#' Predefined Action Effect & Color Usage for It
#'
#' \code{kpred.ext.action.effect} is the \bold{extended} format, which depicts action effect and direction of action together, 
#' and extends action effect to 7 different types.
#'
#' @rdname pred-action-effect
#' @export
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

#' Predefined Action Effect & Color Usage for It
#'
#' \code{kpred.color.ext.effect} is aligned with \code{kpred.ext.action.effect} one-to-one. In default setting,
#' the color used for \bold{extened action effect} in plotting will use colors listed in this variable.
#'
#' @rdname pred-action-effect
#' @export
#'
kpred.color.ext.effect <- c("#8DD3C7", "#FB8072", "#FF5740", "#B3DE69", "#81EF48", "#80B1D3", "#6B6AEA")



#
TgView.formula.onExprs.default <- function(
	data.f,
	data.b
) {
	if (length(data.f) != length(data.b)) {
		stop("Unexpected non-identical length data input!")
	}
	return(data.f * data.b)
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
#' @param pval.log.max The number for replacing the infinite number when \code{log(PVal = 0)}.
#'
TgView.formula.onPVal.default <- function(
	data.f, 
	data.b,
	pval.log.max = 300
) {
	# 'pval.log.max' is 300 in default setting.  
	# The use of 300 as default maximum, because e-303 are usual lowest limit when all log(values) are infinite.
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
#'  The given data need to have columns named 'gene', 'cluster', 'LogFC', 'PVal'. Those each means, 'gene': the 
#'  authorized gene names, 'cluster': cell cluster, 'LogFC': the log transformed fold change of genes, 'PVal': the 
#'  P val for regarding one gene as differentially expressed one. 
#'  The expression level of each gene is also needed when statistical test is required and it should be stored in column named 'Exprs'.
#' @slot database Database stored in \code{InterCellDBPack} object. 
#' @slot formulae The default formulae to be used in analysis. 
#' @slot inter.fullview This is used to store the result of network analysis. 
#' @slot tg.itinfo This stores information about one intercellular communication.
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
		tg.itinfo = "list",
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
	# if (length(object@pred.action$action.mode) == 0 || !all(object@pred.action$action.mode %in% kpred.action.mode)) {
	# 	return(paste0("Using undefined action mode: ", paste0(setdiff(object@pred.action$action.mode, kpred.action.mode), collapse = ", ")))
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
	if (!all(kmusthave.colnames %in% object@misc$musthave.colnames)) {
		return(paste0("Detect manual modification on object internal variables. ", 
			"Program will be failed in unexpected situations. "))
	}

	validInterCellDBPackObject(object@database)
	TRUE
}
setValidity("InterCell", validInterCellObject)



setMethod(
	f = "initialize",
	signature = c("InterCellDBPack"),
	definition = function(.Object, ...) {
		.Object <- methods::callNextMethod(.Object, ...)
		.Object@accessory.db <- list(merge.type.list = Uniprot.key.map.list)
		.Object@misc <- list(TAKEN = "nothing yet")
		methods::validObject(.Object)
		return(.Object)
	}
)

# initialize \code{InterCell-class}
setMethod(
	f = "initialize",
	signature = c("InterCell"),
	definition = function(.Object, ...) {
		.Object <- methods::callNextMethod(.Object, ...)
		.Object@formulae <- list(
			TG.EXPRS = TgView.formula.onExprs.default,
			TG.LOGFC = TgView.formula.onLogFC.default,
			TG.PVAL = TgView.formula.onPVal.default)
		#.Object@pred.action <- list(action.mode = kpred.action.mode, action.effect = kpred.action.effect)
		#.Object@tool.vars <- list(gene.pair.split = "-~-", cluster.split = "~")
		methods::validObject(.Object)
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

		# show the used database
		show(object@database)

		# show if intercellular analysis is processed
		if (!is.null(object@tg.veinfo) && length(object@tg.veinfo) > 0) {
			this.involved.clusters <- getOrigClusterNameTgVEInfo(object)
			cat("Intercellular analysis performed between ", this.involved.clusters$cluster.name.A, " and ", this.involved.clusters$cluster.name.B)
		} else {
			cat("Intercellular analysis not processed yet.")
		}
		cat("\n")
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
#' @return A \code{InterCellDBPack} object.
#'
#' @import methods
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
#' @param add.exprs It decides whether to add gene expression information.
#' @param exprs.data It gives the normalized count matrix as expression information.
#' @param force.write.exprs If there is column 'Exprs' in parameter \code{DEG.table}, users should 
#'  decide whether to overwrite the existed one (set TRUE) or not (set FALSE).
#' @param remap.genes This decides whether to use \pkg{InterCellDB} integrated gene database to
#' standardize all genes to get more perfectly matched with protein/gene pairs database. Default is \code{FALSE}.
#' @param cluster.split The letters used to split 2 cell clusters in one interaction. It can also be modified later by 
#' using \code{setClusterSplit} if not decided yet.
#' @param gene.pair.split The letters used to split 2 gene partners in one gene pair. It can also be modified later by 
#' using \code{setGenePairSplit} if not decided yet.
#'
#'
#' @details
#' The parameter \code{DEG.table} is recommended to be generated by \pkg{Seurat}. Other packages are also applicable, if
#' they could handle scRNA-seq data, do cell clustering and do calculation on cluster-specific differentially expressed genes. The input 
#' format of \code{DEG.table} should be one data.frame with 4 required columns that are named 'cluster', 'gene', 'LogFC', 'PVal'.
#' \itemize{
#'   \item{cluster}: the cell cluster.
#'   \item{gene}: differentially expressed genes, which are grouped by their belonging clusters.
#'   \item{LogFC}: the fold change of genes.
#'   \item{PVal}: the P value for gene being calculated as differentially expressed gene.
#' }
#' 
#' Gene expression data can be added when parameter \code{add.exprs} is set TRUE. It will be stored in column 'Exprs' and 
#' only those genes given in parameter \code{DEG.table[, "gene"]} will be perserved.
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
#' @return A \code{InterCell} object.
#'
#' @importFrom methods new show callNextMethod validObject
#'
#' @export
#'
CreateInterCellObject <- function(
	DEG.table,
	species,
	add.exprs = FALSE,
	exprs.data = NULL,  # could be expression matrix or seurat[["RNA"]]@data, etc
	force.write.exprs = FALSE,
	remap.genes = FALSE,
	cluster.split = "~",
	gene.pair.split = "#~#"
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

	# set database object
	DBPack.obj <- CreateDBPackObject(species)

	# add expression data when needed
	tmp.musthave.new <- kmusthave.colnames
	DEG.align.res <- NA
	if (add.exprs == TRUE) {
		DEG.table <- DataPrep.AddExprs(DEG.table, exprs.data, force.write.exprs)
		tmp.musthave.new <- c(tmp.musthave.new, "Exprs")
		if (remap.genes == TRUE) {
			dummy.gdf <- data.frame(gene = rownames(exprs.data), ID = seq_len(nrow(exprs.data)), cluster = "IT", stringsAsFactors = FALSE)
			DEG.align.res <- DataPrep.RemapClustersMarkers(dummy.gdf, species, final.dup.rm = FALSE)
			colnames(dummy.gdf)[1] <- c("orig.gene")
			dummy.gdf[, "mapped.gene"] <- DEG.align.res$result[order(DEG.align.res$result$ID), "gene"]
			# use this to remap genes in every cluster
			DEG.table$gene <- dummy.gdf[match(DEG.table$gene, dummy.gdf[, "orig.gene"]), "mapped.gene"]
		}
	} else {
		# re-align input table with gene reference database in InterCellDB
		if (remap.genes == TRUE) {
			DEG.align.res <- DataPrep.RemapClustersMarkers(DEG.table, species)
			DEG.table <- DEG.align.res$result
		}
	}

	
	IT.InterCell.Obj <- new(
		Class = "InterCell",
		fgenes = DEG.table,
		database = DBPack.obj,
		#pred.action = list(action.mode = kpred.action.mode, action.effect = kpred.action.effect,
		#	color.mode = kpred.color.mode, color.effect = kpred.color.effect),
		tool.vars = list(gene.pair.split = gene.pair.split, cluster.split = cluster.split), 
		misc = list(musthave.colnames = tmp.musthave.new, if.remap.genes = remap.genes)
	)
	# give the result of alignment in the misc
	IT.InterCell.Obj@misc$input.align.result <- DEG.align.res
	# set object
	IT.InterCell.Obj <- setGenePairSplit(IT.InterCell.Obj, gene.pair.split)
	IT.InterCell.Obj <- setClusterSplit(IT.InterCell.Obj, cluster.split)

	return(IT.InterCell.Obj)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Accessory Function
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Using Reference Database
#'
#' The \code{setRefDatabase} function is to \bold{set} the reference database for \code{InterCell}.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
#'
#' @rdname RefDatabase-InterCell
#' @export
#'
setGeneric(name = "setRefDatabase", def = function(object, ...) {
	standardGeneric("setRefDatabase")
	}
)

#' @param new.ref.database A new reference database.
#'
#' @examples
#' \dontrun{
#'   setRefDatabase(object, some.new.ref.database)
#' }
#'
#' @rdname RefDatabase-InterCell
#'
setMethod(
	f = "setRefDatabase",
	signature = "InterCell",
	definition = function(object, new.ref.database) {
		if (class(new.ref.database) != "InterCellDBPack") {
			stop("Given new reference database is not in right format.")
		}
		print("Change reference database.")
		object@database <- new.ref.database
		return(object)
	}
)

#' Using Reference Database
#'
#' The \code{getRefDatabase} function is to \bold{get} the reference database for \code{InterCell}.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
#'
#' @rdname RefDatabase-InterCell
#' @export
#'
setGeneric(name = "getRefDatabase", def = function(object, ...) {
	standardGeneric("getRefDatabase")
	}
)

#'
#' @examples
#' \dontrun{
#'   getRefDatabase(object)
#' }
#'
#' @rdname RefDatabase-InterCell
#'
setMethod(
	f = "getRefDatabase",
	signature = "InterCell",
	definition = function(object) {
		if (class(object@database) != "InterCellDBPack") {
			stop("Unexpected error: Reference database is invalid!")
		}
		return(object@database)
	}
)


#' Using FullView Result
#' 
#' The \code{setFullViewResult} function is to \bold{set} the result of network analysis.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			warning("Overwrite existed result on network analysis. Former downstream result is cleaned.")
		}
		object@inter.fullview <- new.inter.fullview
		object@tg.action.pairs <- list()
		object@tg.action.comp <- list()
		object@tg.itinfo <- list()
		object@tg.veinfo <- list()
		object@tg.spgenes <- list()
		return(object)
	}
)

#' Using FullView Result
#'
#' The \code{getFullViewResult} function is to \bold{get} the result of network analysis.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			stop("Network analysis is not performed yet. Please use function `AnalyzeInterInFullView` to generate it.")
		}
		return(object@inter.fullview)
	}
)

#' Using Target Plain Information
#' 
#' The \code{setTgItInfo} function is to \bold{set} plain information about gene pairs,
#' which are directly fetched from full view analysis.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
#'
#' @rdname TgItInfo-InterCell
#' @export
#'
setGeneric(name = "setTgItInfo", def = function(object, ...) {
	standardGeneric("setTgItInfo")
	}
)

#' @param new.itinfo A new plain information on gene pairs, which are directly fetched from full view analysis.
#'
#' @examples
#' \dontrun{
#'   setTgItInfo(object, some.new.itinfo)
#' }
#'
#' @rdname TgItInfo-InterCell
#'
setMethod(
	f = "setTgItInfo",
	signature = "InterCell",
	definition = function(object, new.itinfo) {
		object@tg.itinfo <- new.itinfo
		return(object)
	}
)

#' Using Target Plain Information
#' 
#' The \code{getTgItInfo} function is to \bold{get} plain information about gene pairs,
#' which are directly fetched from full view analysis.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
#'
#' @rdname TgItInfo-InterCell
#' @export
#'
setGeneric(name = "getTgItInfo", def = function(object, ...) {
	standardGeneric("getTgItInfo")
	}
)

#' @examples
#' \dontrun{
#'   getTgItInfo(object)
#' }
#' @rdname TgItInfo-InterCell
#'
setMethod(
	f = "getTgItInfo",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@tg.itinfo) || length(object@tg.itinfo) == 0) {
			stop("No plain information between specific 2 cells is fetched. Please use `FetchInterOI` to generate that. ")
		}
		return(object@tg.itinfo)
	}
)

#' Using TargetView Action Pairs
#' 
#' The \code{setTgActionPairs} function is to \bold{set} the result of intercellular analysis 
#' on action properties for one interaction between 2 cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			stop("The action pairs for selected interaction is not generated. Please use `FetchInterOI` to generate that. ")
		}
		return(object@tg.action.pairs)
	}
)

#' Using TargetView VEinfo
#' 
#' The \code{setTgVEInfo} function is to \bold{set} the result of intercellular analysis
#' on detailed gene pairs (forming the interaction network) for one interaction between 2 cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			stop("No interaction between specific 2 cells is fetched. Please use `FetchInterOI` to generate that. ")
		}
		return(object@tg.veinfo)
	}
)


#' Using TargetView Action Composition
#' 
#' The \code{setTgActionComp} function is to \bold{set} the result of intercellular analysis
#' on composition of action mode and action effect for one interaction between 2 cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			stop("The analysis of composition of action properties is not performed. Please use `AnalyzeInterInAction` to generate that. ")
		}
		return(object@tg.action.comp)
	}
)


#' Using TargetView Pair Specificity
#' 
#' The \code{setTgSpGenes} function is to \bold{set} the result of intercellular analysis
#' on gene pair specificity for one interaction between 2 cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
			stop("Analysis of gene pair specificity is not performed. Please use `AnalyzeInterSpecificity` to generate that. ")
		}
		return(object@tg.spgenes)
	}
)


#' Using Gene Pair Split
#' 
#' The \code{setGenePairSplit} function is to \bold{set} the gene pair split, 
#' which is either one charater or string, and used to split 2 gene partners in one gene pair. 
#'
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
		use.clusters <- unique(as.character(object@fgenes$cluster))
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
#' @inheritParams InsideObjectInterCell
#' @param ... Parameters passed to other methods.
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
# Database Selection Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Select Subset of Gene Pair Database
#'
#' This function is to select subset of gene pair databases, and options on evidence sources, 
#' confidence level, action properties (mode and effect) will be used. The result from different evidence sources
#' will be \bold{union} in default settings. The result from action properties could be either intersection result or union.
#' 
#' @param object A \code{InterCellDBPack} or \code{InterCell} object. \code{SelectDBSubset.default} uses
#'  \code{InterCellDBPack} as input.
#' @param combined.score.range The combined score from all evidence sources, which works on the whole database.
#' @param use.exp It adds the data whose evidence is experimentally validatd. 
#' @param exp.score.range It controls the score range when selecting subset by experimentally validatd evidence, i.e. when \code{use.exp = TRUE}. 
#'  The score should be 2 numbers within 1~1000.
#' @param use.know It adds the data whose evidence is pathway curated. 
#' @param know.score.range It controls the score range when selecting subset by pathway curated evidence, i.e. when \code{use.know = TRUE}. 
#'  The score should be 2 numbers within 1~1000.
#' @param use.pred It adds the data whose evidence is predicted.
#' @param pred.score.range It controls the score range when selecting subset by predicted evidence, i.e. when \code{use.pred = TRUE}. 
#'  The score should be 2 numbers within 1~1000.
#' @param sel.physical It selects the subset of gene pairs whose corresponding protein pairs are physical associated. This parameter is 
#'  identical to set \code{sel.action.mode = "binding"} for now due to the database limitation. It may change in future. 
#' @param sel.action.mode Selection by action mode. "ALL" means not use this to select subset.
#'  Other options will be directly select gene pair in that action mode. Supported options are listed in \code{kpred.action.mode}.
#' @param sel.action.effect Selection by action effect. "ALL" means not use this to select subset. Other 
#'  options will be directly select gene pair in that action effect. Supported options are listed in \code{kpred.action.effect}.
#' @param sel.action.merge.option Either 'intersect' or 'union'. The option for merging the result from selection on action mode and action effect.
#' @param slim.along.with.pairs This decides whether to select the corresponding subset of action pair database 
#'  after selecting subset of gene pair database.
#'
#' @details
#' The 3 evidence channels (\code{use.exp}, \code{use.know}, \code{use.pred}) are in identical priority. 
#' For one gene pair, it could have scores from 3 evidence channels at the same time. In default setting, the 
#' result of this function will be \bold{union} of results from 3 channels. If union is not satisfied, the function 
#' \code{\link{MergeDBSubset}} will help.
#
#' The score is within 1~1000. Score (>700) would be consider as high confidence. Score (>400) would be consider over medium confidence.
#' Score (<=400) would be consider low confidence, while Score (<=150) would be the lowest.
#'
#' @return A \code{InterCellDBPack} or \code{InterCell} object, which is the same as given parameter \code{object}.
#'
#' @rdname SelectDBSubset
#' @order 4
#'
#'
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
	sel.action.mode = "ALL",
	sel.action.effect = "ALL",
	sel.action.merge.option = "intersect",
	slim.along.with.pairs = TRUE
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
		not.valid.action.mode <- setdiff(sel.action.mode, kpred.action.mode)
		if (length(not.valid.action.mode) > 0) {
			warning("Given undefined action mode: ", paste0(not.valid.action.mode, collapse = ", ", ". "))
		}
		sel.action.mode <- intersect(sel.action.mode, kpred.action.mode)
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
	# get exp part
	if (use.exp == TRUE) {
		retDB.list <- c(retDB.list, list(intersect(which(this.pairs.db$inter.Experiments.Score >= exp.score.range[1]),
			which(this.pairs.db$inter.Experiments.Score <= exp.score.range[2]))))
	}
	if (use.know == TRUE) {
		retDB.list <- c(retDB.list, list(intersect(which(this.pairs.db$inter.Database.Score >= know.score.range[1]),
			which(this.pairs.db$inter.Database.Score <= know.score.range[2]))))
	}
	if (use.pred == TRUE) {
		retDB.list <- c(retDB.list, list(intersect(which(this.pairs.db$inter.Predicted.Score >= pred.score.range[1]),
			which(this.pairs.db$inter.Predicted.Score <= pred.score.range[2]))))
	}
	# collect result from score selection
	this.pairs.db <- this.pairs.db[Reduce(union, retDB.list), ]
	#this.pairs.db <- DoPartUnique(this.pairs.db, match(c("inter.GeneName.A", "inter.GeneName.B"), colnames(this.pairs.db)))

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
		# add direct way ones (conv)
		tmp.sel.pairs <- paste(object@pairs.db[, "inter.GeneID.A"], object@pairs.db[, "inter.GeneID.B"], sep = use.cut.symbol.slim)
		# add reverse direction ones (rev)
		tmp.sel.pairs.conv <- c(tmp.sel.pairs , paste(object@pairs.db[, "inter.GeneID.B"], object@pairs.db[, "inter.GeneID.A"], sep = use.cut.symbol.slim))
		# slim actions
		tmp.sel.actions.p <- paste(object@actions.db[, "inter.GeneID.A"], object@actions.db[, "inter.GeneID.B"], sep = use.cut.symbol.slim)
		object@actions.db <- object@actions.db[which(tmp.sel.actions.p %in% tmp.sel.pairs), ]
	}

	return(object)
}

#' @param ... Parameters passed to function \code{SelectDBSubset.default}.
#'
#' @rdname SelectDBSubset
#' @order 1
#' @export
#'
setGeneric(name = "SelectDBSubset", def = function(object, ...) {
	standardGeneric("SelectDBSubset")
	}
)

#' @rdname SelectDBSubset
#' @order 2
#' @export
#'
setMethod(
	f = "SelectDBSubset",
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		SelectDBSubset.default(object, ...)
	}
)

#' @rdname SelectDBSubset
#' @order 3
#' @export
#'
setMethod(
	f = "SelectDBSubset",
	signature = "InterCell",
	definition = function(object, ...) {
		setRefDatabase(object, SelectDBSubset.default(object@database, ...))
	}
)


#' Merge Database
#'
#' This function is to merge 2 databases to get user-desired database.
#'
#' @param db.object.1 A \code{InterCellDBPack} object to be merged.
#' @param db.object.2 Another \code{InterCellDBPack} object to be merged.
#' @param merge.option The supported options are 'intersect' and 'union'.
#'
#' @return A \code{InterCellDBPack} object.
#'
#' @export
#'
MergeDBSubset <- function(
	db.object.1,
	db.object.2,
	merge.option = "intersect"
) {
	if (class(db.object.1) != "InterCellDBPack" || class(db.object.2) != "InterCellDBPack") {
		ret.mg <- character()
		if (class(db.object.1) != "InterCellDBPack") ret.mg <- c(ret.mg, "db.object.1")
		if (class(db.object.2) != "InterCellDBPack") ret.mg <- c(ret.mg, "db.object.2")
		stop("Given database subset in parameter ", paste0(ret.mg, collapse = ", ") , " not in `InterCellDBPack-class`. ",
			 "Please use either `CreateDBPackObject` or `SelectDBSubset` to generate.")
	}

	if (length(merge.option) > 1) {
		merge.option <- merge.option[1]
	}
	merge.option <- CheckParamStd(merge.option, c("intersect", "union"), "parameter `merge.option`", stop.on.zero = TRUE)

	# before merging, check species
	if (!identical(db.object.1@species, db.object.2@species)) {
		stop("Merging database subset from different species ", paste(db.object.1@species, db.object.2@species, sep = ", "), ". ")
	} else {
		if (merge.option == "union") {
			# genes.db 
			# $gene.ncbi.db
			addon.ugene.ncbi.db <- setdiff(rownames(db.object.2@genes.db$gene.ncbi.db), rownames(db.object.1@genes.db$gene.ncbi.db))
			db.object.1@genes.db$gene.ncbi.db <- rbind(db.object.1@genes.db$gene.ncbi.db, db.object.2@genes.db$gene.ncbi.db[which(rownames(db.object.2@genes.db$gene.ncbi.db) %in% addon.ugene.ncbi.db), ])
			# $gene.synonyms.db
			addon.ugene.synonyms.db <- setdiff(rownames(db.object.2@genes.db$gene.synonyms.db), rownames(db.object.1@genes.db$gene.synonyms.db))
			db.object.1@genes.db$gene.synonyms.db <- rbind(db.object.1@genes.db$gene.synonyms.db, db.object.2@genes.db$gene.synonyms.db[which(rownames(db.object.2@genes.db$gene.synonyms.db) %in% addon.ugene.synonyms.db), ])
			# $gene.dup.synonyms.db
			db.object.1@genes.db$gene.dup.synonyms.db <- unique(rbind(db.object.1@genes.db$gene.dup.synonyms.db, db.object.2@genes.db$gene.dup.synonyms.db))

			# pairs.db
			addon.pairs.db <- setdiff(rownames(db.object.2@pairs.db), rownames(db.object.1@pairs.db))
			db.object.1@pairs.db <- rbind(db.object.1@pairs.db, db.object.2@pairs.db[which(rownames(db.object.2@pairs.db) %in% addon.pairs.db), ])
			# actions.db
			addon.actions.db <- setdiff(rownames(db.object.2@actions.db), rownames(db.object.1@actions.db))
			db.object.1@actions.db <- rbind(db.object.1@actions.db, db.object.2@actions.db[which(rownames(db.object.2@actions.db) %in% addon.actions.db), ])
			# anno.location.db
			addon.anno.location.db <- setdiff(rownames(db.object.2@anno.location.db), rownames(db.object.1@anno.location.db))
			db.object.1@anno.location.db <- rbind(db.object.1@anno.location.db, db.object.2@anno.location.db[which(rownames(db.object.2@anno.location.db) %in% addon.anno.location.db), ])
			# anno.type.db 
			addon.anno.type.db <- setdiff(rownames(db.object.2@anno.type.db), rownames(db.object.1@anno.type.db))
			db.object.1@anno.type.db <- rbind(db.object.1@anno.type.db, db.object.2@anno.type.db[which(rownames(db.object.2@anno.type.db) %in% addon.anno.type.db), ])
			# go.ref.db 
			addon.go.ref.db <- setdiff(rownames(db.object.2@go.ref.db), rownames(db.object.1@go.ref.db))
			db.object.1@go.ref.db <- rbind(db.object.1@go.ref.db, db.object.2@go.ref.db[which(rownames(db.object.2@go.ref.db) %in% addon.go.ref.db), ])
		}
		if (merge.option == "intersect") {
			# genes.db 
			# $gene.ncbi.db
			addon.ugene.ncbi.db <- intersect(rownames(db.object.2@genes.db$gene.ncbi.db), rownames(db.object.1@genes.db$gene.ncbi.db))
			db.object.1@genes.db$gene.ncbi.db <- db.object.2@genes.db$gene.ncbi.db[which(rownames(db.object.2@genes.db$gene.ncbi.db) %in% addon.ugene.ncbi.db), ]
			# $gene.synonyms.db
			addon.ugene.synonyms.db <- intersect(rownames(db.object.2@genes.db$gene.synonyms.db), rownames(db.object.1@genes.db$gene.synonyms.db))
			db.object.1@genes.db$gene.synonyms.db <- db.object.2@genes.db$gene.synonyms.db[which(rownames(db.object.2@genes.db$gene.synonyms.db) %in% addon.ugene.synonyms.db), ]
			# $gene.dup.synonyms.db
			addon.ugene.dup.synonyms.db <- intersect(rownames(db.object.2@genes.db$gene.dup.synonyms.db), rownames(db.object.1@genes.db$gene.dup.synonyms.db))
			db.object.1@genes.db$gene.dup.synonyms.db <- db.object.2@genes.db$gene.dup.synonyms.db[which(rownames(db.object.2@genes.db$gene.dup.synonyms.db) %in% addon.ugene.dup.synonyms.db), ]
	
			# pairs.db
			addon.pairs.db <- intersect(rownames(db.object.2@pairs.db), rownames(db.object.1@pairs.db))
			db.object.1@pairs.db <- db.object.2@pairs.db[which(rownames(db.object.2@pairs.db) %in% addon.pairs.db), ]
			# actions.db
			addon.actions.db <- intersect(rownames(db.object.2@actions.db), rownames(db.object.1@actions.db))
			db.object.1@actions.db <- db.object.2@actions.db[which(rownames(db.object.2@actions.db) %in% addon.actions.db), ]
			# anno.location.db
			addon.anno.location.db <- intersect(rownames(db.object.2@anno.location.db), rownames(db.object.1@anno.location.db))
			db.object.1@anno.location.db <- db.object.2@anno.location.db[which(rownames(db.object.2@anno.location.db) %in% addon.anno.location.db), ]
			# anno.type.db 
			addon.anno.type.db <- intersect(rownames(db.object.2@anno.type.db), rownames(db.object.1@anno.type.db))
			db.object.1@anno.type.db <- db.object.2@anno.type.db[which(rownames(db.object.2@anno.type.db) %in% addon.anno.type.db), ]
			# go.ref.db 
			addon.go.ref.db <- intersect(rownames(db.object.2@go.ref.db), rownames(db.object.1@go.ref.db))
			db.object.1@go.ref.db <- db.object.2@go.ref.db[which(rownames(db.object.2@go.ref.db) %in% addon.go.ref.db), ]
		}
	}

	return(db.object.1)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Other Related Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Remap Gene Symbols to Their Authorized Genenames
#'
#' @description
#' Use species-specific reference database to remap genes to their authorized genenames.
#'
#' @param markers.all Data.frame. Feature genes that are generated from \code{Seurat::FindAllMarkers()} or 
#' similar functions in other packages.
#' @inheritParams InsideParam.species
#' @param if.used.inside Logic. If used inside, some process will not run.
#' @param final.dup.rm Logic. If set TRUE, duplicated genes of each cluster in final mapping result will be 
#'  removed. If set FALSE, those genes will be applied with function \code{make.unique}.
#'
#' @export
#'
DataPrep.RemapClustersMarkers  <- function(
	markers.all,
	species,
	if.used.inside = FALSE,
	final.dup.rm = TRUE
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
	markers.all <- rbind(markers.raw.match, markers.raw.unmatch)

	# after remapping, genes get to be duplicate with existing ones, 
	# re-check if mapping result has duplicate genes
	if (final.dup.rm == TRUE) {
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
	} else {
		fcheck.result <- lapply(unique(markers.all$cluster), ref.markers = markers.all,
			function(x, ref.markers) {
				this.c.markers <- ref.markers[which(ref.markers$cluster == x), ]
				# in default: make unique of those genes
				this.c.markers$gene <- make.unique(this.c.markers$gene)
				list(dup.genes = character(), markers = this.c.markers)
			})
		names(fcheck.result) <- unique(markers.all$cluster)
	}
	

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



#' Add Expression Data in InterCell Object
#'
#' This function adds expression data to InterCell Object. It is mostly internally used.
#'
#' @param fgenes feature gene table, stored in slot \code{fgenes} of \code{InterCell} object.
#' @param exprs.data normalized count matrix
#' @param force.overwrite If column 'Exprs' already exists, set this parameter \code{TRUE} to overwrite that.
#'
#' @export
#'
DataPrep.AddExprs <- function(
	fgenes,
	exprs.data,
	force.overwrite = FALSE
) {
	if ("Exprs" %in% colnames(fgenes) && !force.overwrite) {
		stop("Column 'Exprs' exists. Use 'force.overwrite = TRUE' to force overwrite existed one.")
	}

	# handle with format of exprs.data
	if ("Seurat" %in% class(exprs.data)) {
		stop("Please use <SeuratObj>[[<assay>]]@<data>, for example SeuratObj[['RNA']]@data, to give the expression data.")
		#exprs.data <- AverageExpression(exprs.data, assays = "RNA", slot = "counts")[[1]]
		#exprs.data <- as.matrix(exprs.data)
	} else {  # transform to matrix
		if (is.null(ncol(exprs.data)) || ncol(exprs.data) == 0 ||
			is.null(nrow(exprs.data)) || nrow(exprs.data) == 0) {
			stop("Please provided valid matrix or data.frame to use.")
		}
		exprs.data <- as.matrix(exprs.data)
	}

	# give warning when duplicated genes are given
	if (anyDuplicated(rownames(exprs.data))) {
		warning("Duplicated gene names in expression data detected. Unexpected result may be generated.")
	}

	# check if cluster are matched with current data
	if (!all(unique(fgenes$cluster) %in% colnames(exprs.data))) {
		stop("Given data lacks some required clusters: ", 
			paste0(setdiff(unique(fgenes$cluster), colnames(exprs.data)), collapse = ", "),
			"!")
	}

	# add expression data in every cluster
	fgenes <- bind_rows(lapply(unique(fgenes$cluster), 
		fgenes = fgenes, exprs.data = exprs.data,
		function(x, fgenes, exprs.data) {
			tmp.exprs <- exprs.data[, which(colnames(exprs.data) == x), drop = FALSE]
			tmp.exprs <- rowSums(tmp.exprs) / ncol(tmp.exprs)
			fgenes <- subset(fgenes, cluster == x)
			inds.exprs <- match(fgenes$gene, names(tmp.exprs))
			fgenes[, "Exprs"] <- as.numeric(tmp.exprs[inds.exprs])
			return(fgenes)
		}
	))

	return(fgenes)
}



#' Generate Gene Pairs in Standard Format
#' 
#' @description
#' This function generates gene pairs in standard format(in data frame), 
#' and gets these pairs easier to be compared with others.
#'
#' @param VEinfos standard storage for one interaction.
#'
#' @details
#' The standard format in this package is that gene pairs are maintained in data.frame, and the 2 genes 
#' participated in each gene pair are recorded in columns named "inter.GeneName.A" and "inter.GeneName.B".
#'
#' @importFrom dplyr left_join
#'
DataPrep.GenStdGenePairs.from.VEinfos <- function(
  VEinfos,
  musthave.colnames
) {
	vertices.infos <- VEinfos$vertices.infos
	edges.infos <- VEinfos$edges.infos
	# pre-check
	if (!all(kmusthave.colnames %in% musthave.colnames)) {
		stop("Provided column names are not matched with InterCellDB requirement!")
	}
	musthave.colnames <- setdiff(musthave.colnames, c("gene", "cluster"))
	#
	tmp.res <- left_join(edges.infos[, c("from", "to")], vertices.infos[, c(c("UID", "ClusterName", "GeneName"), musthave.colnames)], by = c("from" = "UID"))
	colnames(tmp.res)[c(ncol(tmp.res) - (1+length(musthave.colnames)):0)] <- paste("inter", c(c("Cluster", "GeneName"), musthave.colnames), "A", sep = ".")
	tmp.res <- left_join(tmp.res, vertices.infos[, c(c("UID", "ClusterName", "GeneName"), musthave.colnames)], by = c("to" = "UID"))
	colnames(tmp.res)[c(ncol(tmp.res) - (1+length(musthave.colnames)):0)] <- paste("inter", c(c("Cluster", "GeneName"), musthave.colnames), "B", sep = ".")
	# form std data.frame
	align.colnames <- paste("inter", rep(c(c("GeneName"), musthave.colnames, c("Cluster")), each = 2), c("A", "B"), sep = ".")
	tmp.res <- tmp.res[, match(align.colnames, colnames(tmp.res))]

	# match cluster
	# get conv ones
	std.res.conv <- tmp.res[intersect(which(tmp.res$inter.Cluster.A == VEinfos$cluster.name.A), which(tmp.res$inter.Cluster.B == VEinfos$cluster.name.B)), ]
	# get rev ones
	std.res.rev <- tmp.res[intersect(which(tmp.res$inter.Cluster.A == VEinfos$cluster.name.B), which(tmp.res$inter.Cluster.B == VEinfos$cluster.name.A)), ]
	std.res.rev <- std.res.rev[, ReverseOddEvenCols(length(align.colnames))]  # reverse all paired columns
	colnames(std.res.rev) <- colnames(std.res.conv)
	# get the result
	std.res.all <- rbind(std.res.conv, std.res.rev)
	std.res.all <- DoPartUnique(std.res.all, 1:2)

	return(std.res.all)
}


#' List cell clusters
#'
#' This function is to list all cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#'
#' @export
#'
ListAllClusters <- function(
	object
) {
	return(unique(object@fgenes$cluster))
}



#' List Options for Gene Selection
#' 
#' This group of functions are served for giving options on selecting gene subset. It ranges from 
#' subcellular locations, molecular functions and GO terms.
#'
#' @param object Allowed object should be in either class \code{InterCell} or class \code{InterCellDBPack}.
#' @param ... Parameters passed to corresponding function with suffix 'default', like use \code{ListAllGeneLocation} 
#'  and check parameters in \code{ListAllGeneLocation.default}.
#'
#' @return Character. The options.
#'
#' @name ListGeneSelectionProperty
#' @rdname ListGeneSelectionProperty
#'
#'
NULL


#' @rdname ListGeneSelectionProperty
#' @export
#'
ListAllGeneLocation.default <- function(
	object
) {
	return(unique(object@anno.location.db$GO.Term.target))
}

#' @rdname ListGeneSelectionProperty
#' @order 1
#' @export
#'
setGeneric(name = "ListAllGeneLocation", def = function(object, ...) {
	standardGeneric("ListAllGeneLocation")
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneLocation",
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		ListAllGeneLocation.default(object, ...)
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneLocation",
	signature = "InterCell",
	definition = function(object, ...) {
		ListAllGeneLocation.default(object@database, ...)
	}
)


#' @rdname ListGeneSelectionProperty
#' @export
#'
ListAllGeneType.default <- function(
	object
) {
	return(unique(object@anno.type.db$Keyword.Name))
}

#' @rdname ListGeneSelectionProperty
#' @order 2
#' @export
#'
setGeneric(name = "ListAllGeneType", def = function(object, ...) {
	standardGeneric("ListAllGeneType")
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneType",
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		ListAllGeneType.default(object, ...)
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneType",
	signature = "InterCell",
	definition = function(object, ...) {
		ListAllGeneType.default(object@database, ...)
	}
)

#' @rdname ListGeneSelectionProperty
#' @export
#'
ListAllGeneMergeType.default <- function(
	object
) {
	return(unique(object@accessory.db$merge.type.list[, "merged.molecular.function"]))
}
#' @rdname ListGeneSelectionProperty
#' @order 3
#' @export
#'
setGeneric(name = "ListAllGeneMergeType", def = function(object, ...) {
	standardGeneric("ListAllGeneMergeType")
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneMergeType",
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		ListAllGeneMergeType.default(object, ...)
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneMergeType",
	signature = "InterCell",
	definition = function(object, ...) {
		ListAllGeneMergeType.default(object@database, ...)
	}
)


#' @param n.output It gives the first N GO terms in the order of AlphaBet.
#'
#' @rdname ListGeneSelectionProperty
#' @export
#'
ListAllGeneGOTerm.default <- function(
	object,
	n.output = +Inf
) {
	go.terms <- unique(object@go.ref.db$GO_term)
	go.terms <- go.terms[order(go.terms, decreasing = FALSE)]
	if (n.output < 1) {
		n.output <- 1
	}
	if (length(go.terms) > n.output) {
		go.terms <- go.terms[seq_len(n.output)]
	}

	return(go.terms)
}

#' @rdname ListGeneSelectionProperty
#' @order 4
#' @export
#'
setGeneric(name = "ListAllGeneGOTerm", def = function(object, ...) {
	standardGeneric("ListAllGeneGOTerm")
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneGOTerm",
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		ListAllGeneGOTerm.default(object, ...)
	}
)
#' @rdname ListGeneSelectionProperty
#' @export
#'
setMethod(
	f = "ListAllGeneGOTerm",
	signature = "InterCell",
	definition = function(object, ...) {
		ListAllGeneGOTerm.default(object@database, ...)
	}
)



#' Replace Cluster name
#'
#' This is function is to change the names of cell clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param cluster.names.current The currently used name of cell clusters.
#' @param cluster.names.replace The new used name of cell clusters.
#'
#' @return A \code{InterCell} object.
#'
#' @examples
#' \dontrun{
#'   ReplaceClusterName(object, "Macrophage", "Myeloid")
#' }
#'
#' @export
#'
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


#' Fetch Genes of Interest
#'
#' This function is to fetch genes of interest by selecting on subcellular locations, 
#' molecular functions and GO terms.
#'
#' @param object A \code{InterCell} object or \code{InterCellDBPack} object. \code{FetchGeneOI.default}
#'  gets \code{InterCellDBPack} as input.
#' @param sel.location Use subcellular location of gene product to select gene, and options are listed 
#'  in \code{\link{ListAllGeneLocation}}.
#' @param sel.location.score The score of corresponding subcellular location, range from 1 to 5. Consider \code{score = c(4, 5)} 
#'  as high and common usage.
#' @param sel.type Use molecular function of gene product to select gene, and options are listed in 
#'  \code{\link{ListAllGeneType}}.
#' @param sel.merge.type The merged types. Comparing to \code{sel.type}, it has less options, which are given 
#'  in \code{\link{ListAllGeneMergeType}}. See details for help.
#' @param ret.with.property It decides whether the return values are attached with gene properties(location, type, etc).
#' @param sel.go.terms Use GO terms to select gene, and supported options are listed in \code{\link{ListAllGeneGOTerm}}.
#' @param go.use.relative Decide if the go terms in the related GO term tree are used. See details for help.
#' @param go.relative.option Decide which relation to the selected GO term is used. Options are 
#'  'ancestor', 'parents', 'offspring', 'children'. See details for help.
#'
#' @details
#' The parameter \code{sel.merge.type} is the summary for \code{sel.type}. The options in \code{sel.type} are originally 
#' generated directly from Uniprot, which comprises over 100 types. For the convenience of usage, we summarize those types and gather
#' them to 16 merged types, which comprise the common used types: 'Receptor', 'Cytokine', 'Growth Factor', etc.
#'
#' GO terms have 3 basic words: 'cellular_component', 'molecular_function', 'biological_process', and all other GO terms are the offspring of 
#' one of them. As a result, GO terms form 3 family tree. Use parameter \code{go.use.relative}, it can extend one given GO term to all its 
#' related GO terms. There are 4 pre-defined options for selecting specific relative group, which are 'ancestor', 'parents', 'offspring', 'children'.
#' The 'parents' and 'children' are selecting the most close GO terms. For example, 'GO:0006955-immune response' is the nearest level above 
#' 'GO:0002250-adaptive immune response', and in turn, 'adaptive immune response' is the children of 'immune response'. 
#' The 'ancestor' and 'offspring' goes further than 'parents' and 'children'. The 'ancestor' iteratively searches for the 'parents'. 
#' Conversely, the 'offspring' iteratively searches for the 'children'. For example, 'immune response' is the parents of 'adaptive immune response', and 
#' 'adaptive immune response' is the parents of 'adaptive immune effector response'. Then, 'immune response' is the 'ancestor' of 'adaptive immune effector response'.
#' The 'offspring' goes the same way by propagating 'children'. 
#'
#' @return Character. The selected genes.
#'
#' @rdname FetchGeneOI
#' @order 4
#' @export
#'
FetchGeneOI.default <- function(
	object,
	sel.location = NULL,
	sel.location.score = c(1:5),
	sel.type = NULL,
	sel.merge.type = NULL,  # merged type only have 16 options. Prioritize upon sel.type
	ret.with.property = TRUE,  
	# whether the selection should be carried in the result, 
	# for analysis usage it's TRUE, but only to get some genes, it should be false.
	# GO terms are not included by this functions
	sel.go.terms = NULL,  # ID and Term are supported
	go.use.relative = TRUE, 
	go.relative.option = "offspring"
) {
	## input parameter process
	# check location
	if (!is.null(sel.location)) {
		avb.opt.location <- unique(object@anno.location.db$GO.Term.target)
		not.valid.location <- setdiff(sel.location, avb.opt.location)
		if (length(not.valid.location) > 0) {
			warning("Given undefined location: ", paste0(not.valid.location, collapse = ", ", ". "))
		}
		sel.location <- intersect(sel.location, avb.opt.location)
		if (length(sel.location) == 0) {
			sel.location <- NULL
		}
	}
	if ((length(sel.location.score) > 1 && !is.integer(sel.location.score)) || 
		(length(sel.location.score) == 1 && !is.numeric(sel.location.score))) {
		stop("Location score ranges from 1 to 5, and only those 5 integers are supported!")
	}
	if (!is.null(sel.location.score)) {
		avb.location.score <- c(1:5)
		not.valid.location.score <- setdiff(sel.location.score, avb.location.score)
		if (length(not.valid.location.score) > 0) {
			warning("Given invalid location score: ", paste0(not.valid.location.score, collapse = ", ", ". "))
		}
		sel.location.score <- intersect(sel.location.score, avb.location.score)
	}

	# check type
	if (!is.null(sel.type)) {
		avb.opt.type <- unique(object@anno.type.db$Keyword.Name)
		not.valid.type <- setdiff(sel.type, avb.opt.type)
		if (length(not.valid.type) > 0) {
			warning("Given undefined type: ", paste0(not.valid.type, collapse = ", ", ". "))
		}
		sel.type <- intersect(sel.type, avb.opt.type)
		if (length(sel.type) == 0) {
			sel.type <- NULL
		}
	}
	# check merged type (if given, then override type) and overwrite type
	if (!is.null(sel.merge.type)) {
		if (!is.null(sel.type)) {
			warning("Select genes using both parameter 'sel.type' & 'sel.merge.type'. Only the options in 'sel.merge.type' are used!")
		}
		avb.opt.mg.type <- unique(object@accessory.db$merge.type.list$merged.molecular.function)
		sel.merge.type <- CheckParamStd(sel.merge.type, avb.opt.mg.type, "merged type", stop.on.zero = FALSE)
		if (length(sel.merge.type) == 0) {
			sel.type <- NULL
		} else {  # overwrite sel.type
			ref.db <- object@accessory.db$merge.type.list
			sel.type <- ref.db[which(ref.db$merged.molecular.function %in% sel.merge.type), "orig.uniprot.keywords"]
		}
	}

	## save all above limitations
	property.saved <- list(
		location = sel.location,
		location.score = sel.location.score,
		type = sel.type
	)

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
	if (!is.null(sel.location)) {
		gene.oi.from.locs <- object@anno.location.db[intersect(which(object@anno.location.db$GO.Term.target %in% sel.location), 
			which(object@anno.location.db$score %in% sel.location.score)), "Gene.name"]
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.locs)
	}
	if (!is.null(sel.type)) {
		gene.oi.from.type <- object@anno.type.db[which(object@anno.type.db$Keyword.Name %in% sel.type), "Gene.name"]
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.type)
	}
	if (!is.null(sel.go.terms)) {
		gene.oi.from.go <- Tool.FindGenesFromGO(sel.go.terms, object@genes.db, object@go.ref.db, 
			go.use.relative = go.use.relative, go.relative.option = go.relative.option)
		gene.oi.from.go <- as.character(unlist(gene.oi.from.go))
		ret.gene.oi <- inside.set.gene.oi(ret.gene.oi, gene.oi.from.go)
	}
	# get unique result
	ret.gene.oi <- unique(ret.gene.oi)

	print(paste("Fetch", length(ret.gene.oi), "genes of interest."))
	if (ret.with.property == TRUE) {
		return(list(genes = ret.gene.oi, property = property.saved))	
	} else {
		return(ret.gene.oi)
	}
}


#' @param ... Parameters passed to function \code{FetchGeneOI.default}.
#'
#' @rdname FetchGeneOI
#' @order 1
#' @export 
#'
setGeneric(name = "FetchGeneOI", def = function(object, ...) {
		standardGeneric("FetchGeneOI")
	}
)

#' @rdname FetchGeneOI
#' @order 2
#' @export 
#'
setMethod(
	f = "FetchGeneOI", 
	signature = "InterCellDBPack",
	definition = function(object, ...) {
		return(FetchGeneOI.default(object, ...))
	}
)

#' @rdname FetchGeneOI
#' @order 3
#' @export 
#'
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
#' @param gene.pairs.table 2-column table, and each column records one list of genes.
#' @inheritParams InsideParam.species
#' @param extend.reverse Decide whether to extend pairs with the reverse pairs. It is useful
#'  when analytic process requires cell-cluster-aligned gene pairs, e.g. \code{\link{AnalyzeInterInFullView}}.
#'
#' @details
#' The formatting process will remap the genes in 1st column from input to 'inter.*.A' columns in result.
#' The result of 2nd column from input will be put in 'inter.*.B' columns.
#'
#' Explanation on \code{extend.reverse}:
#' For example, given pair C3~C3ar1, if set \code{extend.reverse = TRUE}, then 
#' both C3~C3ar1 and C3ar1~C3 will be generated in result.
#'
#' @return A list with \code{$result} storing the formatted gene pairs.
#'
#' @export
#'
FormatCustomGenePairs <- function(
	gene.pairs.table,
	species,
	extend.reverse = FALSE
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
		tmp.genes <- as.character(gene.pairs.table[, i])
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

	# get genename corresponding IDs
	genes.ref.db <- switch(species, 
		"human" = genes.human.ref.db$gene.ncbi.db,
		"mouse" = genes.mouse.ref.db$gene.ncbi.db)
	result.id.A <- genes.ref.db[match(result.list[[1]], genes.ref.db$Symbol_from_nomenclature_authority), "GeneID"]
	result.id.B <- genes.ref.db[match(result.list[[2]], genes.ref.db$Symbol_from_nomenclature_authority), "GeneID"]

	result.IT <- data.frame(inter.GeneID.A = result.id.A,
		inter.GeneID.B = result.id.B,
		inter.GeneName.A = result.list[[1]],
		inter.GeneName.B = result.list[[2]],
		stringsAsFactors = FALSE)

	if (extend.reverse == TRUE) {
		tmp.rev <- result.IT[, ReverseOddEvenCols(4)]
		colnames(tmp.rev) <- colnames(result.IT)
		result.IT <- rbind(result.IT, tmp.rev)
	}

	# return
	list(result = result.IT, match.status = apx.list)
}



