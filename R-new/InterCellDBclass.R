
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pre-defined default variables and functions of Class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#' Predefined action effect
#'
#' @description
#' This is extracted from databases depicting action effect of actions for functional usages.
#'
kpred.action.effect <- c("positive", "negative", "unspecified", "undirected")



# [inside usage]
# This function evaluates the power of gene-interact pairs by the formula:
#   SUM(abs(LogFC[geneA] * LogFC[geneB]))
FullView.Evaluation.func.default <- function(
  pairs.ref,
  colname.eval = c("inter.LogFC.A", "inter.LogFC.B")
) {
  eval.res <- 0.0
  if (nrow(pairs.ref) > 0) {
    eval.res <- sum(abs(pairs.ref[, colname.eval[1]] * pairs.ref[, colname.eval[2]]))
  }
  eval.res
}



#' Formula on calculating LogFC
#'
#' @description
#' This is the default formula on calculating LogFC values, and is 
#' belong to \code{GetResult.PlotOnepairClusters.GeneCrosstalk()}.
#'
#' @param data.f vector. The operated data.
#' @param data.b vector. The operated data with its length the same as the 1st parameter.
#'
#' @details
#' This is the default formula for calculation. Users can define their own formula,
#' as well as define their own evalution parameter besides LogFC(by changing the 
#' param - \code{colnames.to.cmp} of \code{GetResult.PlotOnepairClusters.GeneCrosstalk()}).
#'
#'
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



#' Formula on calculating p_val_adj
#'
#' @description
#' This is the default formula on calculating p_val_adj values, and is 
#' belong to \code{GetResult.PlotOnepairClusters.GeneCrosstalk()}.
#'
#' @param data.f vector. The operated data.
#' @param data.b vector. The operated data with its length the same as the 1st parameter.
#'
#' @details
#' This is the default formula for calculation. Users can define their own formula,
#' as well as define their own evalution parameter besides p_val_adj(by changing the 
#' param - \code{colnames.to.cmp} of \code{GetResult.PlotOnepairClusters.GeneCrosstalk()}).
#'
#'
#'
TgView.formula.onPVal.default <- function(
	data.f, 
	data.b,
	pval.log.max = 1000
) {
	default.max.replace <- pval.log.max  # use 1000 as default maximum(e-999 are usual lowest limit), when all log(values) are infinite.
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# InterCell Object definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
		accessory.db = "list",  # a list of accessory database
		misc = "list"  # other things, default NULL
	)
)

InterCell <- setClass(
	Class = "InterCell",
	slots = c(
		# const part
		fgenes = "data.frame", 
		database = "InterCellDBPack", 
		formulae = "list",  # formulae on fullview, tg-logFC, Pval
		# analysis result part
		inter.fullview = "list",  # result of fullview of all interactions
		tg.action.pairs = "list",  # gene pairs and their actions in target 2-cell group
		tg.veinfo = "list",  # transformed to vertices and edges
		tg.spgenes = "list",
		# other: global variables, etc
		pred.action = "list",  # kpred.mode, kpred.action.effect
		tool.vars = "list",  # kGeneSplit, kClusterSplit
		misc = "list"  # put gene-align-unmatched results, etc in it
	)
	# prototype # validity
)



validInterCellDBPackObject <- function(object) {
	# check species
	allowed.species <- c("human", "mouse")
	if (length(object@species) != 1 || sum(allowed.species %in% object@species) == 0) {
		return(paste("InterCellDB only supports species in human and mouse."))
	}
	# [TODO] check more things
	TRUE
}
setValidity("InterCellDBPack", validInterCellDBPackObject)

# only the most un-modifiable slots should be checked, 
# and all checked ones should set its initial value in the first time, or get errors
validInterCellObject <- function(object) {
	# check pred.action
	if (length(object@pred.action$action.mode) == 0 || !all(object@pred.action$action.mode %in% kpred.mode)) {
		return(paste0("Using undefined action mode: ", paste0(setdiff(object@pred.action$action.mode, kpred.mode), collapse = ", ")))
	}
	if (length(object@pred.action$action.effect) == 0 || !all(object@pred.action$action.effect %in% kpred.action.effect)) {
		return(paste0("Using undefined action effect: ", paste0(setdiff(object@pred.action$action.effect, kpred.action.effect), collapse = ", ")))
	}

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
	# [TODO] check more things

	validInterCellDBPackObject(object@database)
	TRUE
}
setValidity("InterCell", validInterCellObject)





# create InterCellDBPack Object using new(***), and default values are settled
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

# create InterCell Object using new(***), and default values are settled
# This is to set default values. [TODO]
setMethod(
	f = "initialize",
	signature = c("InterCell"),
	definition = function(.Object, ...) {
		.Object <- callNextMethod(.Object, ...)
		.Object@formulae <- list(FULLVIEW = FullView.Evaluation.func.default,
			TG.LOGFC = TgView.formula.onLogFC.default,
			TG.PVAL = TgView.formula.onPVal.default)
		#.Object@pred.action <- list(action.mode = kpred.mode, action.effect = kpred.action.effect)
		#.Object@tool.vars <- list(gene.pair.split = "-~-", cluster.split = "~")
		validObject(.Object)
		return(.Object)
	}
)



# [TODO] modified the show function
setMethod(
	f = "show",
	signature = "InterCellDBPack",
	definition = function(object) {
		cat("A InterCellDBPack object, covering ", object@species, ".")
	}
)
# [TODO] modified the show function
setMethod(
	f = "show",
	signature = "InterCell",
	definition = function(object) {
		cat("A InterCell object, get fgenes rows: ", nrow(object@fgenes), ".")
	}
)


# S4 setMethod
# S3 format: function.class like     levels<-.character


# %%%%%%%%%%%%%%
# class methods
# %%%%%%%%%%%%%%

# mostly used internally, not exported for users to directly use
CreateDBPackObject <- function(
	species
) {
	# Check species, only support human and mouse yet
	DB.support.species <- c("human", "mouse")
	if (length(species) != 1 || is.na(match(species[1], DB.support.species))) {
		stop("InterCellDB only supports species to be 'human' and 'mouse' for now.")
	}
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


CreateInterCellObject <- function(
	DEG.table,
	species,
	cluster.split = "~",
	gene.pair.split = "-~-"
) {
	# pre-check part
	if (!all(kmusthave.colnames %in% colnames(DEG.table))) {
		stop("Required columns are not given!\n",
			paste("  Column named ", paste0(setdiff(kmusthave.colnames, colnames(DEG.table)), collapse = ", "), 
			" are not included in given feature genes table.",
			"Please use colnames(<var>)[<index>] <- '<name>' to set proper columns corresponding to those.",
			"Contents in <> should be replaced by user definitions!")
		)
	}

	# set object
	DBPack.obj <- CreateDBPackObject(species)
	# re-align input table with gene reference database in InterCellDB
	DEG.align.res <- DataPrep.RemapClustersMarkers(DEG.table, DBPack.obj@genes.db)
	IT.InterCell.Obj <- new(
		Class = "InterCell",
		fgenes = DEG.align.res$result,
		database = DBPack.obj,
		pred.action = list(action.mode = kpred.mode, action.effect = kpred.action.effect),
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



# %%%%%%%%%%%%%%%%%%%%%%
# accessor function
# %%%%%%%%%%%%%%%%%%%%%%

# fullview limit maximum than others
setGeneric(name = "setFullViewResult", def = function(object, ...) {
	standardGeneric("setFullViewResult")
	}
)
setMethod(
	f = "setFullViewResult",
	signature = "InterCell",
	definition = function(object, new.inter.fullview) {
		if (!is.null(object@inter.fullview) && length(object@inter.fullview) != 0) {
			print("Overwrite existed fullview result [ TODO] name alignment")
		}
		object@inter.fullview <- new.inter.fullview
		return(object)
	}
)
setGeneric(name = "getFullViewResult", def = function(object, ...) {
	standardGeneric("getFullViewResult")
	}
)
setMethod(
	f = "getFullViewResult",
	signature = "InterCell",
	definition = function(object) {
		if (is.null(object@inter.fullview) || length(object@inter.fullview) == 0) {
			stop("FullView not run [TODO] name alignment")
		}
		return(object@inter.fullview)
	}
)


setGeneric(name = "setTgActionPairs", def = function(object, ...) {
	standardGeneric("setTgActionPairs")
	}
)
setMethod(
	f = "setTgActionPairs",
	signature = "InterCell",
	definition = function(object, new.action.pairs) {
		object@tg.action.pairs <- new.action.pairs
		return(object)
	}
)
setGeneric(name = "getTgActionPairs", def = function(object, ...) {
	standardGeneric("getTgActionPairs")
	}
)
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


setGeneric(name = "setTgVEInfo", def = function(object, ...) {
	standardGeneric("setTgVEInfo")
	}
)
setMethod(
	f = "setTgVEInfo",
	signature = "InterCell",
	definition = function(object, new.veinfo) {
		object@tg.veinfo <- new.veinfo
		return(object)
	}
)
setGeneric(name = "getTgVEInfo", def = function(object, ...) {
	standardGeneric("getTgVEInfo")
	}
)
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


setGeneric(name = "setTgSpGenes", def = function(object, ...) {
	standardGeneric("setTgSpGenes")
	}
)
setMethod(
	f = "setTgSpGenes",
	signature = "InterCell",
	definition = function(object, new.spgenes) {
		object@tg.spgenes <- new.spgenes
		return(object)
	}
)
setGeneric(name = "getTgSpGenes", def = function(object, ...) {
	standardGeneric("getTgSpGenes")
	}
)
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


# for @tool.vars$gene.pair.split
setGeneric(name = "setGenePairSplit", def = function(object, ...) {
		standardGeneric("setGenePairSplit")
	}
)
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
setGeneric(name = "getGenePairSplit", def = function(object) {
		standardGeneric("getGenePairSplit")
	}
)
setMethod(
	f = "getGenePairSplit", 
	signature = "InterCell",
	definition = function(object) {
		return(object@tool.vars$gene.pair.split)
	}
)

# for @tool.vars$cluster.split
setGeneric(name = "setClusterSplit", def = function(object, ...) {
		standardGeneric("setClusterSplit")
	}
)
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
setGeneric(name = "getClusterSplit", def = function(object) {
		standardGeneric("getClusterSplit")
	}
)
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
	select.physical = FALSE,
	select.action.mode = "ALL",  # "ALL" consider "Other", so pairs not in actions.db will be reserved as well
	select.action.effect = "ALL",
	select.action.merge.option = "intersect",  # or "union"
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
	if (select.action.mode[1] != "ALL") {
		not.valid.action.mode <- setdiff(select.action.mode, kpred.mode)
		if (length(not.valid.action.mode) > 0) {
			warning("Given undefined action mode: ", paste0(not.valid.action.mode, collapse = ", ", ". "))
		}
		select.action.mode <- intersect(select.action.mode, kpred.mode)
		if (length(select.action.mode) == 0) {
			stop("No valid action mode is selected!")
		}
	}
	if (select.action.effect[1] != "ALL") {
		not.valid.action.effect <- setdiff(select.action.effect, kpred.action.effect)
		if (length(not.valid.action.effect) > 0) {
			warning("Given undefined aciton effect: ", paste0(not.valid.action.effect, collapse = ", "), ". ")
		}
		select.action.effect <- intersect(select.action.effect, kpred.action.effect)
		if (length(select.action.effect) == 0) {
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
	if (select.physical == TRUE) {  # currently it's the same to set select.action.mode => 'binding'
		if (select.action.mode[1] == "ALL") {
			select.action.mode <- "binding"
		} else {
			select.action.mode <- unique(c(select.action.mode, "binding"))
		}
	}
	if (select.action.mode[1] != "ALL") {
		tmp.sel.mode.pairs.df <- FastAlignPairs(this.actions.db[which(this.actions.db$mode %in% select.action.mode), 
			c("inter.GeneID.A", "inter.GeneID.B", "inter.GeneName.A", "inter.GeneName.B")], 4)
		tmp.sel.mode.pairs <- paste(tmp.sel.mode.pairs.df[, "inter.GeneID.A"], tmp.sel.mode.pairs.df[, "inter.GeneID.B"], sep = use.cut.symbol)
		use.action.pairs.list <- c(use.action.pairs.list, list(action.mode = tmp.sel.mode.pairs))
	}
	# in effect
	if (select.action.effect[1] != "ALL") {
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
		if ("positive" %in% select.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["activation"]], tmp.col.direction[["t"]])))
		}
		if ("negative" %in% select.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["inhibition"]], tmp.col.direction[["t"]])))	
		}
		if ("unspecified" %in% select.action.effect) {
			tmp.ret.inds <- c(tmp.ret.inds, Reduce(intersect, list(tmp.col.action[["non"]], tmp.col.direction[["t"]])))
		}
		if ("undirected" %in% select.action.effect) {
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
		use.action.pairs <- switch(select.action.merge.option,
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
		involved.genes <- unique(c(object@pairs.db$inter.GeneName.A, object@pairs.db$inter.GeneName.B))
		tmp.sel.pairs.conv <- paste(object@pairs.db[, "inter.GeneID.A"], object@pairs.db[, "inter.GeneID.B"], sep = use.cut.symbol.slim)
		tmp.sel.pairs.rev <- paste(object@pairs.db[, "inter.GeneID.B"], object@pairs.db[, "inter.GeneID.A"], sep = use.cut.symbol.slim)
		# slim genes properties
		object@anno.location.db <- object@anno.location.db[which(object@anno.location.db$Gene.name %in% involved.genes), ]
		object@anno.type.db <- object@anno.type.db[which(object@anno.type.db$Gene.name %in% involved.genes), ]
		object@go.ref.db <- object@go.ref.db[which(object@go.ref.db$Gene.name %in% involved.genes), ]
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






