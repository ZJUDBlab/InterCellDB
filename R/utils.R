

# This indicate the standard data format for InterCellDB database on gene pairs
std.index.colname.end.dual <- 4

# This indicate the data format for processed InterCellDB gene pairs
proc.index.colname.end.dual <- 8

#' Options for Expression Change
#'
#' \code{kexprs.change} defines the options to be used for selecting 
#' gene pairs under specific expression-change status.
#'
#' @rdname pred-exprs-change
#' @export
#'
kexprs.change <- c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn")



#' Check Species
#'
#' This function helps to check species.
#'
#' @param species Allowed species are 'human' and 'mouse'.
#'
#'
CheckSpeciesValidity <- function(
	species
) {
	# Check species, only support human and mouse yet
	DB.support.species <- c("human", "mouse")
	if (length(species) != 1 || is.na(match(species[1], DB.support.species))) {
		stop("InterCellDB only supports species to be 'human' and 'mouse' for now.")
	}
	return(species)
}

# This function is to provide standard format to check some parameters.
CheckParamStd <- function(
	input.opt,
	allowed.opt,
	opt.name = "opt",
	stop.on.zero = FALSE
) {
	not.valid.opt <- setdiff(input.opt, allowed.opt)
	if (length(not.valid.opt) > 0) {
		warning("Given undefined ", opt.name, ": ", paste0(not.valid.opt, collapse = ", ", ". "))
	}
	input.opt <- intersect(input.opt, allowed.opt)
	if (stop.on.zero == TRUE && length(input.opt) == 0) {
		stop("No valid ", opt.name, " is selected!")
	}
	return(input.opt)
}

# Transform 4 action effect to corresponding 7 extended action effect
TransActionEffectToActId <- function(
	action.effect
) {
	ret.act.ids <- as.character(unlist(lapply(action.effect, function(x) {
		ret <- switch(x, 
			"positive" = c("A-->B", "A<--B"),
			"negative" = c("A--|B", "A|--B"),
			"unspecified" = c("A--oB", "Ao--B"),
			"undirected" = "A---B",
			stop("Undefined action effect!")
		)
	})))
	return(ret.act.ids)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# more common used tool function
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Split Character To Data.frame
#' 
#' @description
#' This function generates n-column data frame by spliting character by some specific letters or phrases.
#'
#' @param to.splits.string Character. The string to be splited. 
#' @param to.split.by Character(1). The string will be splited by this parameter, and it will be directly 
#' passed to parameter \code{split} of function \code{strsplit}. 
#' @param res.colnames Character(2). It gives 2 column names in the result table.
#' The splited string gets to have 2 parts of equal length, and 
#' they will be reconstructed to form a table so that 2 column names are needed. 
#'
#'
#'
#' @export
#'
Tool.SplitToGenDataFrame <- function(
  to.splits.string, 
  to.split.by, 
  res.colnames
) {
	tmp.splits <- strsplit(to.splits.string, split = to.split.by, fixed = TRUE)
	tmp.len.splits <- as.integer(unlist(lapply(tmp.splits, FUN = length)))
	tmp.len.splits <- unique(tmp.len.splits)
	if (length(tmp.len.splits) != 1) {
		stop("Given data cannot be uniformly splited, with different splits: ", 
			paste0(tmp.len.splits, collapse = ", "), ".")
	}
	tmp.merge.elems <- as.character(unlist(tmp.splits))
	tmp.merge.base.index <- seq_len(length(tmp.merge.elems) / tmp.len.splits)
	res.prep.list <- list()
	for (i in seq_len(tmp.len.splits)) {
		tmp.indices <- tmp.merge.base.index * tmp.len.splits - (tmp.len.splits - i)
		res.prep.list <- c(res.prep.list, list(tmp.merge.elems[tmp.indices]))
	}
	res.df <- data.frame(res.prep.list, stringsAsFactors = FALSE)
	if (length(res.colnames) == ncol(res.df)) {
		colnames(res.df) <- res.colnames
	} else {
		warning("Given colnames are not matched with the result columns, and unexpected errors may happen!")
		if (length(res.colnames) > ncol(res.df)) {
			colnames(res.df) <- res.colnames[seq_len(ncol(res.df))]
		} else {
			if (length(res.colnames) > 0) {
				colnames(res.df)[seq_along(res.colnames)] <- res.colnames
			}
		}
	}
	return(res.df)
}



# [inside usage]
# get from ?toupper, .simpleCap
Tc.Cap.simple <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
}

# [inside usage]
# parallel version of Tc.Cap.simple
Tc.Cap.simple.vec <- function(to.cap.vec) {
	unlist(lapply(to.cap.vec, FUN = Tc.Cap.simple))
}



#' Speed-up Doing Unique for Data.frame
#'
#' @description
#' This function uses the properties of \code{data.frame} and \code{rownames}, and 
#' makes partial unique process fast and effective.
#'
#' @param xxpairs Data.frame. Any data.frame object.
#' @param cols.select Integer. Specify some columns in integer vector, e.g. c(1:2).
#' Its length \bold{must be >= 2}, or the function will not work properly. 
#'
#' @details
#' When encountering multi-columns tables, e.g. 6 columns or more, \code{unique()} will 
#' be really slow if it is applied on the whole table. However, in most circumstances, 
#' there is no need to apply \code{unique()} on all columns, i.e. only do unique on some 
#' columns, which is exactly the thing this function does.
#'
#' @export
#'
DoPartUnique <- function(
	xxpairs,
	cols.select = c(1:2)
) {
	if (!all(cols.select %in% seq_len(ncol(xxpairs)))) {
		stop("Some columns selected are undefined! Please check again!")
	}
	# rownames(xxpairs) <- NULL
	tmp.uni <- xxpairs[, cols.select]
	rownames(tmp.uni) <- NULL
	tmp.uni <- unique(tmp.uni)
	xxpairs <- xxpairs[as.integer(rownames(tmp.uni)),]
	xxpairs
}



#' Permutate Reverse Order of Odds and Evens
#'
#' @description
#' This function generates reverse premutation of odd values and even values 
#' in continues values, e.g. \code{1:10}.
#'
#' @param ncols Integer. An arbitrary integer.
#'
#' @examples
#' # for 1:10 
#' ReverseOddEvenCols(10)
#' # the reuslt is c(2,1,4,3,6,5,8,7,10,9)
#'
#' @export
#'
ReverseOddEvenCols <- function(
	ncols
) {
	len.all <- ncols %/% 2
	val.even <- 2 * 1:len.all
	val.odds <- (2 * 1:len.all) -1
	new.serial <- NULL
	for (i in 1:len.all) {
		new.serial <- c(new.serial, val.even[i], val.odds[i])
	}
	if (ncols %% 2 != 0) {  # odds cols
		# leave last one unchanged
		new.serial <- c(new.serial, ncols)
	}
	#end# return
	new.serial
}



#' Fast Align Gene Pairs
#'
#' This function helps to align gene pairs by their GeneID or GeneName. Its purpose is
#' to make different gene pairs data comparable in the same settings.
#'
#' @param xxpairs Gene pairs to be compared.
#' @param ind.colname.end.dual Gene Pairs are given as pairs, so as their properties.
#' @param use.cols Use which columns to align those pairs. Column ID or column name is supported.
#' @param use.class Defines what is given in \code{use.cols}.
#'
#' @export
#'
FastAlignPairs <- function(
	xxpairs, 
	ind.colname.end.dual, 
	use.cols = c(1,2),
	use.class = "integer"
) {
	inds.rev <- switch(use.class, 
		"integer" = which(as.integer(xxpairs[, use.cols[1]]) >= as.integer(xxpairs[, use.cols[2]])),
		"character" = which(as.character(xxpairs[, use.cols[1]]) >= as.character(xxpairs[, use.cols[2]])),
		"numeric" = which(as.numeric(xxpairs[, use.cols[1]]) >= as.numeric(xxpairs[, use.cols[2]])),
		stop("Unspported type of input!")
	)
	xxpairs.result <- NULL
	if (length(inds.rev) > 0) {
		xxpairs.rev <- xxpairs[inds.rev, ]
		xxpairs.conv <- xxpairs[-inds.rev, ]
		if (ncol(xxpairs.rev) <= ind.colname.end.dual) {
			if (ncol(xxpairs.rev) < ind.colname.end.dual) {
				stop("Given database has less cols than given parameter: ", ind.colname.end.dual)
			}
			xxpairs.rev.rem <- xxpairs.rev[, c(ReverseOddEvenCols(ind.colname.end.dual))]
		} else {
			xxpairs.rev.rem <- xxpairs.rev[, c(ReverseOddEvenCols(ind.colname.end.dual), (ind.colname.end.dual+1):ncol(xxpairs.rev))]
		}
		colnames(xxpairs.rev.rem) <- colnames(xxpairs.conv)
		xxpairs.result <- rbind(xxpairs.conv, xxpairs.rev.rem)
	} else {
		# do nothing and return	
		xxpairs.result <- xxpairs
	}
	# return
	xxpairs.result
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tool function for specific purpose
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find Genes Annotated in Specific GO Terms
#' 
#' @description
#' This function uses GO terms or GO IDs to get a specific list of genes.
#'
#' @param go.todolist Character. Several GO_terms or GO_IDs or mixed, which will 
#' be used to get subsets of feature genes. 
#' @inheritParams Inside.DummyGenesRefDB
#' @inheritParams Inside.DummyGORefDB
#' @param go.use.relative Character. If set TRUE, GO IDs or Terms will be explored by its relatives of GO term tree. 
#' This is useful as the lite version of GO database will only records the gene and its most accurate GO terms.
#' @param go.relative.option Character. 4 options: ancestor, parents, offspring, children. 
#' With each corresponding to one group of database listed in \pkg{GO.db}.
#' 
#' @return Character. A gene list of given GO IDs or terms.
#'
#' @importFrom GO.db GOCCANCESTOR GOCCPARENTS GOCCOFFSPRING GOCCCHILDREN 
#' @importFrom GO.db GOMFANCESTOR GOMFPARENTS GOMFOFFSPRING GOMFCHILDREN 
#' @importFrom GO.db GOBPANCESTOR GOBPPARENTS GOBPOFFSPRING GOBPCHILDREN
#'
#' @export
#'
Tool.FindGenesFromGO <- function(
	go.todolist,
	genes.ref.db,
	go.ref.db,
	go.use.relative = TRUE, 
	go.relative.option = "offspring"
) {
	# pre-process
	entrez.ref.db <- genes.ref.db$gene.ncbi.db
	### --- doing all checks ---
	# 1. check if GO_IDs or GO_terms given are available
	# 2. auto transform of all go.todolist to be GO_ID
	# 3. get spliting
	inds.ID.given <- grep("^GO:", go.todolist)
	go.ID.given.list <- go.todolist[inds.ID.given]
	go.term.given.list <- go.todolist[setdiff(seq_along(go.todolist), inds.ID.given)]
	## Giving error report
	# term matching
	go.term.transID <- go.term.given.exist <- go.term.given.nonexist <- character()
	if (length(go.term.given.list) > 0) {
		inds.mTerm.s <- match(go.term.given.list, go.ref.db$GO_term)
		go.term.given.nonexist <- go.term.given.list[which(is.na(inds.mTerm.s))]
		go.term.given.exist <- go.term.given.list[which(!is.na(inds.mTerm.s))]
		if (length(go.term.given.nonexist) > 0) {
			warning("The following GO_terms are not found: \n ", paste0(go.term.given.nonexist, collapse = ", "), ".")
		}
		# trans the terms to be IDs
		go.term.transID <- go.ref.db[inds.mTerm.s, "GO_ID"]
	}

	# merge the result to go.ID.given.list
	go.ID.given.list <- c(go.ID.given.list, go.term.transID)
	# then ID matching
	go.ID.given.exist <- go.ID.given.nonexist <- character()
	if (length(go.ID.given.list) > 0) {
		inds.mID.s <- match(go.ID.given.list, go.ref.db$GO_ID)
		go.ID.given.nonexist <- go.ID.given.list[which(is.na(inds.mID.s))]
		go.ID.given.exist <- go.ID.given.list[which(!is.na(inds.mID.s))]
		if (length(go.ID.given.nonexist) > 0) {
			warning("The following GO_IDs are not found: \n ", paste0(go.ID.given.nonexist, collapse = ", "), ".")
		}
	} 
	
	# check if use relative and the option
	if (go.use.relative == TRUE) {
		go.relative.db <- switch(go.relative.option,
			"ancestor"  = list("Component" = as.list(GO.db::GOCCANCESTOR), "Function" = as.list(GO.db::GOMFANCESTOR), "Process" = as.list(GO.db::GOBPANCESTOR)), 
			"parents"   = list("Component" = as.list(GO.db::GOCCPARENTS), "Function" = as.list(GO.db::GOMFPARENTS), "Process" = as.list(GO.db::GOBPPARENTS)), 
			"offspring" = list("Component" = as.list(GO.db::GOCCOFFSPRING), "Function" = as.list(GO.db::GOMFOFFSPRING), "Process" = as.list(GO.db::GOBPOFFSPRING)), 
			"children"  = list("Component" = as.list(GO.db::GOCCCHILDREN), "Function" = as.list(GO.db::GOMFCHILDREN), "Process" = as.list(GO.db::GOBPCHILDREN)), 
			stop("Undefined options! Please recheck the given parameter!")
		)

		# extract the genes
		go.rel.list <- lapply(go.ID.given.exist, go.ref.db = go.ref.db, go.relative.db = go.relative.db, 
			function(x, go.ref.db, go.relative.db) {
				# take the category
				tmp.template <- match(x, go.ref.db$GO_ID)
				this.category <- go.ref.db[tmp.template, "Category"]
				tmp.rel.db <- go.relative.db[[which(names(go.relative.db) == this.category)]]
				# get relative part
				tmp.rel.inds <- which(names(tmp.rel.db) == x)
				tmp.res <- character(0)
				if (length(tmp.rel.inds) != 0) {
					tmp.rel.go.s <- tmp.rel.db[[tmp.rel.inds]]
					tmp.res <- unique(go.ref.db[which(go.ref.db$GO_ID %in% tmp.rel.go.s), "Gene.name"])
				}
				# get plain part
				tmp.plain.res <- unique(go.ref.db[which(go.ref.db$GO_ID %in% x), "Gene.name"])
				return(unique(c(tmp.res, tmp.plain.res)))
			})	
	} else {
		# use the plain GO IDs or Terms
		go.rel.list <- lapply(go.ID.given.exist, go.ref.db = go.ref.db,
			function(x, go.ref.db) {
				unique(go.ref.db[which(go.ref.db$GO_ID %in% x), "Gene.name"])
			})
	}

	## re-align with the given order in parameter go.todolist
	# reload raw parameter segment
	inds.ID.given <- grep("^GO:", go.todolist)
	go.ID.given.list <- go.todolist[inds.ID.given]
	go.term.given.list <- go.todolist[setdiff(seq_along(go.todolist), inds.ID.given)]
	# those IDs
	tmp.valid.IDs <- match(go.ID.given.list, go.ID.given.exist)
	tmp.valid.IDs <- tmp.valid.IDs[which(!is.na(tmp.valid.IDs))]
	# thoes Terms
	go.trans.ID.Term.exist <- go.ref.db[match(go.ID.given.exist, go.ref.db$GO_ID), "GO_term"]
	tmp.valid.Terms <- match(go.term.given.list, go.trans.ID.Term.exist)
	tmp.valid.Terms <- tmp.valid.Terms[which(!is.na(tmp.valid.Terms))]
	# collect result and reorder by their given order
	tmp.names <- c(go.ID.given.exist[tmp.valid.IDs], go.trans.ID.Term.exist[tmp.valid.Terms])
	go.res.final <- go.rel.list[c(tmp.valid.IDs, tmp.valid.Terms)]
	names(go.res.final) <- tmp.names
	tmp.reorder.inds <- match(go.todolist, tmp.names)
	go.res.final <- go.res.final[tmp.reorder.inds[which(!is.na(tmp.reorder.inds))]]

	return(go.res.final)
}



#' Generate Permutations of Count Matrix
#'
#' This function uses cell label permutation to generate one list of normalized count matrix.
#'
#' @param object A \code{InterCell} object.
#' @param count.matrix Normalized count matrix, with columns as cell barcodes and rows as gene names.
#' @param cells.meta cell group labels, which are used to group cells in several groups rather all different barcodes.
#' @param genes.name Readable gene names, which should be the same set as genes given in \code{DEG.table}.
#'  To be noted, if genes are remapped when created \code{InterCell} object, it will be automatically applied to this.
#' @param perm.times Permutation times. Default is set to 1000.
#' @param p.seed Seed for generating permutations.
#'
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply 
#' @import Matrix
#'
#' @export
#'
Tool.GenPermutation <- function(
	object,
	count.matrix,
	cells.meta = NULL,
	genes.name = NULL,
	perm.times = 1000,
	p.seed = 101
) {
	if (!is.null(p.seed)) {
		set.seed(p.seed)
	}

	if (is.null(ncol(count.matrix)) || ncol(count.matrix) == 0 ||
		is.null(nrow(count.matrix)) || nrow(count.matrix) == 0) {
		stop("Please provided valid matrix or data.frame to use.")
	}
	
	if (is.data.frame(count.matrix)) {  # transform data.frame to matrix
		count.matrix <- as.matrix(count.matrix)
	} else {
		count.matrix <- as(count.matrix, "sparseMatrix")
	}
	
	if (!is.null(cells.meta)) {  # use cells.meta to overwrite the columns of count.matrix
		if (length(cells.meta) != ncol(count.matrix)) {
			stop("Provided meta data on cells is not equal to actual cell count (columns in count matrix).")
		}
		cells.meta <- as.character(cells.meta)  # force character
	}
	
	if (!is.null(genes.name)) {
		if (length(genes.name) != nrow(count.matrix)) {
			stop("Provided gene names are not equal to rows in count matrix.")
		}
		genes.name <- as.character(genes.name)  # force character
	}
	# set dimnames
	if (!is.null(cells.meta)) {
		colnames(count.matrix) <- cells.meta
	}
	if (!is.null(genes.name)) {
		rownames(count.matrix) <- genes.name
	}


	# check if matrix are unqiue in colnames and rownames
	if (anyDuplicated(rownames(count.matrix))) {
		warning("Provided gene names or original rownames of count matrix have duplicated items, which may cause unexpected result!")
	}

	# remapping genes
	if (object@misc$if.remap.genes == TRUE) {
		dummy.gdf <- data.frame(gene = rownames(count.matrix), ID = seq_len(nrow(count.matrix)), cluster = "IT", stringsAsFactors = FALSE)
		dummy.remap.res <- DataPrep.RemapClustersMarkers(dummy.gdf, object@database@species, if.used.inside = TRUE, final.dup.rm = FALSE)
		colnames(dummy.gdf)[1] <- c("orig.gene")
		dummy.gdf[, "mapped.gene"] <- dummy.remap.res$result[order(dummy.remap.res$result$ID), "gene"]
		rownames(count.matrix) <- dummy.gdf[match(rownames(count.matrix), dummy.gdf[, "orig.gene"]), "mapped.gene"]
	}


	cellgroup <- unique(colnames(count.matrix))
	run.func <- function(x, counts, cellgroup) {
		ret.data <- matrix(data = NA, nrow = nrow(count.matrix), ncol = length(cellgroup), 
						   dimnames = list(rownames(count.matrix), cellgroup))
		# cell label permutation
		new.cells.meta <- sample(colnames(count.matrix), ncol(count.matrix), replace = FALSE)
		dimnames(count.matrix) <- list(rownames(count.matrix), new.cells.meta)
		
		for (i in cellgroup) {
			tmp.ind.it <- which(colnames(count.matrix) == i)
			ret.data[, i] <- as.numeric(rowSums(count.matrix[, tmp.ind.it, drop = FALSE])) / length(tmp.ind.it)
		}
		return(as(ret.data, "sparseMatrix"))
	}
	# run permutation expression
	if (nbrOfWorkers() == 1) {
		tmp.perm.retlist <- pblapply(seq_len(perm.times), counts = count.matrix, cellgroup = cellgroup,
			FUN = run.func
		)
	} else {
		tmp.perm.retlist <- future_lapply(seq_len(perm.times), counts = count.matrix, cellgroup = cellgroup,
			FUN = run.func, future.seed = p.seed
		)
	}

	# add actual expression matrix and average gene expression over all cell clusters
	tmp.actual.list <- list()
	# actual expression matrix
	actual.exprs.mat <- matrix(data = NA, nrow = nrow(count.matrix), ncol = length(cellgroup),
							   dimnames = list(rownames(count.matrix), cellgroup))
	for (i in cellgroup) {
		tmp.ind.it <- which(colnames(count.matrix) == i)
		actual.exprs.mat[, i] <- as.numeric(rowSums(count.matrix[, tmp.ind.it, drop = FALSE])) / length(tmp.ind.it)
	}
	actual.exprs.mat <- as(actual.exprs.mat, "sparseMatrix")
	# all average
	actual.all.avg.data <- rowSums(count.matrix) / ncol(count.matrix)
	names(actual.all.avg.data) <- rownames(count.matrix)
	tmp.actual.list <- c(tmp.actual.list, actual = list(actual.exprs.mat), allavg = list(actual.all.avg.data))

	# return
	tmp.ret.list <- list(actual = tmp.actual.list, perm = tmp.perm.retlist)
	
	return(tmp.ret.list)
}



