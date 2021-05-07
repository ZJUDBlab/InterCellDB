


#' Add user-defined dataset
#' 
#' @description
#' This function gives a standlized procedure that package users can easily add its 
#' definitions about genes and their functions or types or locations etc.
#'
#' @param user.def.db Data.frame. It should be a data.frame with at least 2 columns, whose 1st column is gene names, and the 2nd is 
#' definitions in users' personal will, and colums from the 3rd will be inherited intactly.
#' @inheritParams Inside.DummyGenesRefDB
#' @param warning.given Character. It applies extra limitation on the types(molecular functions) of A in gene pairs formatted as A-B.
#'
#' @return Data.frame.
#'
#' | <col-1> <col-2> <col-3> ...|
#'
#' | GeneID Gene.name user.type ...|
#'
#' | <data> <data> <data> ...|
#'
#'
#'
#'
#' @importFrom dplyr left_join
#'
#' @export
#'
Tool.AddUserRestrictDB <- function(
	user.def.db,
	genes.ref.db,
	warning.given = "genes",
	if.used.inside = TRUE
) {
	if ("cluster" %in% colnames(user.def.db)) {
		stop("Program internally use 'cluster' for some operations. Please rename columns of the raw data.")
	}
	colnames(user.def.db)[1] <- "gene"
	user.def.db$cluster <- user.def.db[, 2]
	user.def.db <- DataPrep.RemapClustersMarkers(user.def.db, genes.ref.db, warning.given, if.used.inside)$result
	colnames(user.def.db)[1] <- "Gene.name"
	# add GeneID
	entrez.db <- genes.ref.db$gene.ncbi.db
	user.def.db <- left_join(user.def.db, entrez.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
		by = c("Gene.name" = "Symbol_from_nomenclature_authority"))
	tmp.ncol <- length(colnames(user.def.db))
	user.def.db <- user.def.db[, c(tmp.ncol, 1:(tmp.ncol-1))]
	# return
	user.def.db  # GeneID Gene.name user.type *1 *2 *3 ...
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
#' @export
#'
Tool.formula.onLogFC.default <- function(
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
#' @export
#'
Tool.formula.onPValAdj.default <- function(
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
		max.f <- max(tmp.f[inds.tmp.f])  # if length(inds.tmp.f) == 0, max returns -Inf
		max.f	<- ifelse(is.infinite(max.f), default.max.replace, max.f)
		max.b <- max(tmp.b[inds.tmp.b])  # if length(inds.tmp.b) == 0, max returns -Inf
		max.b	<- ifelse(is.infinite(max.b), default.max.replace, max.b)
		tmp.f[which(is.infinite(tmp.f))] <- 10 * max.f
		tmp.b[which(is.infinite(tmp.b))] <- 10 * max.b
	}
	return(tmp.f * tmp.b)
}





#' generate gene pairs in standard format from VEinfos
#' 
#' @description
#' This function generates gene pairs in standard format(in data frame), 
#' and gets these pairs easier to be compared with others.
#'
#' @inheritParams Inside.DummyVEinfos
#'
#' @details
#' The standard format in this package is that gene pairs are maintained in data.frame, and the 2 genes 
#' participated in each gene pair are recorded in columns named "inter.GeneName.A" and "inter.GeneName.B".
#'
#' @importFrom dplyr left_join
#'
#'
#'
#' @export
#'
Tool.GenStdGenePairs.from.VEinfos <- function(
  VEinfos
) {
  vertices.infos <- VEinfos$vertices.infos
  edges.infos <- VEinfos$edges.infos
  #
  tmp.res <- left_join(edges.infos[, c("from", "to")], vertices.infos[, c("UID", "ClusterName", "GeneName", "LogFC", "PValAdj")], by = c("from" = "UID"))
  colnames(tmp.res)[c(ncol(tmp.res) - 3:0)] <- c("inter.Cluster.A", "inter.GeneName.A", "inter.LogFC.A", "inter.PValAdj.A")
  tmp.res <- left_join(tmp.res, vertices.infos[, c("UID", "ClusterName", "GeneName", "LogFC", "PValAdj")], by = c("to" = "UID"))
  colnames(tmp.res)[c(ncol(tmp.res) - 3:0)] <- c("inter.Cluster.B", "inter.GeneName.B", "inter.LogFC.B", "inter.PValAdj.B")
  # form std data.frame
  align.colnames <- c("inter.GeneName.A", "inter.GeneName.B", "inter.LogFC.A", "inter.LogFC.B", "inter.PValAdj.A", "inter.PValAdj.B", "inter.Cluster.A", "inter.Cluster.B")
  tmp.res <- tmp.res[, match(align.colnames, colnames(tmp.res))]
  # result
  std.df <- DoPartUnique(tmp.res, 1:2)
  # match cluster
  # get conv ones
  std.res.conv <- std.df[intersect(which(std.df$inter.Cluster.A == VEinfos$cluster.name.A), which(std.df$inter.Cluster.B == VEinfos$cluster.name.B)), ]
  # get rev ones
  std.res.rev <- std.df[intersect(which(std.df$inter.Cluster.A == VEinfos$cluster.name.B), which(std.df$inter.Cluster.B == VEinfos$cluster.name.A)), ]
  std.res.rev <- std.res.rev[, ReverseOddEvenCols(length(align.colnames))]  # reverse all paired columns
  colnames(std.res.rev) <- colnames(std.res.conv)
  # get the result
  std.res.all <- rbind(std.res.conv, std.res.rev)

  return(std.res.all)
}





#' Find genes annotated in specific GO terms
#' 
#' @description
#' This function uses GO terms or GO IDs to get a specific list of genes.
#'
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
#'
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






