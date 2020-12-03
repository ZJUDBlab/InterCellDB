


#' generate gene pairs from VEinfos in standard format
#' 
#' @description
#' This function generates gene pairs in standard format(in data frame), 
#' and gets these pairs easier to be compared with others.
#'
#' @param VEinfos List. It contains informations about vertices and edges, and is exactly return value of
#' \code{GenerateVEinfos()} or \code{TrimVEinfos()}.
#'
#' @details
#' The standard format in this package is that gene pairs are maintained in data.frame, and the 2 genes 
#' participated in each gene pair are recorded in columns named "inter.GeneName.A" and "inter.GeneName.B".
#'
#'
#'
#' @export
#'
Tool.GenStdGenePairs <- function(
  VEinfos
) {
  vertices.infos <- VEinfos$vertices.infos
  edges.infos <- VEinfos$edges.infos
  #
  inds.e.from.match <- match(edges.infos$from, vertices.infos$UID)
  inds.e.to.match <- match(edges.infos$to, vertices.infos$UID)
  # from to data.frame
  std.df <- data.frame("inter.GeneName.A" = vertices.infos$GeneName[inds.e.from.match], 
    "inter.GeneName.B" = vertices.infos$GeneName[inds.e.to.match], stringsAsFactors = FALSE)
  std.df <- DoPartUnique(std.df, 1:2)
  return(std.df)
}





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
	warning.given = "genes"
) {
	colnames(user.def.db)[1] <- "gene"
	user.def.db <- DataPrep.RemapClustersMarkers(user.def.db, genes.ref.db, warning.given)
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
#' 
#' @return Character. A gene list of given GO IDs or terms.
#'
#'
#'
#' @export
#'
Tool.FindGenesFromGO <- function(
	go.todolist,
	genes.ref.db,
	go.ref.db
) {
	# pre-process
	entrez.ref.db <- genes.ref.db$gene.ncbi.db
	### --- doing all checks ---
	# 1. check if GO_IDs or GO_terms given are available
	# 2. auto transform of all go.todolist to be GO_ID
	# 3. get spliting
	inds.ID.given <- grep("^GO:", go.todolist)
	if (length(inds.ID.given) > 0) {
		go.ID.given.list <- go.todolist[inds.ID.given]
		go.term.given.list <- go.todolist[-inds.ID.given]
	} else {
		go.ID.given.list <- NULL
		go.term.given.list <- go.todolist
	}
	## Giving error report
	# ID matching
	go.ID.given.nonexist <- character()
	go.ID.given.exist <- character()
	if (length(go.ID.given.list) > 0) {
		inds.mID.s <- match(go.ID.given.list, go.ref.db$GO_ID)
		go.ID.given.nonexist <- go.ID.given.list[which(is.na(inds.mID.s))]
		go.ID.given.exist <- go.ID.given.list[which(!is.na(inds.mID.s))]
		if (length(go.ID.given.nonexist) > 0) {
			warning("The following GO_IDs are not found: \n ", paste0(go.ID.given.nonexist, collapse = ", "), ".")
		}
	}
  res.go.id.match.list <- go.ref.db[which(go.ref.db$GO_ID %in% go.ID.given.exist), "GeneName"]
	# term matching
	go.term.given.nonexist <- character()
	go.term.given.exist <- character()
	if (length(go.term.given.list) > 0) { 
		inds.mTerm.s <- match(go.term.given.list, go.ref.db$GO_term)
		go.term.given.nonexist <- go.term.given.list[which(is.na(inds.mTerm.s))]
		go.term.given.exist <- go.term.given.list[which(!is.na(inds.mTerm.s))]
		if (length(go.term.given.nonexist) > 0) {
			warning("The following GO_terms are not found: \n ", paste0(go.term.given.nonexist, collapse = ", "), ".")
		}
	}
  res.go.term.match.list <- go.ref.db[which(go.ref.db$GO_term %in% go.term.given.exist), "GeneName"]
	# finish transformation
	res.go.rel.genes <- unique(as.character(c(res.go.id.match.list, res.go.term.match.list)))
	res.go.rel.genes <- res.go.rel.genes[order(res.go.rel.genes)]
  # return
	res.go.rel.genes
}






