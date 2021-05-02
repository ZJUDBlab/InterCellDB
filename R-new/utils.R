

# This indicate the standard data format for InterCellDB database on gene pairs
std.index.colname.end.dual <- 4

# This indicate the data format for processed InterCellDB gene pairs
proc.index.colname.end.dual <- 8

# For utility, illustrate the expression change options
kexprs.change <- c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn")

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
