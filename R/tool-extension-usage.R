
# 1 GeneID Gene.name user.type
# 2 mapping Gene.name to be authorized and recorded in our database

# user.def.db, 2col-data.frame. 1-gene name, 2-user type. 2nd+ colname will be inherited.
# user database colname should add user.* to avoid anything unexpected runtime error.

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
	user.def.db  # GeneID Gene.name user.type1 *2 *3 ...
}



# @param go.todolist Character. Several GO_terms or GO_IDs or mixed, which will 
# be used to get subsets of feature genes.
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
		for (i in 1:length(go.ID.given.list)) {
			ind.tmp <- match(go.ID.given.list[i], go.ref.db$GO_ID)
			if (is.na(ind.tmp)) { go.ID.given.nonexist <- append(go.ID.given.nonexist, go.ID.given.list[i]) }
			else { go.ID.given.exist <- append(go.ID.given.exist, go.ID.given.list[i]) }
		}
		if (length(go.ID.given.nonexist) > 0) {
			warning("The following GO_IDs are not found: \n ", paste0(go.ID.given.nonexist, collapse = ", "), ".")
		}
	}
	# term matching
	go.term.given.nonexist <- character()
	go.term.given.exist <- character()
	if (length(go.term.given.list) > 0) { 
		for (i in 1:length(go.term.given.list)) {
			ind.tmp <- match(go.term.given.list[i], go.ref.db$GO_term)
			if (is.na(ind.tmp)) { go.term.given.nonexist <- append(go.term.given.nonexist, go.term.given.list[i]) }
			else {  # here doing the term-ID tranformation
				go.term.given.exist <- append(go.term.given.exist, go.ref.db[ind.tmp, "GO_ID"])
			}
		}
		if (length(go.term.given.nonexist) > 0) {
			warning("The following GO_terms are not found: \n ", paste0(go.term.given.nonexist, collapse = ", "), ".")
		}
	}
	# finish transformation
	go.final.useID <- c(go.ID.given.exist, go.term.given.exist)
	res.go.rel.genes <- unique(go.ref.db[which(go.ref.db$GO_ID %in% go.final.useID), "GeneName"])  # [TODO] GeneName -> Gene.name
	# return
	res.go.rel.genes
}




