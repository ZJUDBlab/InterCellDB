

# created at 2021.01.24, used to split commplex and cofactor to 
# make gene pairs

# copy from the source code in InterCellDB by remove the prefix "Tool."
SplitToGenDataFrame <- function(
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

# fetch database
cellchat.human.db <- readRDS("../../peer-db-packed/orig-CellChatDB-human-dblist-20210124.rds")
cellchat.mouse.db <- readRDS("../../peer-db-packed/orig-CellChatDB-mouse-dblist-20210124.rds")

# global variable definition
kDB.gsplit <- "->"  # split 2 genes in gene pairs like CCL2->CCR5


to.modf.db <- cellchat.mouse.db
ref.genes.db <- InterCellDB::genes.mouse.ref.db

# pre-processing cofactor and complex
cofactor.mat <- as.matrix(to.modf.db$cofactor)
cofactor.list <- apply(cofactor.mat, MARGIN = 1, function(x) {
	names(x) <- NULL
	x[which(x != "")]
	})
#
complex.mat <- as.matrix(to.modf.db$complex)
complex.list <- apply(complex.mat, MARGIN = 1, function(x) {
	names(x) <- NULL
	x[which(x != "")]
	})

# process interaction, split out cofactor and complex
# --- illustration of dealing cofactor and complex ---
# get 4 items in one interaction, ligand, receptor, cofactor, co-receptor
# The 4 parts get to interact with each other though, but no detailed info is given
# The solution is expand all parts, and regard them to interact between 2 groups, and 
# within their group. The 2 group are:
#   group 1: ligand + cofactor
#   group 2: receptor + co-receptor
prog.bar.m1 <- progress::progress_bar$new(total = nrow(to.modf.db$interaction))
prog.bar.m1$tick(0)
multiplex.interaction.result <- apply(to.modf.db$interaction, MARGIN = 1, 
	cofactor.list = cofactor.list, complex.list = complex.list, 
	function(x, cofactor.list, complex.list) {
		this.ligand <- x["ligand"]
		this.receptor <- x["receptor"]
		this.agonist <- x["agonist"]
		this.antagonist <- x["antagonist"]
		this.coA.receptor <- x["co_A_receptor"]
		this.coI.receptor <- x["co_I_receptor"]

		# get identity
		is.ligand.complex <- this.ligand %in% names(complex.list)
		if (is.ligand.complex == TRUE) {
			this.ligand <- complex.list[which(names(complex.list) == this.ligand)][[1]]
		}
		is.receptor.complex <- this.receptor %in% names(complex.list)
		if (is.receptor.complex == TRUE) {
			this.receptor <- complex.list[which(names(complex.list) == this.receptor)][[1]]
		}

		# get identity for cofactor (those are all recorded in cofactor.list)
		#if (this.agonist %in% names(cofactor.list)) {
		#	this.agonist <- cofactor.list[which(names(cofactor.list) == this.agonist)][[1]]
		#} else {
		#	this.agonist <- character(0)
		#}
		#if (this.antagonist %in% names(cofactor.list)) {
		#	this.antagonist <- cofactor.list[which(names(cofactor.list) == this.antagonist)][[1]]
		#} else {
		#	this.antagonist <- character(0)
		#}
		#if (this.coA.receptor %in% names(cofactor.list)) {
		#	this.coA.receptor <- cofactor.list[which(names(cofactor.list) == this.coA.receptor)][[1]]
		#} else {
		#	this.coA.receptor <- character(0)
		#}
		#if (this.coI.receptor %in% names(cofactor.list)) {
		#	this.coI.receptor <- cofactor.list[which(names(cofactor.list) == this.coI.receptor)][[1]]
		#} else {
		#	this.coI.receptor <- character(0)
		#}

		# group those genes
		group.lig.some <- c(this.ligand)    #, this.agonist, this.antagonist)
		group.rec.some <- c(this.receptor)  #, this.coA.receptor, this.coI.receptor)

		## fetch result 1
		res.inter.all <- as.character(outer(group.lig.some, group.rec.some, FUN = paste, sep = kDB.gsplit))
		res.df.inter <- SplitToGenDataFrame(res.inter.all, kDB.gsplit, res.colnames = c("gene.A", "gene.B"))
		## fetch result 2
		res.inner.lig.all <- as.character(outer(group.lig.some, group.lig.some, FUN = paste, sep = kDB.gsplit))
		res.inner.rec.all <- as.character(outer(group.rec.some, group.rec.some, FUN = paste, sep = kDB.gsplit))
		# remove itself cross gene pairs
		res.df.inner.lig <- SplitToGenDataFrame(res.inner.lig.all, kDB.gsplit, res.colnames = c("gene.A", "gene.B"))
		res.df.inner.rec <- SplitToGenDataFrame(res.inner.rec.all, kDB.gsplit, res.colnames = c("gene.A", "gene.B"))
		res.df.inner <- rbind(res.df.inner.lig, res.df.inner.rec)
		inds.self.cross <- which(res.df.inner$gene.A == res.df.inner$gene.B)
		res.df.inner <- res.df.inner[setdiff(seq_len(nrow(res.df.inner)), inds.self.cross), ]

		prog.bar.m1$tick()
		## result
		cbind(cbind(rbind(res.df.inter, res.df.inner), evidence = x["evidence"]), annotation = x["annotation"])
	})
multiplex.interaction.result <- dplyr::bind_rows(multiplex.interaction.result)

## further process
# use remapping to get gene pairs right
extra.cols <- c("evidence", "annotation")
inter.A.it <- multiplex.interaction.result[, c("gene.A", extra.cols), drop = FALSE]
inter.A.it <- cbind(inter.A.it, row.id = seq_len(nrow(inter.A.it)))
inter.B.it <- multiplex.interaction.result[, "gene.B", drop = FALSE]  # extra infos only by A, it's once enough
inter.B.it <- cbind(inter.B.it, row.id = seq_len(nrow(inter.B.it)))

new.inter.A.it <- InterCellDB::Tool.AddUserRestrictDB(inter.A.it, ref.genes.db)
new.inter.B.it <- InterCellDB::Tool.AddUserRestrictDB(inter.B.it, ref.genes.db)
new.inter.A.it <- new.inter.A.it[order(new.inter.A.it$row.id), ]
new.inter.B.it <- new.inter.B.it[order(new.inter.B.it$row.id), ]
# in running Human database, no warning get
# in running Mouse database, it get 1 gene cannot be remapped, is H2-BI

new.inter.db <- dplyr::left_join(new.inter.A.it, new.inter.B.it, by = ("row.id" = "row.id"))
new.inter.db$row.id <- NULL
# rearrange the cols
new.inter.db <- new.inter.db[, c(1,5,2,6,3,4)]
colnames(new.inter.db)[1:4] <- c("GeneID.A", "GeneID.B", "GeneName.A", "GeneName.B")

# use FastAlignPairs, DoPartUnique, ReverseOddEvenCols to get unique gene pairs.
# Get FastAlignPairs, DoPartUnique, ReverseOddEvenCols function from 
# /Users/jinziyang/Sc-RNAsequence/interaction-database/mus-10090-workflow/R-work/LRE-extend-annotation-v1.R
new.inter.db <- FastAlignPairs(new.inter.db, 4)
new.inter.db <- DoPartUnique(new.inter.db, c(1,2))
#
# here, to save RDS
# turnback, and run on mouse genes
#
# Human database, get 2147 gene pairs
# Mouse database, get 2230 gene pairs




