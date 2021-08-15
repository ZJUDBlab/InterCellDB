
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database

# fetch database
cellchat.human.db <- readRDS("../peer-db-packed/orig-CellChatDB-human-dblist-20210124.rds")
cellchat.mouse.db <- readRDS("../peer-db-packed/orig-CellChatDB-mouse-dblist-20210124.rds")
### panel for extra informations
# 0. human and mouse
# 1. get interaction annotated as 3 categories: "Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling".
# 2. it gives the source of interaction pairs, which is KEGG, PMID, PMC. 
#    The latter 2 are seen experimentally confirmed as the article said.
#    So, like InterCellDB, it gets to split to Experiments and Knowledge database.
#
# warn1: The structure readin is list of data.frame, here to noted.
#


# fetch processecd database (list constructure removed, only gene pairs)
cellchat.proc.human.db <- readRDS("../peer-db-packed/processed-CellChatDB-human-20210125.rds")
cellchat.proc.mouse.db <- readRDS("../peer-db-packed/processed-CellChatDB-mouse-20210125.rds")

# fetch complex sublist
# human
complex.human.mat <- as.matrix(cellchat.human.db$complex)
complex.human.list <- apply(complex.human.mat, MARGIN = 1, function(x) {
	names(x) <- NULL
	x[which(x != "")]
	})
rm(complex.human.mat)
# mouse
complex.mouse.mat <- as.matrix(cellchat.mouse.db$complex)
complex.mouse.list <- apply(complex.mouse.mat, MARGIN = 1, function(x) {
	names(x) <- NULL
	x[which(x != "")]
	})
rm(complex.mouse.mat)

## shared part
tmp.used.db <- cellchat.human.db
# #1# extract Gene.group.name information from its sub-DB geneInfo
tmp.split.gene.group <- as.character(unlist(strsplit(tmp.used.db$geneInfo$Gene.group.name, split = "|", fixed = TRUE)))
tmp.split.gene.group <- unique(tmp.split.gene.group)
tmp.split.gene.group <- tmp.split.gene.group[order(tmp.split.gene.group)]
# result: human get 1415 gene groups


##### For Describing GENE - Human
	### --- 1st process ---
	##  CellChatDB has gene reference database, but count on genes participating gene pairs
	cellchat.human.genes <- unique(c(cellchat.proc.human.db$GeneName.A, cellchat.proc.human.db$GeneName.B))
	# mapping to ref
	cellchat.human.gdf <- data.frame(dummy.col = "dummy", gene = cellchat.human.genes, stringsAsFactors = FALSE)
	cellchat.human.gdf.modf <- DataPrep.RemapClustersMarkers(cellchat.human.gdf, genes.human.ref.db, if.used.inside = TRUE)
	#
	cellchat.human.genes.res <- unique(cellchat.human.gdf.modf$result[, "gene"])
	rm(cellchat.human.gdf, cellchat.human.gdf.modf, cellchat.human.genes)
	#write.table(cellchat.human.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-human-cellchat.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	### --- 2nd process ---
	# by function definition
	# 

	# get all diff Gene.group.name
	cellchat.human.gene.groups <- unique(as.character(unlist(strsplit(cellchat.human.db$geneInfo$Gene.group.name, split = "|", fixed = TRUE))))


	# option 1.1: get all Cytokine genes by InterCellDB definition
	if (FALSE) {  # classified by Gene.group.name
	# old way to extract cytokine genes, which is wrong.
	# cytokine is not actually defined by CellChatDB, InterCellDB definition is used here
		tmp.sel.gene.groups.cytokine <- c("Chemokine ligands", "Interferons", "Interleukins", "Interleukin 6 type cytokine family", 
			"Tumor necrosis factor superfamily")
		tmp.cytokine.cellchat <- as.character(unlist(lapply(tmp.sel.gene.groups.cytokine, gene.info = cellchat.human.db$geneInfo, function(x, gene.info) {
				gene.info[grep(x, gene.info$Gene.group.name, fixed = TRUE), "Symbol"]
				})))
	}
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find cytokine from ligand), but to be safe, finding cytokine in receptor is also runned
	tmp.cytokine.cellchat.m1 <- intersect(tmp.ligand.cellchat.res, tmp.cytokine.intercelldb.res)
	tmp.cytokine.cellchat.m2 <- intersect(tmp.receptor.cellchat.res, tmp.cytokine.intercelldb.res)
	# *.m2 get 1 gene only, is "CD40LG", which is in the *.m1
	tmp.cytokine.cellchat.res <- unique(c(tmp.cytokine.cellchat.m1, tmp.cytokine.cellchat.m2))
	write.table(tmp.cytokine.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-cellchat-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.2: get all Growth factor genes by InterCellDB definition
	if (FALSE) {  # classified by Gene.group.name
	# old way to extract growfactor genes, which is wrong.
	# growfactor is not actually defined by CellChatDB, InterCellDB definition is used here
		# sub.fetch
		tmp.growfactor.explore <- unique(grep("growth factor", tolower(cellchat.human.db$geneInfo$Gene.group.name), value = TRUE, fixed= TRUE))
		tmp.growfactor.explore.res <- unique(as.character(unlist(strsplit(tmp.growfactor.explore, split = "|", fixed = TRUE))))
		#
		tmp.sel.gene.groups.growfactor <- c("Transforming growth factor beta family", "Transforming growth factor beta superfamily",
			"Fibroblast growth factor family", "Heparin binding growth factor family")
		tmp.growfactor.cellchat <- as.character(unlist(lapply(tmp.sel.gene.groups.growfactor, gene.info = cellchat.human.db$geneInfo, function(x, gene.info) {
				gene.info[grep(x, gene.info$Gene.group.name, fixed = TRUE), "Symbol"]
				})))
		tmp.growfactor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.growfactor.cellchat, stringsAsFactors = FALSE)
		tmp.growfactor.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.growfactor.cmp.df, genes.human.ref.db)$result[, "gene"])
	}
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find growfactor from ligand), but to be safe, finding growfactor in receptor is also runned
	tmp.growfactor.cellchat.m1 <- intersect(tmp.ligand.cellchat.res, tmp.growfactor.intercelldb.res)
	tmp.growfactor.cellchat.m2 <- intersect(tmp.receptor.cellchat.res, tmp.growfactor.intercelldb.res)
	# *.m2 get 0 gene, all in the *.m1
	tmp.growfactor.cellchat.res <- unique(c(tmp.growfactor.cellchat.m1, tmp.growfactor.cellchat.m2))
	write.table(tmp.growfactor.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-cellchat-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)


	# option 1.3: get all Receptor genes
	if (FALSE) {  # classified by Gene.group.name  
	# old way to extract receptor genes, which is wrong.
	# receptor should be extracted by CellChatDB interaction defnition
		# sub.fetch
		tmp.receptor.explore <- unique(grep("receptor", tolower(cellchat.human.db$geneInfo$Gene.group.name), value = TRUE, fixed= TRUE))
		tmp.receptor.explore.res <- unique(as.character(unlist(strsplit(tmp.receptor.explore, split = "|", fixed = TRUE))))
		tmp.receptor.explore.res <- unique(grep("receptor", tmp.receptor.explore.res, value = TRUE, fixed = TRUE))
		# suffix receptor / receptors
		tmp.receptor.by.suffix.1 <- grep("receptors$", tmp.receptor.explore.res, value = TRUE)
		tmp.receptor.by.suffix.2 <- grep("receptor$", tmp.receptor.explore.res, value = TRUE)
		#
		tmp.receptor.rm.suffix.res <- setdiff(tmp.receptor.explore.res, union(tmp.receptor.by.suffix.1, tmp.receptor.by.suffix.2))
		# manually remove some 
		tmp.receptor.rm.inds <- c(14,15,54,56:58,60)   # ? 34, 
		tmp.receptor.rm.suffix.rescue <- tmp.receptor.rm.suffix.res[setdiff(seq_along(tmp.receptor.rm.suffix.res), tmp.receptor.rm.inds)]
		.simpleCap <- function(x) {
				paste(toupper(substring(x, 1, 1)), substring(x, 2),
					sep = "", collapse = " ")
		}
		#
		tmp.sel.gene.groups.receptor <- sapply(c(tmp.receptor.by.suffix.1, tmp.receptor.by.suffix.2, tmp.receptor.rm.suffix.rescue), USE.NAMES = FALSE, .simpleCap)
		tmp.receptor.cellchat <- as.character(unlist(lapply(tmp.sel.gene.groups.receptor, gene.info = cellchat.human.db$geneInfo, function(x, gene.info) {
				gene.info[grep(x, gene.info$Gene.group.name, fixed = TRUE), "Symbol"]
				})))
		tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.cellchat, stringsAsFactors = FALSE)
		tmp.receptor.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db)$result[, "gene"])
		write.table(tmp.receptor.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/receptor-cellchat-IT.csv", 
			quote = FALSE, row.names = FALSE, col.names = FALSE)
	}
	# sub.fetch
	tmp.receptor.cellchat.raw <- cellchat.human.db$interaction$receptor
	# split complex and gene
	tmp.receptor.cellchat.inds.is.complex <- which(tmp.receptor.cellchat.raw %in% names(complex.human.list))
	tmp.receptor.cellchat.raw.gene <- tmp.receptor.cellchat.raw[setdiff(seq_along(tmp.receptor.cellchat.raw), tmp.receptor.cellchat.inds.is.complex)]
	tmp.receptor.cellchat.raw.complex <- tmp.receptor.cellchat.raw[tmp.receptor.cellchat.inds.is.complex]
	# get complex corresponding genes
	tmp.receptor.cellchat.genes.from.complex <- as.character(unlist(lapply(tmp.receptor.cellchat.raw.complex, complex.list = complex.human.list, 
		function(x, complex.list) {
			complex.list[which(names(complex.list) == x)][[1]]
			})))
	# merge to get the result
	tmp.receptor.cellchat <- c(tmp.receptor.cellchat.raw.gene, tmp.receptor.cellchat.genes.from.complex)
	tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.cellchat, stringsAsFactors = FALSE)
	sel.receptor.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
	#write.table(tmp.receptor.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/receptor-cellchat-IT.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.4: get all Ligand genes
	if (FALSE) {  # classified by Gene.group.name
	# old way to extract ligand genes, which is wrong.
	# ligand should be extracted by CellChatDB interaction defnition
		tmp.sel.gene.groups.ligand <- c("Receptor ligands", "Wnt family", "Chemokine ligands", 
			"Ephrins", "Interferons", "Interleukins", "Interleukin 6 type cytokine family", 
			"Tumor necrosis factor superfamily", "Growth hormone family", "VEGF family", 
			"Neurotrophins", "GDNF family ligands", "R-spondins", "Neuropeptides", 
			"Tachykinin precursors", "Endothelins", "Bone morphogenetic proteins", 
			"Inhibin subunits", "Transforming growth factor beta family", "Transforming growth factor beta superfamily",
			"Fibroblast growth factor family", "Heparin binding growth factor family")
		tmp.ligand.cellchat <- as.character(unlist(lapply(tmp.sel.gene.groups.ligand, gene.info = cellchat.human.db$geneInfo, function(x, gene.info) {
				gene.info[grep(x, gene.info$Gene.group.name, fixed = TRUE), "Symbol"]
				})))
		tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.cellchat, stringsAsFactors = FALSE)
		tmp.ligand.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db)$result[, "gene"])
	}
	# sub.fetch
	tmp.ligand.cellchat.raw <- cellchat.human.db$interaction$ligand
	# split complex and gene
	tmp.ligand.cellchat.inds.is.complex <- which(tmp.ligand.cellchat.raw %in% names(complex.human.list))
	tmp.ligand.cellchat.raw.gene <- tmp.ligand.cellchat.raw[setdiff(seq_along(tmp.ligand.cellchat.raw), tmp.ligand.cellchat.inds.is.complex)]
	tmp.ligand.cellchat.raw.complex <- tmp.ligand.cellchat.raw[tmp.ligand.cellchat.inds.is.complex]
	# get complex corresponding genes
	tmp.ligand.cellchat.genes.from.complex <- as.character(unlist(lapply(tmp.ligand.cellchat.raw.complex, complex.list = complex.human.list, 
		function(x, complex.list) {
			complex.list[which(names(complex.list) == x)][[1]]
			})))
	# merge to get the result
	tmp.ligand.cellchat <- c(tmp.ligand.cellchat.raw.gene, tmp.ligand.cellchat.genes.from.complex)
	tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.cellchat, stringsAsFactors = FALSE)
	sel.ligand.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
	#write.table(tmp.ligand.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/ligand-cellchat-IT.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.5: get all Integrin
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find integrin from ligand), but to be safe, finding integrin in receptor is also runned
	tmp.integrin.cellchat.m1 <- intersect(tmp.ligand.cellchat.res, tmp.integrin.intercelldb.res)
	tmp.integrin.cellchat.m2 <- intersect(tmp.receptor.cellchat.res, tmp.integrin.intercelldb.res)
	# *.m1 get 5 genes, *.m2 has 25 genes and includes those found in *.m1
	tmp.integrin.cellchat.res <- unique(c(tmp.integrin.cellchat.m1, tmp.integrin.cellchat.m2))
	write.table(tmp.integrin.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/integrin-cellchat-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
	tmp.gpcR.cellchat.m1 <- intersect(tmp.ligand.cellchat.res, tmp.gpcR.intercelldb.res)
	tmp.gpcR.cellchat.m2 <- intersect(tmp.receptor.cellchat.res, tmp.gpcR.intercelldb.res)
	# *.m1 get 1 genes, "ADGRE5"(NOTIN m2), and all other in *.m2
	tmp.gpcR.cellchat.res <- unique(c(tmp.gpcR.cellchat.m1, tmp.gpcR.cellchat.m2))
	write.table(tmp.gpcR.cellchat.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-cellchat-IT.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)


	### --- 3rd process ---
	# by definition of subcellular locations and molecular functions
	# [NOTE] inferred from "Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling"
	# - subcellular locations -
	splits.cellchat.scloc <- c("Plasma Membrane", "ECM", "Secreted")

	# - molecular function -
	splits.cellchat.mfunc <- c("ligand", "receptor", "agonist", "antagonist", 
		"co_A_receptor", "co_I_receptor")

##### _END_ GENE - Human



##### For Describing PAIRS - Human
### process <1> 
# get remapped and aligned gene pairs
remap.cellchat.pairs.human <- Ctrl.RemapGenePairs(cellchat.proc.human.db, genes.human.ref.db, c("GeneName.A", "GeneName.B"))
# position-align
pos.cellchat.pairs.human <- list(gene.A = remap.cellchat.pairs.human$align.GeneName.A, gene.B = remap.cellchat.pairs.human$align.GeneName.B)
# pasted style
remap.cellchat.pairs.human <- paste(remap.cellchat.pairs.human$align.GeneName.A, remap.cellchat.pairs.human$align.GeneName.B, sep = kgp.cut)
#write.table(remap.cellchat.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cat-pairs-cellchat-human.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)






##### For Describing GENE - Mouse
	### --- 1st process ---
	##  CellChatDB has gene reference database, directly use it
	cellchat.mouse.genes <- unique(c(cellchat.proc.mouse.db$GeneName.A, cellchat.proc.mouse.db$GeneName.B))
	# mapping to ref
	cellchat.mouse.gdf <- data.frame(dummy.col = "dummy", gene = cellchat.mouse.genes, stringsAsFactors = FALSE)
	cellchat.mouse.gdf.modf <- DataPrep.RemapClustersMarkers(cellchat.mouse.gdf, genes.mouse.ref.db, if.used.inside = TRUE)
	#
	cellchat.mouse.genes.res <- unique(cellchat.mouse.gdf.modf$result[, "gene"])
	rm(cellchat.mouse.gdf, cellchat.mouse.gdf.modf)
	#write.table(cellchat.mouse.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-mouse-cellchat.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)


	### --- 2nd process ---
	# by function definition
	# 

	# option 1.1: get all Cytokine genes by InterCellDB definition
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find cytokine from ligand), but to be safe, finding cytokine in receptor is also runned
	tmp.mouse.cytokine.cellchat.m1 <- intersect(tmp.mouse.ligand.cellchat.res, tmp.mouse.cytokine.intercelldb.res)
	tmp.mouse.cytokine.cellchat.m2 <- intersect(tmp.mouse.receptor.cellchat.res, tmp.mouse.cytokine.intercelldb.res)
	# *.m2 get 0 gene, all in the *.m1
	tmp.mouse.cytokine.cellchat.res <- unique(c(tmp.mouse.cytokine.cellchat.m1, tmp.mouse.cytokine.cellchat.m2))
	write.table(tmp.mouse.cytokine.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/cytokine-cellchat-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.2: get all Growth factor genes by InterCellDB definition
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find growfactor from ligand), but to be safe, finding growfactor in receptor is also runned
	tmp.mouse.growfactor.cellchat.m1 <- intersect(tmp.mouse.ligand.cellchat.res, tmp.mouse.growfactor.intercelldb.res)
	tmp.mouse.growfactor.cellchat.m2 <- intersect(tmp.mouse.receptor.cellchat.res, tmp.mouse.growfactor.intercelldb.res)
	# *.m2 get 0 gene, all in the *.m1
	tmp.mouse.growfactor.cellchat.res <- unique(c(tmp.mouse.growfactor.cellchat.m1, tmp.mouse.growfactor.cellchat.m2))
	write.table(tmp.mouse.growfactor.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/growfactor-cellchat-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.3: get all Receptor genes
	# sub.fetch
	tmp.mouse.receptor.cellchat.raw <- cellchat.mouse.db$interaction$receptor
	# split complex and gene
	tmp.mouse.receptor.cellchat.inds.is.complex <- which(tmp.mouse.receptor.cellchat.raw %in% names(complex.mouse.list))
	tmp.mouse.receptor.cellchat.raw.gene <- tmp.mouse.receptor.cellchat.raw[setdiff(seq_along(tmp.mouse.receptor.cellchat.raw), tmp.mouse.receptor.cellchat.inds.is.complex)]
	tmp.mouse.receptor.cellchat.raw.complex <- tmp.mouse.receptor.cellchat.raw[tmp.mouse.receptor.cellchat.inds.is.complex]
	# get complex corresponding genes
	tmp.mouse.receptor.cellchat.genes.from.complex <- as.character(unlist(lapply(tmp.mouse.receptor.cellchat.raw.complex, complex.list = complex.mouse.list, 
		function(x, complex.list) {
			complex.list[which(names(complex.list) == x)][[1]]
			})))
	# merge to get the result
	tmp.mouse.receptor.cellchat <- c(tmp.mouse.receptor.cellchat.raw.gene, tmp.mouse.receptor.cellchat.genes.from.complex)
	tmp.mouse.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.receptor.cellchat, stringsAsFactors = FALSE)
	sel.mouse.receptor.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.receptor.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
	#write.table(tmp.mouse.receptor.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/receptor-cellchat-mouse.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.4: get all Ligand genes
	# sub.fetch
	tmp.mouse.ligand.cellchat.raw <- cellchat.mouse.db$interaction$ligand
	# split complex and gene
	tmp.mouse.ligand.cellchat.inds.is.complex <- which(tmp.mouse.ligand.cellchat.raw %in% names(complex.mouse.list))
	tmp.mouse.ligand.cellchat.raw.gene <- tmp.mouse.ligand.cellchat.raw[setdiff(seq_along(tmp.mouse.ligand.cellchat.raw), tmp.mouse.ligand.cellchat.inds.is.complex)]
	tmp.mouse.ligand.cellchat.raw.complex <- tmp.mouse.ligand.cellchat.raw[tmp.mouse.ligand.cellchat.inds.is.complex]
	# get complex corresponding genes
	tmp.mouse.ligand.cellchat.genes.from.complex <- as.character(unlist(lapply(tmp.mouse.ligand.cellchat.raw.complex, complex.list = complex.mouse.list, 
		function(x, complex.list) {
			complex.list[which(names(complex.list) == x)][[1]]
			})))
	# merge to get the result
	tmp.mouse.ligand.cellchat <- c(tmp.mouse.ligand.cellchat.raw.gene, tmp.mouse.ligand.cellchat.genes.from.complex)
	tmp.mouse.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.ligand.cellchat, stringsAsFactors = FALSE)
	sel.mouse.ligand.cellchat.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.ligand.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
	#write.table(tmp.mouse.ligand.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/ligand-cellchat-mouse.csv", 
	#	quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.5: get all Integrin
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find integrin from ligand), but to be safe, finding integrin in receptor is also runned
	tmp.mouse.integrin.cellchat.m1 <- intersect(tmp.mouse.ligand.cellchat.res, tmp.mouse.integrin.intercelldb.res)
	tmp.mouse.integrin.cellchat.m2 <- intersect(tmp.mouse.receptor.cellchat.res, tmp.mouse.integrin.intercelldb.res)
	# *.m1 get 7 genes, *.m2 has 25 genes and includes those found in *.m1
	tmp.mouse.integrin.cellchat.res <- unique(c(tmp.mouse.integrin.cellchat.m1, tmp.mouse.integrin.cellchat.m2))
	write.table(tmp.mouse.integrin.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/integrin-cellchat-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

	# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
	# [NOTE] need after InterCellDB runs out
	# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
	tmp.mouse.gpcR.cellchat.m1 <- intersect(tmp.mouse.ligand.cellchat.res, tmp.mouse.gpcR.intercelldb.res)
	tmp.mouse.gpcR.cellchat.m2 <- intersect(tmp.mouse.receptor.cellchat.res, tmp.mouse.gpcR.intercelldb.res)
	# *.m1 get 1 genes, "Adgre5"(NOTIN m2), and all other in *.m2
	tmp.mouse.gpcR.cellchat.res <- unique(c(tmp.mouse.gpcR.cellchat.m1, tmp.mouse.gpcR.cellchat.m2))
	write.table(tmp.mouse.gpcR.cellchat.res, file = "../anal-peer-db-Result/GeneFunc-Mouse/gpcR-cellchat-mouse.csv", 
		quote = FALSE, row.names = FALSE, col.names = FALSE)

##### _END_ GENE - Mouse


##### For Describing PAIRS - Mouse
### process <1> 
# get remapped and aligned gene pairs
remap.cellchat.pairs.mouse <- Ctrl.RemapGenePairs(cellchat.proc.mouse.db, genes.mouse.ref.db, c("GeneName.A", "GeneName.B"))
# position-align
pos.cellchat.pairs.mouse <- list(gene.A = remap.cellchat.pairs.mouse$align.GeneName.A, gene.B = remap.cellchat.pairs.mouse$align.GeneName.B)
# pasted style
remap.cellchat.pairs.mouse <- paste(remap.cellchat.pairs.mouse$align.GeneName.A, remap.cellchat.pairs.mouse$align.GeneName.B, sep = kgp.cut)
#write.table(remap.cellchat.pairs.mouse, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/cat-pairs-cellchat-mouse.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

