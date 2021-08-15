
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database

library(nichenetr)

## get 3 database
# lr_network
# sig_network
# gr_network
# get mouse mapping lr_network
if (FALSE) {
	tmp.m.from <- convert_human_to_mouse_symbols(lr_network$from)
	tmp.m.to <- convert_human_to_mouse_symbols(lr_network$to)
	tmp.mouse.db <- data.frame(from = tmp.m.from, to = tmp.m.to, stringsAsFactors = FALSE)
	# saveRDS
}
nichenet.mouse.db <- readRDS("../peer-db-packed/NicheNet-mouse-db-20210303.rds")
### panel for extra informations
# 0. human only, but can use orthologs according to Ensembl to get mouse database
# 1. give the source of gene pairs, like Ramilowski, kegg, etc.


##### For Describing GENE
### --- 1st process ---
##  gene extraction
nichenet.genes <- unique(c(lr_network$from, lr_network$to))
# mapping to ref
nichenet.gdf <- data.frame(dummy.col = "dummy", gene = nichenet.genes, stringsAsFactors = FALSE)
nichenet.gdf.modf <- DataPrep.RemapClustersMarkers(nichenet.gdf, genes.human.ref.db, if.used.inside = TRUE)
#
nichenet.genes.res <- unique(nichenet.gdf.modf$result[, "gene"])
rm(nichenet.gdf, nichenet.gdf.modf, nichenet.genes)
#write.table(nichenet.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-nichenet.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)


### --- 2nd process ---
# by function definition
# [NOTE] NicheNet don't give gene classification about their function

# option 1.1: get all Cytokine genes by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find cytokine from ligand), but to be safe, finding cytokine in receptor is also runned
tmp.cytokine.nichenet.m1 <- intersect(tmp.ligand.nichenet.res, tmp.cytokine.intercelldb.res)
tmp.cytokine.nichenet.m2 <- intersect(tmp.receptor.nichenet.res, tmp.cytokine.intercelldb.res)
# *.m2 get 1 gene only, is "CD40LG", which is in the *.m1
tmp.cytokine.nichenet.res <- unique(c(tmp.cytokine.nichenet.m1, tmp.cytokine.nichenet.m2))
write.table(tmp.cytokine.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-nichenet-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.2: get all Growth factor genes by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find growfactor from ligand), but to be safe, finding growfactor in receptor is also runned
tmp.growfactor.nichenet.m1 <- intersect(tmp.ligand.nichenet.res, tmp.growfactor.intercelldb.res)
tmp.growfactor.nichenet.m2 <- intersect(tmp.receptor.nichenet.res, tmp.growfactor.intercelldb.res)
# *.m2 get 0 gene, all in *.m1
tmp.growfactor.nichenet.res <- unique(c(tmp.growfactor.nichenet.m1, tmp.growfactor.nichenet.m2))
write.table(tmp.growfactor.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-nichenet-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)


# option 1.3: get all Receptor by this defition
tmp.receptor.nichenet <- unique(lr_network$to)
tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.nichenet, stringsAsFactors = FALSE)
sel.receptor.nichenet.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.receptor.nichenet, tmp.receptor.cmp.df)
#write.table(tmp.receptor.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/receptor-nichenet-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.4: get all Ligand by this definition
tmp.ligand.nichenet <- unique(lr_network$from)
tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.nichenet, stringsAsFactors = FALSE)
sel.ligand.nichenet.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.ligand.nichenet, tmp.ligand.cmp.df)
#write.table(tmp.ligand.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/ligand-nichenet-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.5: get all Integrin by InterCellDB definition
tmp.integrin.nichenet.m1 <- intersect(tmp.ligand.nichenet.res, tmp.integrin.intercelldb.res)
tmp.integrin.nichenet.m2 <- intersect(tmp.receptor.nichenet.res, tmp.integrin.intercelldb.res)
# *.m2 get 26 genes, *.m1 get 7 genes
tmp.integrin.nichenet.res <- unique(c(tmp.integrin.nichenet.m1, tmp.integrin.nichenet.m2))
write.table(tmp.integrin.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/integrin-nichenet-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
tmp.gpcR.nichenet.m1 <- intersect(tmp.ligand.nichenet.res, tmp.gpcR.intercelldb.res)
tmp.gpcR.nichenet.m2 <- intersect(tmp.receptor.nichenet.res, tmp.gpcR.intercelldb.res)
# *.m1 get 3 genes, "CELSR1"(NOTIN m2), "CELSR2"(IN m2), "CELSR3"(IN m2), and all other in *.m2
tmp.gpcR.nichenet.res <- unique(c(tmp.gpcR.nichenet.m1, tmp.gpcR.nichenet.m2))
write.table(tmp.gpcR.nichenet.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-nichenet-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)


### --- 3rd process ---
# by definition of subcellular locations and molecular functions
# - subcellular locations -
splits.nichenet.scloc <- character(0)

# - molecular function -
splits.nichenet.mfunc <- c("Ligand", "Receptor")  # at least get 2 splits



##### For Describing PAIRS
### process <1> 
# get remapped and aligned gene pairs
remap.nichenet.pairs.human <- Ctrl.RemapGenePairs(lr_network, genes.human.ref.db, c("from", "to"))
# position-align
pos.nichenet.pairs.human <- list(gene.A = remap.nichenet.pairs.human$align.GeneName.A, gene.B = remap.nichenet.pairs.human$align.GeneName.B)
# pasted style
remap.nichenet.pairs.human <- paste(remap.nichenet.pairs.human$align.GeneName.A, remap.nichenet.pairs.human$align.GeneName.B, sep = kgp.cut)
#write.table(remap.nichenet.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cat-pairs-nichenet-human.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)




##### For Describing GENE - Mouse
### --- 1st process ---
##  gene extraction
nichenet.mouse.genes <- unique(c(nichenet.mouse.db$from, nichenet.mouse.db$to))
nichenet.mouse.genes <- nichenet.mouse.genes[which(!is.na(nichenet.mouse.genes))]
# mapping to ref
nichenet.mouse.gdf <- data.frame(dummy.col = "dummy", gene = nichenet.mouse.genes, stringsAsFactors = FALSE)
nichenet.mouse.gdf.modf <- DataPrep.RemapClustersMarkers(nichenet.mouse.gdf, genes.mouse.ref.db, if.used.inside = TRUE)
#
nichenet.mouse.genes.res <- unique(nichenet.mouse.gdf.modf$result[, "gene"])
rm(nichenet.mouse.gdf, nichenet.mouse.gdf.modf)
#write.table(nichenet.mouse.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-mouse-nichenet.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)


# Mouse MF
# option 1.3: get all Receptor by this defition
tmp.mouse.receptor.nichenet <- unique(nichenet.mouse.db$to)
tmp.mouse.receptor.nichenet <- tmp.mouse.receptor.nichenet[which(!is.na(tmp.mouse.receptor.nichenet))]
tmp.mouse.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.receptor.nichenet, stringsAsFactors = FALSE)
sel.mouse.receptor.nichenet.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.receptor.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.mouse.receptor.nichenet, tmp.mouse.receptor.cmp.df)

# option 1.4: get all Ligand by this definition
tmp.mouse.ligand.nichenet <- unique(nichenet.mouse.db$from)
tmp.mouse.ligand.nichenet <- tmp.mouse.ligand.nichenet[which(!is.na(tmp.mouse.ligand.nichenet))]
tmp.mouse.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.ligand.nichenet, stringsAsFactors = FALSE)
sel.mouse.ligand.nichenet.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.ligand.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.mouse.ligand.nichenet, tmp.mouse.ligand.cmp.df)



##### For Describing GENE PAIRS - Mouse
# remove NAs
remap.nichenet.ref.db <- nichenet.mouse.db[intersect(which(!is.na(nichenet.mouse.db$from)), which(!is.na(nichenet.mouse.db$to))), ]
remap.nichenet.pairs.mouse <- Ctrl.RemapGenePairs(remap.nichenet.ref.db, genes.mouse.ref.db, c("from", "to"))
# position-align
pos.nichenet.pairs.mouse <- list(gene.A = remap.nichenet.pairs.mouse$align.GeneName.A, gene.B = remap.nichenet.pairs.mouse$align.GeneName.B)
# pasted style
remap.nichenet.pairs.mouse <- paste(remap.nichenet.pairs.mouse$align.GeneName.A, remap.nichenet.pairs.mouse$align.GeneName.B, sep = kgp.cut)
rm(remap.nichenet.ref.db)

