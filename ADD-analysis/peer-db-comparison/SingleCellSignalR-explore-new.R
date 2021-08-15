
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database

# fetch database
scsr.db <- readRDS("../peer-db-packed/SingleCellSignalR-LRdb-20210124.rds")
# use orthologs (SingleCellSignalR::mm2Hs) to generate mouse db, and the result is put here
if (FALSE) {  # [NOTE] the code is runned under R 4.0.4
	scsr.mouse.db <- left_join(scsr.db, mm2Hs, by = c("ligand" = "Gene name"))
	scsr.mouse.db <- scsr.mouse.db[, -1]
	colnames(scsr.mouse.db)[4] <- "ligand"
	scsr.mouse.db <- left_join(scsr.mouse.db, mm2Hs, by = c("receptor" = "Gene name"))
	scsr.mouse.db <- scsr.mouse.db[, -1]
	colnames(scsr.mouse.db)[4] <- "receptor"
	scsr.mouse.db <- scsr.mouse.db[, c(3, 4, 1, 2)]
}
scsr.mouse.db <- readRDS("../peer-db-packed/SingleCellSignalR-LRdb-mouse-20210220.rds")
### panel for extra informations
# 0. human only, but can use orthologs according to Ensembl to get mouse database
# 1. give many pairs PMID reference, with 2013 pairs having but 1238 not
#
# warn1:
#


##### For Describing GENE - Human
### --- 1st process ---
##  gene extraction
scsr.genes <- unique(c(scsr.db$ligand, scsr.db$receptor))
# mapping to ref
scsr.gdf <- data.frame(dummy.col = "dummy", gene = scsr.genes, stringsAsFactors = FALSE)
scsr.gdf.modf <- DataPrep.RemapClustersMarkers(scsr.gdf, genes.human.ref.db, if.used.inside = TRUE)
#
scsr.genes.res <- unique(scsr.gdf.modf$result[, "gene"])
rm(scsr.gdf, scsr.gdf.modf, scsr.genes)
#write.table(scsr.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-scsr.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)


### --- 2nd process ---
# by function definition
# [NOTE] No classification about gene function

# option 1.1: get all Cytokine genes by InterCellDB definition
# [NOTE] need after InterCellDB runs out, as well as Receptor remapped from SingleCellSignalR
# The normal routine is to get *.m1(find cytokine from ligand), but to be safe, finding cytokine in receptor is also runned
tmp.cytokine.scsr.m1 <- intersect(tmp.ligand.scsr.res, tmp.cytokine.intercelldb.res)
tmp.cytokine.scsr.m2 <- intersect(tmp.receptor.scsr.res, tmp.cytokine.intercelldb.res)
# *.m2 get 1 gene only, is "IL36RN", which is NOT in the *.m1
tmp.cytokine.scsr.res <- unique(c(tmp.cytokine.scsr.m1, tmp.cytokine.scsr.m2))
write.table(tmp.cytokine.scsr.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-scsr-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.2: get all Growth factor genes by InterCellDB definition
# The normal routine is to get *.m1(find cytokine from ligand), but to be safe, finding cytokine in receptor is also runned
tmp.growfactor.scsr.m1 <- intersect(tmp.ligand.scsr.res, tmp.growfactor.intercelldb.res)
tmp.growfactor.scsr.m2 <- intersect(tmp.receptor.scsr.res, tmp.growfactor.intercelldb.res)
# *.m2 get 0 gene, all in *.m1
tmp.growfactor.scsr.res <- unique(c(tmp.growfactor.scsr.m1, tmp.growfactor.scsr.m2))
write.table(tmp.growfactor.scsr.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-scsr-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)


# option 1.3: get all Receptor genes by this definition
tmp.receptor.scsr <- unique(scsr.db$receptor)
tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.scsr, stringsAsFactors = FALSE)
sel.receptor.scsr.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.receptor.scsr, tmp.receptor.cmp.df)
#write.table(tmp.receptor.scsr.res, file = "../anal-peer-db-Result/GeneFunc/receptor-scsr-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.4: get all Ligand genes by this definition
tmp.ligand.scsr <- unique(scsr.db$ligand)
tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.scsr, stringsAsFactors = FALSE)
sel.ligand.scsr.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.ligand.scsr, tmp.ligand.cmp.df)
#write.table(tmp.ligand.scsr.res, file = "../anal-peer-db-Result/GeneFunc/ligand-scsr-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.5: get all Integrin by InterCellDB definition
tmp.integrin.scsr.m1 <- intersect(tmp.ligand.scsr.res, tmp.integrin.intercelldb.res)
tmp.integrin.scsr.m2 <- intersect(tmp.receptor.scsr.res, tmp.integrin.intercelldb.res)
# *.m1 get 1 gene only, which is "ITGB2", and *.m2 get it included
tmp.integrin.scsr.res <- unique(c(tmp.integrin.scsr.m1, tmp.integrin.scsr.m2))
write.table(tmp.integrin.scsr.res, file = "../anal-peer-db-Result/GeneFunc/integrin-scsr-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
tmp.gpcR.scsr.m1 <- intersect(tmp.ligand.scsr.res, tmp.gpcR.intercelldb.res)
tmp.gpcR.scsr.m2 <- intersect(tmp.receptor.scsr.res, tmp.gpcR.intercelldb.res)
# *.m1 get 0 genes, and all in *.m2
tmp.gpcR.scsr.res <- unique(c(tmp.gpcR.scsr.m1, tmp.gpcR.scsr.m2))
write.table(tmp.gpcR.scsr.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-scsr-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)



### --- 3rd process ---
# by definition of subcellular locations and molecular functions
# [NOTE]
# It defines ligand from GO:0005615 [extracellular space] and GO:0005576 [extracellular region]
# It defines receptor from GO:0043235 [receptor complex]
# - subcellular locations -
splits.scsr.scloc <- character(0)

# - molecular function -
splits.scsr.mfunc <- c("ligand", "receptor")


##### For Describing PAIRS
### process <1> 
# get remapped and aligned gene pairs
remap.scsr.pairs.human <- Ctrl.RemapGenePairs(scsr.db, genes.human.ref.db, c("ligand", "receptor"))
# position-align
pos.scsr.pairs.human <- list(gene.A = remap.scsr.pairs.human$align.GeneName.A, gene.B = remap.scsr.pairs.human$align.GeneName.B)
# pasted style
remap.scsr.pairs.human <- paste(remap.scsr.pairs.human$align.GeneName.A, remap.scsr.pairs.human$align.GeneName.B, sep = kgp.cut)
#write.table(remap.scsr.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cat-pairs-scsr-human.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)






##### For Describing GENE - Mouse
### --- 1st process ---
##  gene extraction
scsr.mouse.genes <- unique(c(scsr.mouse.db$ligand, scsr.mouse.db$receptor))
scsr.mouse.genes <- scsr.mouse.genes[which(!is.na(scsr.mouse.genes))]
# mapping to ref
scsr.mouse.gdf <- data.frame(dummy.col = "dummy", gene = scsr.mouse.genes, stringsAsFactors = FALSE)
scsr.mouse.gdf.modf <- DataPrep.RemapClustersMarkers(scsr.mouse.gdf, genes.mouse.ref.db, if.used.inside = TRUE)
#
scsr.mouse.genes.res <- unique(scsr.mouse.gdf.modf$result[, "gene"])
rm(scsr.mouse.gdf, scsr.mouse.gdf.modf)
#write.table(scsr.mouse.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-mouse-scsr.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# Mouse MF
# option 1.3: get all Receptor genes by this definition
tmp.mouse.receptor.scsr <- unique(scsr.mouse.db$receptor)
tmp.mouse.receptor.scsr <- tmp.mouse.receptor.scsr[which(!is.na(tmp.mouse.receptor.scsr))]
tmp.mouse.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.receptor.scsr, stringsAsFactors = FALSE)
sel.mouse.receptor.scsr.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.receptor.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.mouse.receptor.scsr, tmp.mouse.receptor.cmp.df)

# option 1.4: get all Ligand genes by this definition
tmp.mouse.ligand.scsr <- unique(scsr.mouse.db$ligand)
tmp.mouse.ligand.scsr <- tmp.mouse.ligand.scsr[which(!is.na(tmp.mouse.ligand.scsr))]
tmp.mouse.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.mouse.ligand.scsr, stringsAsFactors = FALSE)
sel.mouse.ligand.scsr.res <- unique(DataPrep.RemapClustersMarkers(tmp.mouse.ligand.cmp.df, genes.mouse.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.mouse.ligand.scsr, tmp.mouse.ligand.cmp.df)


##### For Describing GENE PAIRS - Mouse
###
#
# remove NAs
remap.scsr.ref.db <- scsr.mouse.db[intersect(which(!is.na(scsr.mouse.db$ligand)), which(!is.na(scsr.mouse.db$receptor))), ]
remap.scsr.pairs.mouse <- Ctrl.RemapGenePairs(remap.scsr.ref.db, genes.mouse.ref.db, c("ligand", "receptor"))
# position-align
pos.scsr.pairs.mouse <- list(gene.A = remap.scsr.pairs.mouse$align.GeneName.A, gene.B = remap.scsr.pairs.mouse$align.GeneName.B)
# pasted style
remap.scsr.pairs.mouse <- paste(remap.scsr.pairs.mouse$align.GeneName.A, remap.scsr.pairs.mouse$align.GeneName.B, sep = kgp.cut)
rm(remap.scsr.ref.db)


