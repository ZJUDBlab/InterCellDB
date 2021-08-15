
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database

# fetch database
italk.db <- readRDS("../peer-db-packed/iTalk-lrdb-20210124.rds")
### panel for extra informations
# 0. human only
# 1. get classification of 4 ligand function
#
# warn1: database has some NA values in Receptor.ApprovedSymbol, be careful!
# warn2: CCL2 is not existing anymore by HGNC. iTALK records it anyway
#



##### For Describing GENE
### --- 1st process ---
##  gene extraction
italk.genes <- c(italk.db$Ligand.ApprovedSymbol, italk.db$Receptor.ApprovedSymbol)
# get some Receptor.ApprovedSymbol NA
italk.genes <- c(italk.genes, italk.db[which(is.na(italk.db$Receptor.ApprovedSymbol)), "Receptor.Name"])
italk.genes <- unique(italk.genes[which(!is.na(italk.genes))])
# mapping to ref
italk.gdf <- data.frame(dummy.col = "dummy", gene = italk.genes, stringsAsFactors = FALSE)
italk.gdf.modf <- DataPrep.RemapClustersMarkers(italk.gdf, genes.human.ref.db, if.used.inside = TRUE)
#
italk.genes.res <- unique(italk.gdf.modf$result[, "gene"])
rm(italk.gdf, italk.gdf.modf, italk.genes)
#write.table(italk.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-italk.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

### --- 2nd process ---
# by function definition
# ligand
#   checkpoint      cytokine    growth factor         other 
#           32           327              227          2063 
# receptor: NO diff on it.
#

# option 1.1: get all Cytokine genes
tmp.cytokine.italk <- unique(italk.db[which(italk.db$Classification == "cytokine"), "Ligand.ApprovedSymbol"])
tmp.cytokine.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.cytokine.italk, stringsAsFactors = FALSE)
tmp.cytokine.italk.res <- unique(DataPrep.RemapClustersMarkers(tmp.cytokine.cmp.df, genes.human.ref.db)$result[, "gene"])
write.table(tmp.cytokine.italk.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-italk-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.2: get all Growth factor genes
tmp.growfactor.italk <- unique(italk.db[which(italk.db$Classification == "growth factor"), "Ligand.ApprovedSymbol"])
tmp.growfactor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.growfactor.italk, stringsAsFactors = FALSE)
tmp.growfactor.italk.res <- unique(DataPrep.RemapClustersMarkers(tmp.growfactor.cmp.df, genes.human.ref.db)$result[, "gene"])
write.table(tmp.growfactor.italk.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-italk-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.3: get all Receptor, iTalk treats every pair to be ligand-receptor
tmp.receptor.italk <- unique(italk.db$Receptor.ApprovedSymbol)
# some get no approved symbol
tmp.receptor.italk <- tmp.receptor.italk[which(!is.na(tmp.receptor.italk))]
tmp.receptor.italk <- c(tmp.receptor.italk, unique(italk.db[which(is.na(italk.db$Receptor.ApprovedSymbol)), "Receptor.Name"]))
# doing on merged result 
tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.italk, stringsAsFactors = FALSE)
sel.receptor.italk.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.receptor.italk, tmp.receptor.cmp.df)
#write.table(tmp.receptor.italk.res, file = "../anal-peer-db-Result/GeneFunc/receptor-italk-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.4: get all Ligand, iTalk treats every pair to be ligand-receptor
tmp.ligand.italk <- unique(italk.db$Ligand.ApprovedSymbol)
# some get no approved symbol
tmp.ligand.italk <- tmp.ligand.italk[which(!is.na(tmp.ligand.italk))]
tmp.ligand.italk <- c(tmp.ligand.italk, unique(italk.db[which(is.na(italk.db$Ligand.ApprovedSymbol)), "Ligand.Name"]))
# doing on merged result 
tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.italk, stringsAsFactors = FALSE)
sel.ligand.italk.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.ligand.italk, tmp.ligand.cmp.df)
#write.table(tmp.ligand.italk.res, file = "../anal-peer-db-Result/GeneFunc/ligand-italk-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.5: get all Integrin by InterCellDB definition
tmp.integrin.italk.m1 <- intersect(tmp.ligand.italk.res, tmp.integrin.intercelldb.res)
tmp.integrin.italk.m2 <- intersect(tmp.receptor.italk.res, tmp.integrin.intercelldb.res)
# *.m1 get 0 genes, *.m2 get 26 genes
tmp.integrin.italk.res <- unique(c(tmp.integrin.italk.m1, tmp.integrin.italk.m2))
write.table(tmp.integrin.italk.res, file = "../anal-peer-db-Result/GeneFunc/integrin-italk-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
tmp.gpcR.italk.m1 <- intersect(tmp.ligand.italk.res, tmp.gpcR.intercelldb.res)
tmp.gpcR.italk.m2 <- intersect(tmp.receptor.italk.res, tmp.gpcR.intercelldb.res)
# *.m1 get 0 genes, and all in *.m2
tmp.gpcR.italk.res <- unique(c(tmp.gpcR.italk.m1, tmp.gpcR.italk.m2))
write.table(tmp.gpcR.italk.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-italk-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)


### --- 3rd process ---
# by definition of subcellular locations and molecular functions
# - subcellular locations -
splits.italk.scloc <- character(0)

# - molecular function -
splits.italk.mfunc <- c("checkpoint", "cytokine", "growth factor", "other",  # ligand splits to 4 categories
	"receptor")





##### For Describing PAIRS
### process <1> 
# get remapped and aligned gene pairs
tmp.italk.db <- italk.db
tmp.inds.italk.rec.NA <- which(is.na(tmp.italk.db$Receptor.ApprovedSymbol))
tmp.italk.db[tmp.inds.italk.rec.NA, "Receptor.ApprovedSymbol"] <- tmp.italk.db[tmp.inds.italk.rec.NA, "Receptor.Name"]
remap.italk.pairs.human <- Ctrl.RemapGenePairs(tmp.italk.db, genes.human.ref.db, c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol"))
# position-align
pos.italk.pairs.human <- list(gene.A = remap.italk.pairs.human$align.GeneName.A, gene.B = remap.italk.pairs.human$align.GeneName.B)
# pasted style
remap.italk.pairs.human <- paste(remap.italk.pairs.human$align.GeneName.A, remap.italk.pairs.human$align.GeneName.B, sep = kgp.cut)
rm(tmp.italk.db, tmp.inds.italk.rec.NA)
#write.table(remap.italk.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cat-pairs-italk-human.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)








