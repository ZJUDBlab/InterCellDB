
# created at 2021.01.24

library(InterCellDB)  # use it to mapping to ref gene database


# get result 
cellphone.db <- readRDS("../peer-db-packed/CellPhoneDB-gene-pairs.rds")
### panel for extra informations
# 0. human only. It said to use orthologs to get mouse database(Method not shown)
# 1. give many attribute about the complex but not the gene
#	"transmembrane", "peripheral", "secreted", "secreted_desc", "secreted_highlight", 
#	"receptor", "receptor_desc", "integrin", "other", "other_desc")
# 2. it gives interaction pairs PUBMED IDs or PMC <url> or just from what database 
#  2.ex For exploring its source and annotation_strategy, it can be devided into 2 kinds
#       From IMEx (annotation_strategy != "curated" plus "curated" & Uniprot) 
#       From curated ("curated" and not "uniprot")
#  But, finally it still the Experimentally validated pairs.
#
# warn1:
#


##### For Describing GENE
### --- 1st process ---
##  gene extraction
cellphone.genes <- unique(c(cellphone.db$genename.A, cellphone.db$genename.B))
# mapping to ref
cellphone.gdf <- data.frame(dummy.col = "dummy", gene = cellphone.genes, stringsAsFactors = FALSE)
cellphone.gdf.modf <- DataPrep.RemapClustersMarkers(cellphone.gdf, genes.human.ref.db, if.used.inside = TRUE)
#
cellphone.genes.res <- unique(cellphone.gdf.modf$result[, "gene"])
rm(cellphone.gdf, cellphone.gdf.modf, cellphone.genes)
#write.table(cellphone.genes.res, file = "../anal-peer-db-Result/GeneItself/allgene-cellphone.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)


### --- 2nd process ---
# by function definition
# ligand - get from tapply on protein.cpdb$secreted_desc
#
#            	 <empty>          CellSignal_WNT                cytokine 
#                    651                      11                       4 
#               Cytokine      Cytokine | Hormone           Cytokine_like 
#                    128                       1                       1 
#           Growthfactor Growthfactor | Cytokine  Growthfactor | Hormone 
#                     65                      43                       6 
#                Hormone          Immune-related 
#                     67                       1 
# Get to be Cytokine, Growthfactor, Hormone. [IN Protein level]
#              177         114        74
# receptor
#    get 2 main definitions: Receptor and Integrin
tmp.protein.cpdb <- read.csv("/Users/jinziyang/Sc-RNAsequence/interaction-database/Ver2-prep-pics/R-v2-pic/extdata/CellPhoneDB/protein_input.csv", stringsAsFactors = FALSE)
tmp.gene.cpdb <- read.csv("/Users/jinziyang/Sc-RNAsequence/interaction-database/Ver2-prep-pics/R-v2-pic/extdata/CellPhoneDB/gene_input.csv", stringsAsFactors = FALSE)
tmp.complex.cpdb <- read.csv("/Users/jinziyang/Sc-RNAsequence/interaction-database/Ver2-prep-pics/R-v2-pic/extdata/CellPhoneDB/complex_input.csv", stringsAsFactors = FALSE)

# generate protein~complex reference for complex [Code copied from CellPhoneDB-explore.R]
cpdb.attr.all.names <- c("transmembrane", "peripheral", 
	"secreted", "secreted_desc", "secreted_highlight", 
	"receptor", "receptor_desc", 
	"integrin", "other", "other_desc")
tmp.uniprot.cols <- grep("^uniprot_", colnames(tmp.complex.cpdb), value = TRUE)
tmp.df.list <- list()
for (i in 1:nrow(tmp.complex.cpdb)) {
	this.row.info <- tmp.complex.cpdb[i, c("complex_name", cpdb.attr.all.names)]
	this.unip.s <- as.character(tmp.complex.cpdb[i, tmp.uniprot.cols])  # this make NA to be "NA"
	this.unip.s <- this.unip.s[intersect(which(this.unip.s != ""), which(this.unip.s != "NA"))]
	tmp.df <- data.frame(uniprot = this.unip.s, complex_name = tmp.complex.cpdb[i, "complex_name"], stringsAsFactors = FALSE)
	this.df <- left_join(tmp.df, this.row.info)
	tmp.df.list <- c(tmp.df.list, list(this.df))
}
tmp.complex.cpdb.splits <- bind_rows(tmp.df.list)
rm(tmp.df.list)

# option 1.1: get all Cytokine genes by this definition
tmp.cytokine.desc <- c("cytokine", "Cytokine", "Cytokine | Hormone", "Cytokine_like", "Growthfactor | Cytokine")
tmp.cytokine.uniprotID <- tmp.protein.cpdb[which(tmp.protein.cpdb$secreted_desc %in% tmp.cytokine.desc), "uniprot"]
tmp.cytokine.cpdb <- unique(tmp.gene.cpdb[which(tmp.gene.cpdb$uniprot %in% tmp.cytokine.uniprotID), "gene_name"])
tmp.cytokine.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.cytokine.cpdb, stringsAsFactors = FALSE)
tmp.cytokine.cpdb.res <- unique(DataPrep.RemapClustersMarkers(tmp.cytokine.cmp.df, genes.human.ref.db)$result[, "gene"])
write.table(tmp.cytokine.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/cytokine-cellphone-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.2: get all Growth factor by this definition
tmp.growfactor.desc <- c("Growthfactor", "Growthfactor | Cytokine", "Growthfactor | Hormone")
tmp.growfactor.uniprotID <- tmp.protein.cpdb[which(tmp.protein.cpdb$secreted_desc %in% tmp.growfactor.desc), "uniprot"]
tmp.growfactor.cpdb <- unique(tmp.gene.cpdb[which(tmp.gene.cpdb$uniprot %in% tmp.growfactor.uniprotID), "gene_name"])
tmp.growfactor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.growfactor.cpdb, stringsAsFactors = FALSE)
tmp.growfactor.cpdb.res <- unique(DataPrep.RemapClustersMarkers(tmp.growfactor.cmp.df, genes.human.ref.db)$result[, "gene"])
write.table(tmp.growfactor.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/growfactor-cellphone-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.3: get all Receptor by this defition
tmp.receptor.cols <- c("receptor", "integrin")  # this 2 cols = TRUE, represents the receptor
tmp.receptor.uniprotID.from.pro <- tmp.protein.cpdb[union(which(tmp.protein.cpdb$receptor == "True"), which(tmp.protein.cpdb$integrin == "True")), "uniprot"]
tmp.receptor.uniprotID.from.complex <- tmp.complex.cpdb.splits[union(which(tmp.complex.cpdb.splits$receptor == "True"), which(tmp.complex.cpdb.splits$integrin == "True")), "uniprot"]
tmp.receptor.uniprotID <- c(tmp.receptor.uniprotID.from.pro, tmp.receptor.uniprotID.from.complex)
tmp.receptor.cpdb <- unique(tmp.gene.cpdb[which(tmp.gene.cpdb$uniprot %in% tmp.receptor.uniprotID), "gene_name"])
tmp.receptor.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.receptor.cpdb, stringsAsFactors = FALSE)
sel.receptor.cpdb.res <- unique(DataPrep.RemapClustersMarkers(tmp.receptor.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.receptor.cols, tmp.receptor.uniprotID.from.pro, tmp.receptor.uniprotID.from.complex, tmp.receptor.uniprotID, tmp.receptor.cpdb, tmp.receptor.cmp.df)
#write.table(tmp.receptor.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/receptor-cellphone-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.4: get all Ligand by this definition
tmp.ligand.cols <- c("secreted")  # ligand defined by their subcellular location
tmp.ligand.uniprotID.from.pro <- tmp.protein.cpdb[which(tmp.protein.cpdb$secreted == "True"), "uniprot"]
tmp.ligand.uniprotID.from.complex <- tmp.complex.cpdb.splits[which(tmp.complex.cpdb.splits$secreted == "True"), "uniprot"]
tmp.ligand.uniprotID <- c(tmp.ligand.uniprotID.from.pro, tmp.ligand.uniprotID.from.complex)
tmp.ligand.cpdb <- unique(tmp.gene.cpdb[which(tmp.gene.cpdb$uniprot %in% tmp.ligand.uniprotID), "gene_name"])
tmp.ligand.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.ligand.cpdb, stringsAsFactors = FALSE)
sel.ligand.cpdb.res <- unique(DataPrep.RemapClustersMarkers(tmp.ligand.cmp.df, genes.human.ref.db, if.used.inside = TRUE)$result[, "gene"])
rm(tmp.ligand.cols, tmp.ligand.uniprotID.from.pro, tmp.ligand.uniprotID.from.complex, tmp.ligand.uniprotID, tmp.ligand.cpdb, tmp.ligand.cmp.df)
#write.table(tmp.ligand.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/ligand-cellphone-IT.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.5: get all Integrin by this definition
tmp.integrin.cols <- c("integrin")  # this 2 cols = TRUE, represents the integrin
tmp.integrin.uniprotID.from.pro <- tmp.protein.cpdb[which(tmp.protein.cpdb$integrin == "True"), "uniprot"]
tmp.integrin.uniprotID.from.complex <- tmp.complex.cpdb.splits[which(tmp.complex.cpdb.splits$integrin == "True"), "uniprot"]
tmp.integrin.uniprotID <- c(tmp.integrin.uniprotID.from.pro, tmp.integrin.uniprotID.from.complex)
tmp.integrin.cpdb <- unique(tmp.gene.cpdb[which(tmp.gene.cpdb$uniprot %in% tmp.integrin.uniprotID), "gene_name"])
tmp.integrin.cmp.df <- data.frame(dummy.col = "dummy", gene = tmp.integrin.cpdb, stringsAsFactors = FALSE)
tmp.integrin.cpdb.res <- unique(DataPrep.RemapClustersMarkers(tmp.integrin.cmp.df, genes.human.ref.db)$result[, "gene"])
write.table(tmp.integrin.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/integrin-cellphone-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# option 1.6: get all G-protein Coupled Receptor by InterCellDB definition
# [NOTE] need after InterCellDB runs out
# The normal routine is to get *.m1(find gpcR from ligand), but to be safe, finding gpcR in receptor is also runned
tmp.gpcR.cpdb.m1 <- intersect(tmp.ligand.cpdb.res, tmp.gpcR.intercelldb.res)
tmp.gpcR.cpdb.m2 <- intersect(tmp.receptor.cpdb.res, tmp.gpcR.intercelldb.res)
# *.m1 get 1 genes, "ADGRE5"(NOTIN m2), and all other in *.m2
tmp.gpcR.cpdb.res <- unique(c(tmp.gpcR.cpdb.m1, tmp.gpcR.cpdb.m2))
write.table(tmp.gpcR.cpdb.res, file = "../anal-peer-db-Result/GeneFunc/gpcR-cpdb-IT.csv", 
	quote = FALSE, row.names = FALSE, col.names = FALSE)



### --- 3rd process ---
# by definition of subcellular locations and molecular functions
# - subcellular locations -
splits.cellphone.scloc <- c("transmembrane", "peripheral", "secreted")

# - molecular function -
splits.cellphone.mfunc <- c("Receptor", "Integrin", "Cytokine", "Growthfactor", "Hormone",
	"Immune-related", "CellSignal_WNT")
	# [NOTE] "Immune-related" get 1 record, while "CellSignal_WNT" get 11 records
	# those are consider not used as splits?





##### For Describing PAIRS
### process <1> 
# get remapped and aligned gene pairs
remap.cellphone.pairs.human <- Ctrl.RemapGenePairs(cellphone.db, genes.human.ref.db, c("genename.A", "genename.B"))
# position-align
pos.cellphone.pairs.human <- list(gene.A = remap.cellphone.pairs.human$align.GeneName.A, gene.B = remap.cellphone.pairs.human$align.GeneName.B)
# pasted style
remap.cellphone.pairs.human <- paste(remap.cellphone.pairs.human$align.GeneName.A, remap.cellphone.pairs.human$align.GeneName.B, sep = kgp.cut)
#write.table(remap.cellphone.pairs.human, file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cat-pairs-cellphone-human.csv", 
#	quote = FALSE, row.names = FALSE, col.names = FALSE)


