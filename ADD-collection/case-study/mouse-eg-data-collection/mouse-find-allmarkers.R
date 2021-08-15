# This file is to process the data to find all markers for subsequent analysis

library(Seurat)
library(dplyr)

# TSNE colour
colour.ref <- c("#33a3dc", "#ef5b9c", "#1d953f", "#00a6ac", "#a67c52", "#d71345", 
	"#6a6da9", "#c7b299", "#f3715c", "#FFE3B6",  "#ba9bc9", "#f15a22", "#5c7a29")  # "#fdb933"
colour.add <- c(colour.ref, "#8A110C", "#FBB719", "#232C68", "#778698")


# ensembl human ref database
ensembl.10090.db <- readRDS("/Users/jinziyang/Sc-RNAsequence/interaction-database/prepFinal-1/peer-db/res-db-packed/ensembl-10090-20201220.rds")

# transform the raw rows to standard 10X data format that could be read by Seurat v3
raw.rows <- read.table("../E-EHCA-2-quantification-raw-files/E-EHCA-2.aggregated_filtered_counts.mtx_rows", 
	sep = "\t", stringsAsFactors = FALSE)
if (length(which(raw.rows[,1] != raw.rows[,2])) != 0) {
	stop("Wrong data!")
}
colnames(raw.rows)[1:2] <- "ensemblID"
#
# trs.rows <- left_join(raw.rows[, 1, drop = FALSE], ensembl.10090.db, by = c("ensemblID" = "Gene.stable.ID"))
# Failed to use the merged Ensembl.db to process the data
# use the all versions of Ensembl mouse database to map these genes
	# collect all the unknown Ensembl IDs
	tmp.unk.target <- raw.rows[, 1]

# tmply change the workstation
setwd("/Users/jinziyang/Sc-RNAsequence/interaction-database/prepFinal-1/peer-db/R-peer-db-check")
	# formatted colnames
	std.colnames.GRCm <- c("Gene.stable.ID", "Protein.stable.ID", "Gene.name", "NCBI.gene.ID", "MGI.ID")

	# GRCmp6 realease 101 as the base
	ensembl.p6.db <- read.csv("../raw-db-set/GRCm38p6-release101-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
	colnames(ensembl.p6.db) <- std.colnames.GRCm
	ensembl.p6.db <- ensembl.p6.db[which(ensembl.p6.db$Protein.stable.ID != ""), ]
	# 115173 rows -> 69022 rows

	# GRCmp5 release 91 as the backup 1
	ensembl.p5.db <- read.csv("../raw-db-set/GRCm38p5-release91-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
	colnames(ensembl.p5.db) <- std.colnames.GRCm
	ensembl.p5.db <- ensembl.p5.db[which(ensembl.p5.db$Protein.stable.ID != ""), ]
	# 108738 rows -> 65693 rows

	# GRCmp4 release 86 as the backup 2
	ensembl.p4.db <- read.csv("../raw-db-set/GRCm38p4-release86-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
	colnames(ensembl.p4.db) <- std.colnames.GRCm
	ensembl.p4.db <- ensembl.p4.db[which(ensembl.p4.db$Protein.stable.ID != ""), ]
	# 97912 rows -> 61122 rows

	# GRCmp3 release 80 as the backup 3
	ensembl.p3.db <- read.csv("../raw-db-set/GRCm38p3-release80-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
	colnames(ensembl.p3.db) <- std.colnames.GRCm
	ensembl.p3.db <- ensembl.p3.db[which(ensembl.p3.db$Protein.stable.ID != ""), ]
	# 88998 rows -> 56550 rows

	# GRCmp2 release 77 as the backup 4
	ensembl.p2.db <- read.csv("../raw-db-set/GRCm38p2-release77-mouse.txt", header = TRUE, stringsAsFactors = FALSE)
	colnames(ensembl.p2.db) <- std.colnames.GRCm
	ensembl.p2.db <- ensembl.p2.db[which(ensembl.p2.db$Protein.stable.ID != ""), ]
	# 90368 rows -> 62269 rows

	### match process
	# 1. p6 r101 matching
	inds.match.p6.101 <- which(tmp.unk.target %in% ensembl.p6.db$Gene.stable.ID)
	# 2. p5 r91 matching
	inds.match.p5.r91 <- which(tmp.unk.target %in% ensembl.p5.db$Gene.stable.ID)
	inds.match.p5.r91 <- setdiff(inds.match.p5.r91, inds.match.p6.101)  # p6 get the highest priority
	# 3. p4 r91 matching
	inds.match.p4.r86 <- which(tmp.unk.target %in% ensembl.p4.db$Gene.stable.ID)
	inds.match.p4.r86 <- setdiff(inds.match.p4.r86, inds.match.p6.101)
	inds.match.p4.r86 <- setdiff(inds.match.p4.r86, inds.match.p5.r91)
	# 4. p3 r80 matching
	inds.match.p3.r80 <- which(tmp.unk.target %in% ensembl.p3.db$Gene.stable.ID)
	inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p6.101)
	inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p5.r91)
	inds.match.p3.r80 <- setdiff(inds.match.p3.r80, inds.match.p4.r86)
	# 5. p2 r77 matching
	inds.match.p2.r77 <- which(tmp.unk.target %in% ensembl.p2.db$Gene.stable.ID)
	inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p6.101)
	inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p5.r91)
	inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p4.r86)
	inds.match.p2.r77 <- setdiff(inds.match.p2.r77, inds.match.p3.r80)
	# collect all matched result
	mouse.ensembl.patch.list <- list(p5 = ensembl.p5.db[which(ensembl.p5.db$Gene.stable.ID %in% tmp.unk.target[inds.match.p5.r91]), ], 
		p4 = ensembl.p4.db[which(ensembl.p4.db$Gene.stable.ID %in% tmp.unk.target[inds.match.p4.r86]), ],
		p3 = ensembl.p3.db[which(ensembl.p3.db$Gene.stable.ID %in% tmp.unk.target[inds.match.p3.r80]), ],
		p2 = ensembl.p2.db[which(ensembl.p2.db$Gene.stable.ID %in% tmp.unk.target[inds.match.p2.r77]), ])
	trs.res.db <- rbind(ensembl.p6.db[which(ensembl.p6.db$Gene.stable.ID %in% tmp.unk.target[inds.match.p6.101]), ], dplyr::bind_rows(mouse.ensembl.patch.list))
	# as Protein.stable.ID is not cared
	trs.res.db <- InterCellDB::DoPartUnique(trs.res.db, c(1,3))
	## check if Gene.name -- Gene.stable.ID has duplicates in different Gene.name but the same Gene.stable.ID
	tmp.a <- tapply(1:nrow(trs.res.db), trs.res.db$Gene.stable.ID, length)
	if (length(which(tmp.a > 1)) != 0) {
		stop("Duplicates exist in Gene.stable.ID, which could map to different Gene.name!")
	}


#
# setwd back
#
setwd("/Users/jinziyang/Sc-RNAsequence/interaction-database/Peer-Work-2a/UseForData-Mouse/optionxf1-data/R-xf1")

# continue the process
trs.rows <- left_join(raw.rows[, 1, drop = FALSE], trs.res.db, by = c("ensemblID" = "Gene.stable.ID"))
## with 6614 rows has no detailed Gene.name, 
## after search manually on Ensembl website, get those part of pseudogene,
## so as to ignore these genes.


# tmply solution strategy is to keep it ENSG*
trs.rows[which(is.na(trs.rows$Gene.name)), "Gene.name"] <- trs.rows[which(is.na(trs.rows$Gene.name)), "ensemblID"]
#
trs.rows <- trs.rows[, c("ensemblID", "Gene.name")]
trs.rows[, "Types"] <- "Gene Expression"

# save to file
write.table(trs.rows, file = "../E-EHCA-2-quantification-raw-files/features.tsv", 
	sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

###### prepare final data
## read matrix
raw.counts <- Matrix::readMM("../E-EHCA-2-quantification-raw-files/matrix.sparse.mtx")
raw.features <- read.table("../E-EHCA-2-quantification-raw-files/features.tsv", sep = "\t", stringsAsFactors = FALSE)
#
raw.barcodes <- read.table("../E-EHCA-2-quantification-raw-files/barcodes.tsv", sep = "\t", stringsAsFactors = FALSE)
raw.barcodes <- as.character(raw.barcodes[, 1])
# set the dimnames
rownames(raw.counts) <- raw.features[, 2]
colnames(raw.counts) <- raw.barcodes
# remove those genes unmapped, which start with ENSMUSG
ind.to.rm.rows <- grep("^ENSMUSG", rownames(raw.counts))
raw.counts <- raw.counts[setdiff(1:nrow(raw.counts), ind.to.rm.rows), ]

## read meta data
raw.meta.data <- read.table("../E-EHCA-2-experiment-metadata-files/E-EHCA-2.sdrf.txt", 
	header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
# reserved cols
raw.meta.col.reserve <- c("Source.Name", "Characteristics.individual.", "Characteristics.sex.", 
	"Characteristics.age.", "Unit.time.unit.", 
	"Characteristics..organism.part.", "Characteristics..sampling.site.", 
	"Comment.submitted.inferred.cell.type.", "Characteristics.inferred.cell.type...authors.labels.", 
	"Characteristics.immunophenotype.", "Description")  # Description is about _ctrl or _day5/8/11
	
# replace colnames
raw.meta.col.replace <- c("barcodes", "it.ID", "it.sex", "it.age", "it.age.unit", 
	"sample.organism", "sample.site",  # <organism> is skin and lymph node, but <site> adds one kind more: tumor.
	"cluster.short", "cluster.name", 
	"cluster.immunophenotype", "sample.identity")  # sample.identity includes the information about if it is ctrl or the experiment part
#
raw.meta.data <- raw.meta.data[, raw.meta.col.reserve]
colnames(raw.meta.data) <- raw.meta.col.replace

raw.meta.data <- InterCellDB::DoPartUnique(raw.meta.data, which(colnames(raw.meta.data) %in% c("barcodes", "it.ID", "cluster.immunophenotype")))
# 6638 rows at 2020.12.23

tmp.barcodes.splits <- strsplit(raw.meta.data[, "barcodes"], split = "#", fixed = TRUE)
for(i in 1:length(tmp.barcodes.splits)) {  # check if it is aligned
	if (length(tmp.barcodes.splits[[i]]) != 2) {
		stop("Error!")
	}
}
tmp.barcodes.splits <- as.character(unlist(tmp.barcodes.splits))
tmp.barcodes.merge <- paste(tmp.barcodes.splits[1:(length(tmp.barcodes.splits) / 2) * 2 - 1], tmp.barcodes.splits[1:(length(tmp.barcodes.splits) / 2) * 2], sep = "_")
if ((length(which(colnames(raw.counts) %in% tmp.barcodes.merge)) == ncol(raw.counts)) != TRUE) {
	stop("NOT all IDs are included in META.DATA")  # further check
}
raw.meta.data[, "barcodes"] <- tmp.barcodes.merge
## find some cluster.name is NA
to.rm.barcodes <- raw.meta.data[which(is.na(raw.meta.data[, "cluster.short"])), "barcodes"]

#
 # raw counts changed here!!!!!! once more!
#
raw.counts <- raw.counts[, setdiff(colnames(raw.counts), to.rm.barcodes)]
raw.meta.data <- raw.meta.data[which(!is.na(raw.meta.data[, "cluster.short"])), ]
# 4626 cells, lost 1 [Endo LN] cell, a slight difference from the official website
# <official website> https://melanoma.cellgeni.sanger.ac.uk


### Seurat process 
tmp.obj <- CreateSeuratObject(counts = raw.counts, project = "mouse-xf1", min.cells = 3, min.features = 200)
#
tmp.obj[["percent.mt"]] <- PercentageFeatureSet(tmp.obj, pattern = "^mt-")
tmp.obj <- NormalizeData(tmp.obj, normalization.method = "LogNormalize", scale.factor = 10000)
tmp.obj <- FindVariableFeatures(tmp.obj, selection.method = "vst", nfeatures = 3000)
tmp.obj <- ScaleData(tmp.obj)
if (length(which(raw.meta.data$barcodes != rownames(tmp.obj@meta.data))) == 0) {
	tmp.obj@meta.data$seurat_clusters <- raw.meta.data$cluster.short
}
Idents(tmp.obj) <- tmp.obj@meta.data$seurat_clusters
#
#
#
# PCA top 24, resolution 0.4 to get 17 clusters
#
#
# run parallel FindAllMarkers
library(future)

plan("multiprocess", workers = 4)
tmp.markers <- FindAllMarkers(tmp.obj, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)



