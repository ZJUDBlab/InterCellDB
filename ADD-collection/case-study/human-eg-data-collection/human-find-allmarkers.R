# human Hepatic Seurat params
# Use Seurat v3.2

nice.colour <- c("#33a3dc",  # blue, slight light
	"#ef5b9c",  # pink colour
	"#1d953f",  # green, more thick
	"#00a6ac",  # blue, slight greenish
	"#a67c52",  # brown, light
	"#d71345",  # red, bit dark
	"#6a6da9",  # grey purple
	"#c7b299",  # white brown
	"#f3715c",  # slight orange
	"#fdb933",  # thick yellow
	"#ba9bc9",  # light purple
	"#f15a22",  # thick orange
	"#5c7a29"   # deep green
	)  

# readin
tmp.relative.path <- "../CountData/"
tmp.filenames <- list.files(tmp.relative.path)
tmp.filenames <- paste0(tmp.relative.path, tmp.filenames)
# avoid importing these data once together
# file size is listed as below
#   The 1 extra column records the gene name not the cell name
# @1: "../CountData/GSM4116579_ICC_18_Adjacent_UMI.csv", ncol: 10318(1 extra), nrow: 14647
# @2: "../CountData/GSM4116580_ICC_18_Tumor_UMI.csv", ncol: 4964(1 extra), nrow: 16423
# @3: "../CountData/GSM4116581_ICC_20_Tumor_UMI.csv", ncol: 2715(1 extra), nrow: 14670
# @4: "../CountData/GSM4116582_ICC_23_Adjacent_UMI.csv", ncol: 4426(1 extra), nrow: 15308
# @5: "../CountData/GSM4116583_ICC_23_Tumor_UMI.csv", ncol: 3505(1 extra), nrow: 16566
# @6: "../CountData/GSM4116584_ICC_24_Tumor1_UMI.csv", ncol: 3469(1 extra), nrow: 17353
# @7: "../CountData/GSM4116585_ICC_24_Tumor2_UMI.csv", ncol: 2442(1 extra), nrow: 16570
# @8: "../CountData/GSM4116586_ICC_25_Adjacent_UMI.csv", ncol: 2160(1 extra), nrow: 13884
# as it is big data
tmp.cnt.1 <- read.csv(tmp.filenames[1], header = TRUE)
# .....
merge.cnt <- full_join(tmp.cnt.1, tmp.cnt.2, by = c("X" = "X"))
merge.cnt <- full_join(merge.cnt, tmp.cnt.3, by = c("X" = "X"))
merge.cnt <- full_join(merge.cnt, tmp.cnt.4, by = c("X" = "X"))
#
merge.cnt <- full_join(merge.cnt, tmp.cnt.5, by = c("X" = "X"))
merge.cnt <- full_join(merge.cnt, tmp.cnt.6, by = c("X" = "X"))
merge.cnt <- full_join(merge.cnt, tmp.cnt.7, by = c("X" = "X"))
merge.cnt <- full_join(merge.cnt, tmp.cnt.8, by = c("X" = "X"))
# change to matrix
rownames(merge.cnt) <- merge.cnt[, 1]
merge.cnt <- merge.cnt[, -1]
merge.cnt <- as.matrix(merge.cnt)
# count NAs, and remove those genes [deprecated]
{
	prog.bar.p.p <- progress::progress_bar$new(total = nrow(merge.cnt))
	prog.bar.p.p$tick(0)
	tmp.hasNA <- logical()
	for(i in 1:nrow(merge.cnt)) {
		prog.bar.p.p$tick()
		tmp.hasNA <- c(tmp.hasNA, !is.na(match(NA, merge.cnt[i, ])))
	}
	merge.cnt <- merge.cnt[which(tmp.hasNA == FALSE), ]  # finally, get 11737 rows
}
#
# use all merged cnt, some may have NAs but it doesn't matter, and will be treated as 0
merge.cnt <- readRDS("../backup-COUNT-ALL.rds")
# [TO NOTE] the following several steps for counts always needed re-run when reading
prog.bar.p.p <- progress::progress_bar$new(total = ncol(merge.cnt))
prog.bar.p.p$tick(0)
for(i in 1:ncol(merge.cnt)) {
	prog.bar.p.p$tick()
	tmp.onecol <- merge.cnt[, i]
	tmp.onecol[which(is.na(tmp.onecol))] <- 0
	tmp.onecol <- as.numeric(tmp.onecol)
	merge.cnt[, i] <- tmp.onecol
}
merge.cnt <- Matrix(merge.cnt, sparse = TRUE)
# Seurat ver2
tmp.obj <- CreateSeuratObject(raw.data = merge.cnt, project = "hepatic-human-z1", min.cells = 3, min.genes = 200)

# calculate percent of mt gene
# e.g. kpattern.mtgene <- "^mt-"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = tmp.obj@data), value = TRUE)
percent.mito <- Matrix::colSums(tmp.obj@raw.data[mito.genes, ])/Matrix::colSums(tmp.obj@raw.data)
# QC
tmp.obj <- AddMetaData(object = tmp.obj, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = tmp.obj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# May find help from FeatureScatter() in Seurat pkg
tmp.obj <- FilterCells(object = tmp.obj, subset.names = c("nGene", "nUMI", "percent.mito"), 
    low.thresholds = c(500, 2000, -Inf), high.thresholds = c(6000, +Inf, 0.20))
# 30728 cells, 19813 genes
tmp.obj <- NormalizeData(object = tmp.obj, normalization.method = "LogNormalize", 
    scale.factor = 10000)
# find highly variable features
tmp.obj <- FindVariableGenes(object = tmp.obj, 
	mean.function = ExpMean, dispersion.function = LogVMR,  # the default settings
    x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 0.5) 
# scale the data [NOTE] need huge memory
tmp.obj <- ScaleData(object = tmp.obj, 
	genes.use = tmp.obj@var.genes,
	vars.to.regress = c("nUMI", "percent.mito"))
tmp.obj <- RunPCA(object = tmp.obj, pc.genes = tmp.obj@var.genes, 
	do.print = TRUE, pcs.print = 1:5,  genes.print = 5)
#
kdims.pca <- 20   # [TO NOTE!!!]
kresolution <- 0.3
tmp.obj <- FindClusters(object = tmp.obj, reduction.type = "pca", dims.use = 1:kdims.pca, 
    resolution = kresolution, print.output = 0, save.SNN = TRUE)
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 5336 4733 3923 3075 2825 2394 2308 1251 1110  739  713  597  513  454  453  304 
# run non-linear dimensional reduction
tmp.obj <- RunTSNE(object = tmp.obj, dims.use = 1:kdims.pca, do.fast = TRUE)
#              B_cell        Cholangiocyte          Endothelial 
#                1167                  597                 1110 
#          Fibroblast           Hepatocyte Macrophage&Dendritic 
#                 513                  453                 4966 
#           Malignant            T_cell&NK 
#               12260                 9662 
#
# TO further clustering fibroblast
# get merge.cnt
merge.cnt
# get fibroblast cells
celltype.m.list <- read.csv("../cell-type-matching-list.csv", header = TRUE, stringsAsFactors = FALSE)
celltype.m.list <- celltype.m.list[, c(2:3)]
celltype.fibroblast.names <- celltype.m.list[which(celltype.m.list$type == "Fibroblast"), "cell"]
# subset merge.cnt
fibro.cnt <- merge.cnt[, colnames(merge.cnt) %in% celltype.fibroblast.names]
# following use Seurat ver3
tmp.fibro <- CreateSeuratObject(counts = fibro.cnt, project = "fibroblast-recalc", 
			min.cells = 3, min.features = 200)
tmp.fibro <- NormalizeData(tmp.fibro, normalization.method = "LogNormalize", scale.factor = 10000)
# find highly variable features
tmp.fibro <- FindVariableFeatures(tmp.fibro, selection.method = "vst", nfeatures = 3000)
tmp.fibro <- ScaleData(tmp.fibro)
# perform linear dimensional reduction
tmp.fibro <- RunPCA(tmp.fibro, npcs = 50, features = VariableFeatures(object = tmp.fibro))
# determine the dimensionality of the data

tmp.fibro <- JackStraw(tmp.fibro, num.replicate = 50, dims = 50)
tmp.fibro <- ScoreJackStraw(tmp.fibro, dims = 1:50)  # the max val here decided by the JackStraw @param dims
JackStrawPlot(tmp.fibro, dims = 1:50)
#
kdims.pca.fibro <- 8
kresolution.fibro <- 0.3
# cluster the cells
tmp.fibro <- FindNeighbors(tmp.fibro, dims = 1:kdims.pca.fibro)
# e.g. kresolution <- 0.4 ~ 1.2  (recommended)  # self
tmp.fibro <- FindClusters(tmp.fibro, resolution = kresolution.fibro)

# run non-linear dimensional reduction
tmp.fibro <- RunTSNE(tmp.fibro, dims = 1:kdims.pca.fibro)
DimPlot(tmp.fibro, reduction = "tsne")
# find markers
fibro.markers <- FindAllMarkers(tmp.fibro, only.pos=TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# check clusters
top10 <- fibro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(tmp.fibro, features = top10$gene) + NoLegend()
# saveRDS(tmp.fibro@meta.data, "../meta-fibro-fromall-subset.rds")
# saveRDS(tmp.fibro.clusters, "../fibro-fromall-subset-clusters.rds")
#
# so all sub-clustering is done
# reload merge.cnt
merge.cnt
# align the new celltype
celltype.m.list <- read.csv("../cell-type-matching-list.csv", header = TRUE, stringsAsFactors = FALSE)
celltype.m.list <- celltype.m.list[, c(2:3)]
fibro.sub.list <- readRDS("../fibro-fromall-subset-clusters.rds")
fibro.meta <- readRDS("../meta-fibro-fromall-subset.rds")
fibro.meta$seurat_clusters <- factor(fibro.meta$seurat_clusters)
levels(fibro.meta$seurat_clusters) <- fibro.sub.list
celltype.fibro.list <- data.frame(cell = rownames(fibro.meta), type = fibro.meta$seurat_clusters, stringsAsFactors = FALSE)
#
celltype.m.list <- celltype.m.list[-which(celltype.m.list$cell %in% celltype.fibro.list$cell), ]
#
celltype.final <- rbind(celltype.m.list, celltype.fibro.list)
#
tmp.new <- CreateSeuratObject(counts = merge.cnt, project = "opz1-final", min.cells = 3, min.features = 200)
tmp.new[["percent.mt"]] <- PercentageFeatureSet(tmp.new, pattern = "^MT-")
tmp.new <- subset(tmp.new, subset = 
			nFeature_RNA > 500
			& nFeature_RNA < 6000
			& percent.mt < 20
			& nCount_RNA > 2000)
# 30732 cells
tmp.new <- NormalizeData(tmp.new, normalization.method = "LogNormalize", scale.factor = 10000)
# find highly variable features
tmp.new <- FindVariableFeatures(tmp.new, selection.method = "vst", nfeatures = 3000)
tmp.new <- ScaleData(tmp.new)
# dummy run TSNE
if (FALSE) {
	tmp.human.obj <- RunPCA(tmp.human.obj, npcs = 50,  features = VariableFeatures(object = tmp.human.obj))
	tmp.human.obj <- FindNeighbors(tmp.human.obj, dims = 1:22)
	tmp.human.obj <- FindClusters(tmp.human.obj, resolution=0.8)
	tmp.human.obj <- RunTSNE(tmp.human.obj, dims=1:22)
	DimPlot(tmp.human.obj, reduction = "tsne")
}
# 
celltype.final <- celltype.final[match(names(Idents(tmp.new)), celltype.final[, 1]), ]
Idents(tmp.new) <- celltype.final[, 2]
# new markers
tmp.final.markers <- FindAllMarkers(tmp.new, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)
# A good result with IL6 as DEG in vCAFs, and good result for downstream analysis




## added at 2021.03.18
# further splits the definition of cluster T cell and NK, & Macrophage & Dendritic
# readin
tmp.obj  # Seurat obj
celltype.m.df <- read.csv("../cell-type-matching-list.csv", header = TRUE, stringsAsFactors = FALSE)
celltype.m.df <- celltype.m.df[, c(2,3)]
#
## splits T cell and NK cell
	cells.tnk <- celltype.m.df[which(celltype.m.df[, "type"] == "T_cell&NK"), "cell"]
	tmp.tnk.obj <- subset(tmp.obj, cells = cells.tnk)
	tmp.tnk.obj <- RunPCA(tmp.tnk.obj, npcs = 50, features = VariableFeatures(object = tmp.tnk.obj))
	# determine the dimensionality of the data
	ElbowPlot(tmp.tnk.obj, ndims = 30)
	# cluster the cells
	tmp.tnk.obj <- FindNeighbors(tmp.tnk.obj, dims = 1:19)
	# e.g. kresolution <- 0.4 ~ 1.2  (recommended)  # self
	tmp.tnk.obj <- FindClusters(tmp.tnk.obj, resolution = 0.3)

	# run non-linear dimensional reduction
	tmp.tnk.obj <- RunTSNE(tmp.tnk.obj, dims = 1:19)
	#
	DimPlot(tmp.tnk.obj, reduction = "tsne", label = TRUE)
	#
	FeaturePlot(tmp.tnk.obj, features = c("KLRF1", "CD8A", "CD7", "FGFBP2", "CD3D", "GZMK", "GZMB"), label = TRUE)
	VlnPlot(tmp.tnk.obj, features = c("KLRF1", "CD8A", "CD7", "FGFBP2", "CD3D", "GZMK", "GZMB"), ncol = 3)
	VlnPlot(tmp.tnk.obj, features = c("FGFBP2", "GZMB", "FCGR3A", "GZMH"), ncol = 2)
	VlnPlot(tmp.tnk.obj, features = c("TYROBP", "FCER1G", "GZMK", "CD160"), ncol = 2)
	# get reclustered
	toc.nk.cells <- rownames(tmp.tnk.obj@meta.data[which(tmp.tnk.obj@meta.data$seurat_clusters %in% c(5,7,8)), ])
	toc.t.cells <- setdiff(rownames(tmp.tnk.obj@meta.data), toc.nk.cells)
	#
#
# splits Macrophage and Dendritic celll
	cells.mpd <- celltype.m.df[which(celltype.m.df[, "type"] == "Macrophage&Dendritic"), "cell"]
	tmp.mpd.obj <- subset(tmp.obj, cells = cells.mpd)
	tmp.mpd.obj <- RunPCA(tmp.mpd.obj, npcs = 50, features = VariableFeatures(object = tmp.mpd.obj))
	# determine the dimensionality of the data
	ElbowPlot(tmp.mpd.obj, ndims = 30)
	# cluster the cells
	tmp.mpd.obj <- FindNeighbors(tmp.mpd.obj, dims = 1:19)
	# e.g. kresolution <- 0.4 ~ 1.2  (recommended)  # self
	tmp.mpd.obj <- FindClusters(tmp.mpd.obj, resolution = 0.3)

	# run non-linear dimensional reduction
	tmp.mpd.obj <- RunTSNE(tmp.mpd.obj, dims = 1:19)
	#
	DimPlot(tmp.mpd.obj, reduction = "tsne", label = TRUE)
	#
	FeaturePlot(tmp.mpd.obj, features = c("CD14", "CD1C", "CLEC9A", "IDO1"), label = TRUE)
	VlnPlot(tmp.mpd.obj, features = c("CD14", "CD1C", "CLEC9A", "IDO1"), ncol = 2)
	# get reclustered
	toc.dendritic.cells <- rownames(tmp.mpd.obj@meta.data[which(tmp.mpd.obj@meta.data$seurat_clusters %in% c(3,6,7,8)), ])
	toc.macrophage.cells <- setdiff(rownames(tmp.mpd.obj@meta.data), toc.dendritic.cells)

#
# replace all new clustr name
	celltype.m.new.df <- celltype.m.df
	tmp.new.type <- celltype.m.new.df[, "type"]
	tmp.new.type[which(celltype.m.new.df$cell %in% toc.nk.cells)] <- "NK"
	tmp.new.type[which(celltype.m.new.df$cell %in% toc.t.cells)] <- "T_cell"
	tmp.new.type[which(celltype.m.new.df$cell %in% toc.dendritic.cells)] <- "Dendritic"
	tmp.new.type[which(celltype.m.new.df$cell %in% toc.macrophage.cells)] <- "Macrophage"
	celltype.m.new.df[, "type"] <- tmp.new.type
	# save new type list
	write.csv(celltype.m.new.df, "../new-celltype-allsplited.csv")

# rerun the process
	tmp.obj <- RunPCA(tmp.obj, npcs = 50, features = VariableFeatures(object = tmp.obj))
	# determine the dimensionality of the data
	ElbowPlot(tmp.obj, ndims = 50)
	# cluster the cells
	tmp.obj <- FindNeighbors(tmp.obj, dims = 1:22)
	# e.g. kresolution <- 0.4 ~ 1.2  (recommended)  # self
	tmp.obj <- FindClusters(tmp.obj, resolution = 0.3)

	# run non-linear dimensional reduction
	tmp.obj <- RunTSNE(tmp.obj, dims = 1:22)
	# readin clusters
	celltype.m.new.df <- read.csv("../new-celltype-allsplited.csv", header = TRUE, stringsAsFactors = FALSE)
	celltype.m.new.df <- celltype.m.new.df[, 2:3]
	#
	check.to.replace <- sum(names(Idents(tmp.obj)) == celltype.m.new.df[, "cell"]) == nrow(celltype.m.new.df)
	if (check.to.replace == FALSE) {
		stop("Not aligned data!")
	} else {
		Idents(tmp.obj) <- celltype.m.new.df[, "type"]
	}
	#
	DimPlot(tmp.obj, reduction = "tsne", label = TRUE)
library(future)
plan("multiprocess", workers = 2)
all.markers <- FindAllMarkers(tmp.obj, only.pos=FALSE, min.pct = 0.1, logfc.threshold = 0)



## merge in Seurat ver3 
#tmp.clusters <- readRDS("../clusters-1st-16-groups.rds")
#tmp.clusters <- tmp.clusters[order(tmp.clusters)]
#levels(tmp.obj@meta.data$seurat_clusters) <- names(tmp.clusters)
#Idents(tmp.obj) <- tmp.obj@meta.data$seurat_clusters
## merge in Seurat ver2
#tmp.clusters <- readRDS("../new-cluster-Seuratv2.rds")  # has ordered
#tmp.obj@meta.data$seurat_clusters <- tmp.obj@meta.data$res.0.3  # this colname will change!!!
#tmp.obj@meta.data$seurat_clusters <- as.integer(tmp.obj@meta.data$seurat_clusters)
#tmp.obj@meta.data$seurat_clusters <- factor(tmp.obj@meta.data$seurat_clusters)
#for(i in 1:16) {
	#tmp.obj <- RenameIdent(tmp.obj, old.ident.name= i-1, new.ident.name= tmp.clusters[i])
#}


# added at 2021.03.13 for drawing nice TSNE gragh
tmp.article.cols <- c("#AEC7E8",  # T_cell&NK
"#BCBD21",  # Malignant
"#1F77B4",  # B_cell
"#8C564C",  # Hepatocyte
"#D62727",  # Endothelial
"#9467BD",  # Fibroblast
"#E377C3",  # Macrophage&Dendritic
"#FF7F0F"   # Cholangiocyte
)
names(tmp.article.cols) <- c("T_cell&NK", "Malignant", "B_cell", 
	"Hepatocyte", "Endothelial", "Fibroblast", "Macrophage&Dendritic", "Cholangiocyte")




