### CellTalkDB step 1
## download example data from website
# option.1: Tutorial in https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# gives one example data. The data links are as below:
# https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# what you get after downloading?
# ONE file named pbmc3k_filtered_gene_bc_matrices.tar.gz
# By decompress it, you will find 3 files in path "./filtered_gene_bc_matrices/hg19/".
# (1)barcodes.tsv; (2)genes.tsv; (3)matrix.mtx;
# Those are all data we need in this example, and move the directory under the "code-examples".
# So, those will be in "CellTalkDB/code-examples/filtered_gene_bc_matrices/hg19/".


### CellTalkDB step 2
## process the raw scRNA-seq data with Seurat
# The codes below are the same as codes provided in tutorial in Seurat website,
# (https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)
# but given as a easy-to-use referece codes.

# start R and run all codes below!
library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc  # This is the Seurat object.

# add MT- DNAs into consideration
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# estimate the cutoff for selection of one data subset
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# select part of the data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize data
pbmc <- NormalizeData(pbmc)

# find feature genes 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# scale the data
pbmc <- ScaleData(pbmc)

# run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# use this to estimate the dimensionality
ElbowPlot(pbmc)

# find dimensionality in more accurate way
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
# !attention!
# A sharp drop-off in significance after the first 10-12 PCs in this data
# For other dataset, this step needs to judge by yourself!
this.choose.dims <- 10

# cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:this.choose.dims)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# run tSNE
pbmc <- RunTSNE(pbmc, dims = 1:this.choose.dims)
DimPlot(pbmc, reduction = "tsne")

# FindAllMarkers: Finding differentially expressed features
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)

# when reaching here, the process is done
pbmc.markers  # it is the result CellTalkDB most need!
pbmc  # it is the Seurat object, which contains all results as processed above.

# choose to save it
saveRDS(pbmc.markers, "./fgenes-example.rds")  #!!!must saved
saveRDS(pbmc, "./seurat-out-example.rds")  # (optional) as sometimes taken too large disk space



