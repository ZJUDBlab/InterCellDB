
# This file is the workflow for analyzing scRNA-seq data on human cholangiocarcinoma
# The steps for fetching marker genes are omitted.

# load library
library(InterCellDB)

# read differentially expressed genes
tmp.markers <- read.csv("./DEGs-human-clusters.csv", stringsAsFactors = FALSE)
# select the subset of passing bonferroni correction
fgenes.all <- tmp.markers[which(tmp.markers$p_val_adj < 0.05), ]  # 53147 rows left
# format the input data, to match the input requirement of InterCellDB
#  There are 4 required columns named 'LogFC', 'PVal', 'gene', 'cluster'.
colnames(fgenes.all)[3] <- "LogFC"
colnames(fgenes.all)[6] <- "PVal"

# process
tmp.obj <- CreateInterCellObject(fgenes.all, species = "human")


## step 1: show the gene pairs selected by the original article
# selected gene pairs
tmp.s1.gpairs <- data.frame(gene.a = c("MIF", "GRN", "C3", "CD99", "NRP2", "GRN", "TNFSF10", "EREG", "HBEGF", "AREG", "MIF", "JAG1", "TNFSF10", "CLEC2B", "BST2", "IL6R", "JAG1", "JAG2", "PDGFD", "NRP2"), 
	gene.b = c("TNFRSF14", "TNFRSF1B", "C3AR1", "PILRA", "VEGFA", "TNFRSF1A", "TNFRSF11B", "EGFR", "EGFR", "EGFR", "TNFRSF10D", "NOTCH4", "TNFRSF10D", "KLRF1", "LILRA4", "IL6", "NOTCH3", "NOTCH3", "PDGFRB", "PGF"),
	stringsAsFactors = FALSE)
# as the gene pairs selected is indirectional in original article
tmp.sel.gene.pairs.set <- FormatCustomGenePairs(tmp.s1.gpairs, 
	species = "human", 
	extend.reverse = FALSE)
# get the formatted gene pairs
tmp.sel.gene.pairs <- tmp.sel.gene.pairs.set$result

# do network analysis
tmp.obj <- AnalyzeInterInFullView(tmp.obj, 
	sel.gene.pairs = tmp.sel.gene.pairs)
# fetch the target cluster group
tmp.obj <- FetchInterOI(tmp.obj, "Malignant", "Fibroblast")
# evaluate gene pair specificity
tmp.to.cmp.cluster.groups <- ListClusterGroups(tmp.obj, use.former = TRUE, cluster.former = c("Malignant"))
tmp.obj <- AnalyzeInterSpecificity(tmp.obj, 
	to.cmp.cluster.groups = tmp.to.cmp.cluster.groups,
	extended.search = TRUE)
result.specificity <- GetResultTgSpecificity(tmp.obj,
	sel.uq.cnt.options = 1:10,  # this should be within the 1:(length(tmp.target.cluster.groups) + 1)
	prioritize.cluster.groups = tmp.to.cmp.cluster.groups[c(1,3,6,4,2,5,7,8,9,10)],
	plot.uq.cnt.merged = TRUE, 
	plot.name.X.to.Y = TRUE,
	dot.size = c(2, 8),
	dot.range.to.use = list(LogFC = c(0, +Inf), PVal = c(-Inf, 120000)),
	dot.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
	dot.colour.value.seq = c(0.0, 0.5, 1.0),
	dot.y.order.in.alphabet = FALSE
	)


## step 2: show all interleukin genes involved gene pairs
# get interleukin-related genes (collect from gene reference database)
tmp.IL.genes <- read.csv("./interleukin-human-genes.csv", stringsAsFactors = FALSE)
tmp.IL.R.genes <- read.csv("./interleukin-receptor-human-genes.csv", stringsAsFactors = FALSE)

# renewing the object and select the physical associated gene pairs
tmp.obj <- CreateInterCellObject(fgenes.all, species = "human")
tmp.obj <- SelectDBSubset(tmp.obj, sel.action.mode = "binding")

tmp.obj <- AnalyzeInterInFullView(tmp.obj, 
	sel.some.genes.X = tmp.IL.R.genes$gene.name,
	sel.some.genes.Y = tmp.IL.genes$gene.name)

# get the target cluster group
tmp.obj <- FetchInterOI(tmp.obj, "Malignant", "Fibroblast")

# evaluate gene pair specificity
tmp.to.cmp.cluster.groups <- ListClusterGroups(tmp.obj, use.former = TRUE, cluster.former = c("Malignant"))
tmp.obj <- AnalyzeInterSpecificity(tmp.obj, 
	to.cmp.cluster.groups = tmp.to.cmp.cluster.groups,
	extended.search = TRUE)
result.specificity <- GetResultTgSpecificity(tmp.obj,
	sel.uq.cnt.options = 1:10,  # this should be within the 1:(length(tmp.target.cluster.groups) + 1)
	prioritize.cluster.groups = tmp.to.cmp.cluster.groups[c(1,3,6,4,2,5,7,8,9,10)],
	plot.uq.cnt.merged = TRUE, 
	plot.name.X.to.Y = FALSE,
	dot.size = c(2, 8),
	dot.range.to.use = list(LogFC = c(0, +Inf), PVal = c(-Inf, 120000)),
	dot.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
	dot.colour.value.seq = c(0.0, 0.5, 1.0),
	dot.y.order.in.alphabet = TRUE
	)


## step 3: show all IL6-* gene pairs
# renewing the object and select the physical associated gene pairs
tmp.obj <- CreateInterCellObject(fgenes.all, species = "human")
tmp.obj <- SelectDBSubset(tmp.obj, sel.action.mode = "binding")

# interacting partner in plasma membrane
genes.in.pmem <- FetchGeneOI(tmp.obj, 
	sel.location = "Plasma Membrane",
	sel.location.score = c(4:5)
	)
# network analysis
tmp.obj <- AnalyzeInterInFullView(tmp.obj, 
	sel.some.genes.X = genes.in.pmem,
	sel.some.genes.Y = c("IL6"))

# show the result
Tool.ShowGraph(
	GetResultFullView(tmp.obj, power.max.limit = 3, hide.power.label = TRUE)
)

# get the target cluster group
tmp.obj <- FetchInterOI(tmp.obj, "Malignant", "Fibroblast")

# get the Fibroblast-acting-on-Malignant pairs
tmp.obj <- SelectInterSubset(tmp.obj, direction.X.to.Y = FALSE)


# show the spatial pattern result
tmp.hide.locations <- setdiff(ListAllGeneLocation(tmp.obj), c("Plasma Membrane", "Extracellular Region"))
result.il6.sptialpattern <- GetResultTgCellPlot(tmp.obj, 
	area.extend.times = 10,
	hide.other.area = TRUE,
	hide.locations.X = c(tmp.hide.locations, "Extracellular Region"),
	hide.locations.Y = c(tmp.hide.locations, "Plasma Membrane"),
	link.colour = c("#D70051", "#00913A", "#1296D4", "#956134", "#C8DC32", "#B5B5B6", "#0A0AFF", "#762A83"),
	legend.manual.left.spacing = grid::unit(0.1, "cm"))
Tool.ShowGraph(result.il6.sptialpattern)
