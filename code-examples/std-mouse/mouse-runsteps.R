
# This file is the workflow for analyzing scRNA-seq data on mouse melanoma
# The steps for fetching marker genes are omitted.

# load library
library(InterCellDB)

# read differentially expressed genes
tmp.markers <- read.csv("./DEGs-mouse-clusters.csv", stringsAsFactors = FALSE)
# select the subset of passing bonferroni correction
fgenes.all <- tmp.markers[which(tmp.markers$p_val_adj < 0.05), ]  # 66332 rows left
# format the input data, to match the input requirement of InterCellDB
#  There are 4 required columns named 'LogFC', 'PVal', 'gene', 'cluster'.
colnames(fgenes.all)[3] <- "LogFC"
colnames(fgenes.all)[6] <- "PVal"

# [TODO] HERE to add [AddExprs to InterCell  structure]


# process
tmp.obj <- CreateInterCellObject(fgenes.all, species = "mouse")
# [TODO] HERE, inside add the expression data to fgenes

# prep for permutation of expression
# [TODO] This is the stand-alone step, will not influence the main process

# settle cluster names
ListAllClusters(tmp.obj)
tmp.obj <- ReplaceClusterName(tmp.obj, "Macrophage/Monocyte", "Myeloid")

# determine the subset of database to be used
tmp.obj <- SelectDBSubset(tmp.obj, sel.action.mode = "binding")

# fetch genes of interest
genes.receiver <- FetchGeneOI(tmp.obj, 
	sel.location = "Plasma Membrane",
	sel.location.score = c(4:5),
	sel.type = "Receptor",
	sel.go.terms = "GO:0006955"
	)
genes.sender <- FetchGeneOI(tmp.obj, 
	sel.location = "Extracellular Region",
	sel.location.score = c(4:5),
	sel.go.terms = "GO:0006955"
	)

# [TODO] add permutation here
tmp.permlist <- Tool.GenPermutation(seurat.obj, cells.meta, genes.name, perm.times = 1000)


# do network analysis
# [TODO] add switch here, whether to use the permutation test to get valid interactions
# [TODO] or just use the differentially expressed genes to quickly explore the interactions
tmp.obj <- AnalyzeInterInFullView(tmp.obj, 
	sel.some.genes.X = genes.receiver,
	sel.some.genes.Y = genes.sender,
	sel.exprs.change = c("Xup.Yup"))
# the result
tmp.res.fullview <- GetResultFullView(tmp.obj, 
	power.max.limit = 100, 
	hide.power.label = TRUE,
	plot.axis.x.name = "Cells involved in signal reception",
	plot.axis.y.name = "Cells involved in signal generation"
	)
Tool.ShowGraph(tmp.res.fullview)
Tool.WriteTables(tmp.res.fullview)

# explore deeper for intercellular analysis (between specific 2 cell clusters)
# usage of action mode and corresponding colors (mode 'other' is not concerned and 'binding' is pre-filtered by database selction)
used.action.mode <- c("activation", "inhibition", "catalysis", "reaction", "expression", "ptmod")
used.color.mode <- c("#FB8072", "#80B1D3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDB462")

# check Myeloid ~ CAF 1
tmp.obj <- FetchInterOI(tmp.obj, "Myeloid", "CAF 1")
tmp.obj <- AnalyzeInterInAction(tmp.obj)
Tool.ShowGraph(GetResultPieActionMode(tmp.obj, limits.exprs.change = c("Xup.Yup"), 
	limits.action.mode = used.action.mode,
	color.action.mode = used.color.mode))

# check Myeloid ~ CAF 2
tmp.obj <- FetchInterOI(tmp.obj, "Myeloid", "CAF 2")
tmp.obj <- AnalyzeInterInAction(tmp.obj)
Tool.ShowGraph(GetResultPieActionMode(tmp.obj, limits.exprs.change = c("Xup.Yup"), 
	limits.action.mode = used.action.mode,
	color.action.mode = used.color.mode))

# check Myeloid ~ CAF 3
tmp.obj <- FetchInterOI(tmp.obj, "Myeloid", "CAF 3")
tmp.obj <- AnalyzeInterInAction(tmp.obj)
Tool.ShowGraph(GetResultPieActionMode(tmp.obj, limits.exprs.change = c("Xup.Yup"), 
	limits.action.mode = used.action.mode,
	color.action.mode = used.color.mode))


# decide to use Myeloid ~ CAF 1 to further explore the intercellular communications
tmp.obj <- FetchInterOI(tmp.obj, "Myeloid", "CAF 1")

# select the postive and activation part of intercellular communications
tmp.obj <- SelectInterSubset(tmp.obj, 
	sel.action.mode = "activation",
	sel.action.effect = c("positive"),
	sel.action.merge.option = "union")

# filter and evaluate important gene pairs 
GetResultTgCrosstalk(tmp.obj, 
	direction.X.to.Y = NULL,
	plot.X.to.Y = TRUE,
	colnames.to.cmp = c("LogFC", "PVal"),
	axis.order.xy = c("LogFC", "LogFC"), 
	axis.order.xy.decreasing = c(TRUE, TRUE), 
	range.to.use = list("LogFC" = c(-Inf, Inf), "PVal" = c(-Inf, 40000)), 
	nodes.size.range = c(1, 6))$plot


# collect gene pairs that have higher power and higher confidence
tmp.obj <- SelectInterSubset(tmp.obj, 
	sel.some.genes.X = c("C5ar1","Ccr1","C3ar1","Itgam","Ccr2","Ccr5"),
	sel.some.genes.Y = c("C3","Ccl11","C4b","Cxcl1","Ccl7","Cxcl12","Ccl2")
	)

# evaluate the specificity of gene pairs (by explore the co-occurence of specific gene pairs among all interactions)
tmp.target.cluster.groups <- ListClusterGroups(tmp.obj, use.former = TRUE, cluster.former = c("Myeloid"))
tmp.obj <- AnalyzeInterSpecificity(tmp.obj, to.cmp.cluster.groups = tmp.target.cluster.groups)

# show the result of exploration on specificity
result.specificity <- GetResultTgSpecificity(tmp.obj,
	sel.uq.cnt.options = 1:10,  # this should be within the 1:(length(tmp.target.cluster.groups) + 1)
	sel.gene.pairs = NULL,
	plot.uq.cnt.merged = TRUE, 
	plot.name.X.to.Y = FALSE,
	dot.size = c(2, 8),
	dot.range.to.use = list(LogFC = c(-Inf, +Inf), PVal = c(-Inf, 40000)),
	dot.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
	dot.colour.value.seq = c(0.0, 0.5, 1.0)
	)
Tool.ShowGraph(result.specificity)

# Show the result of spatial pattern of selected gene pairs (Not used)
tmp.hide.locations <- setdiff(ListAllGeneLocation(tmp.obj), c("Plasma Membrane", "Extracellular Region"))
result.sptialpattern <- GetResultTgCellPlot(tmp.obj, 
	area.extend.times = 10,
	hide.other.area = TRUE,
	hide.locations.X = tmp.hide.locations,
	hide.locations.Y = tmp.hide.locations,
	link.size = 0.3,
	link.alpha = 0.8,
	legend.manual.left.spacing = grid::unit(0.1, "cm"))
Tool.ShowGraph(result.sptialpattern)

# draw ideal graph (only for result.sptialpattern)
pdf('./spatial-1.pdf', height = 12, width = 16)
Tool.ShowGraph(result.sptialpattern)
dev.off()


