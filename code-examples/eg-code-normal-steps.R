# This file gives the example of codes for step 3~15
# The data is got from the result from "./eg-code-Seurat.R"
library(CellTalkDB)

### get result from CellTalkDB step 2
## feature genes
fgenes.all <- readRDS("./fgenes-example.rds")
fgenes.all$cluster <- as.character(fgenes.all$cluster)  # fgenes.all$cluster in default is factor(), so change it to be character().
## Seurat object
this.Seurat.obj <- readRDS("./seurat-out-example.rds")


### CellTalkDB step 3 (optional)
## This step should be strictly run before step 4 as step 4 changes the gene names, 
## which will cause runtime error.
fgenes.all <- DataPrep.AddExprsVals(fgenes.all, this.Seurat.obj)

### CellTalkDB step 4
fgenes.all <- DataPrep.RemapClustersMarkers(fgenes.all, genes.human.ref.db)

### CellTalkDB step 5 (optional)
tmp.cluster.names.current <- c("0", "1", "2")
tmp.cluster.names.replace <- c("x0", "x1", "x2")
fgenes.all <- DataPrep.AddClusterName(fgenes.all, tmp.cluster.names.current, tmp.cluster.names.replace)

### CellTalkDB step 6
## set some parameters
tmp.sub.sel.exprs.changes <- c("Xup.Yup")
tmp.sub.sel.X.Location <- c("Plasma Membrane", "Endoplasmic Reticulum")
tmp.sub.sel.Y.Location <- c("Nucleus", "Extracellular Region")

### CellTalkDB step 7
## run the analysis in fullview, use human experiments database, and human annotation databases
this.fullview.analysis <- AnalyzeClustersInteracts(fgenes.all, 
	pairs.ref = pairs.human.experiments.db,
	anno.location.ref = anno.location.human.ref.db,
	anno.type.ref = anno.type.human.ref.db,
	restricted.some.genes = NULL,  # currently no restriction applied on genes
	restricted.gene.pairs = NULL,  # currently no restriction applied on gene pairs
	sub.sel.exprs.changes = tmp.sub.sel.exprs.changes,
	sub.sel.X.Location = tmp.sub.sel.X.Location,
	sub.sel.Y.Location = tmp.sub.sel.Y.Location,
	sub.sel.X.Type = NULL,  # currently no restriction applied on gene types in X
	sub.sel.Y.Type = NULL   # currently no restriction applied on gene types in Y
	)

### CellTalkDB step 8
## get the result of fullview analysis
this.res.fullview <- GetResult.SummaryClustersInteracts(this.fullview.analysis, power.max.limit=1.2, cnt.max.limit=150)
# show graph
this.res.fullview$plot  # shortcuts
Tool.ShowGraph(this.res.fullview)  # normal usage
# write result tables
Tool.WriteTables(this.res.fullview, dir.path = "./")


### CellTalkDB step 9
## select the 2nd most interesting cluster pairs
this.target.cp <- ExtractTargetOnepairClusters(this.fullview.analysis, "x2", "8")

### CellTalkDB step 10
## use human actions database
this.target.gmoc <- GenerateMapDetailOnepairClusters(this.target.cp, actions.human.ref.db)

### CellTalkDB step 11
this.res.targetview <- GetResult.SummaryOnepairClusters(this.target.gmoc, 
	show.exprs.change = c("Cup.Dup", "Cup.Ddn"),
	show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B")
	)
## as up-up gene pairs are pre-selected, so here only 1 column will be shown
library(ggplot2)  # as obviously use ggplot2 here
Tool.ShowGraph(this.res.targetview, scale_fill_brewer(palette = "Set3"))
## write result tables
Tool.WriteTables(this.res.targetview, dir.path = "./")


### CellTalkDB step 12
this.target.VEinfos <- GenerateVEinfos(this.target.gmoc, fgenes.all,
	is.directional = FALSE,
	sel.mode.val = NULL,  # select all kinds of modes
	sel.action.effect.val = c("positive", "negative")
	)

### CellTalkDB step 13
this.res.tgenes <- GetResult.PlotOnepairClusters.CellPlot(this.target.VEinfos,
	hide.locations.A = NULL,
	hide.locations.B = NULL,
	hide.sole.vertices = TRUE,
	nodes.label.size = 3,
	nodes.label.colour = "white"
	)
## show graph
Tool.ShowGraph(this.res.tgenes)
## write result tables
Tool.WriteTables(this.res.tgenes, dir.path = "./")


### CellTalkDB step 14
## get result from step 7: this.fullview.analysis
## choose 2 cluster pairs, 1: x2~8, 2: 8~8
tmp.select.interacts <- c(paste0("x2", kClustersSplit, "8"), paste0("8", kClustersSplit, "8"))

### CellTalkDB step 15
## the plot will be drawn, and the result tables will be saved by "<-"
this.tables.circos <- PlotClusters.Circos(this.fullview.analysis, fgenes.all, actions.human.ref.db,
	return.table.vals = TRUE,
	select.interacts = tmp.select.interacts,
	is.directional = TRUE,
	sel.mode.val = NULL,  # select all kinds of modes
	sel.action.effect.val = c("positive", "negative"),
	show.legend = TRUE
	)
# can write result tables in the similar way
Tool.WriteTables(this.tables.circos, dir.path = "./")

