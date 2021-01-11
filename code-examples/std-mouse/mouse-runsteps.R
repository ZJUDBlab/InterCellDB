
# This file is the workflow for analyzing scRNA-seq data on mouse melanoma
# The steps for fetching marker genes are omitted.

# load library
library(InterCellDB)

# read differentially expressed genes
tmp.markers <- read.csv("./DEGs-all-clusters.csv", stringsAsFactors = FALSE)
# select the subset of passing bonferroni correction
fgenes.all <- tmp.markers[which(tmp.markers$p_val_adj < 0.05), ]  # 66332 rows left

# remapping gene names to InterCellDB built-in standard gene information database
# some genes cannot be remapped, but by manuall check, those genes are mostly non-determining-name.
fgenes.all <- DataPrep.RemapClustersMarkers(fgenes.all, genes.mouse.ref.db)$result
# rename some cluster
fgenes.all <- DataPrep.AddClusterName(fgenes.all, "Macrophage/Monocyte", "Myeloid")

###### fullview analysis
# set the expression changes
tmp.sub.sel.exprs.changes <- c("Xup.Yup")

# one more step - select the high confidence interaction database
tmp.db.exp.only <- pairs.mouse.db[which(pairs.mouse.db$inter.Experiments.Score != 0), ]
tmp.db.rest.high.conf <- pairs.mouse.db[intersect(which(pairs.mouse.db$inter.Combined.Score >= 700),
	which(pairs.mouse.db$inter.Experiments.Score == 0)), ]
ref.data <- rbind(tmp.db.exp.only, tmp.db.rest.high.conf)

# fetch genes of GO:0006955 [immune response] from GO
ref.GO.immune.response <- Tool.FindGenesFromGO("GO:0006955", genes.mouse.ref.db, go.mouse.ref.db, 
	go.use.relative = TRUE, go.relative.option = "offspring")
tmp.go.im <- data.frame(gene.name = ref.GO.immune.response[[1]], user.type = "immune", stringsAsFactors = FALSE)
tmp.go.im.df <- Tool.AddUserRestrictDB(tmp.go.im, genes.mouse.ref.db)

# analyze the immune-response-related process in these cell groups
this.0.anal <- AnalyzeClustersInteracts(fgenes.all, 
	pairs.ref = ref.data, 
	anno.location.ref = anno.location.mouse.ref.db,
	anno.type.ref = anno.type.mouse.ref.db,
	sub.sel.X.clusters = NULL,
	sub.sel.Y.clusters = NULL,
	sub.sel.exprs.changes = tmp.sub.sel.exprs.changes,
	user.type.database = tmp.go.im.df,
	sub.sel.user.type.colname = "user.type",
	sub.sel.X.user.type = "immune",
	sub.sel.Y.user.type = "immune",
	sub.sel.user.merge.option = "union"
	)
# pics- id.1.left: indicate the immune response concentrating part
GetResult.SummaryClustersInteracts(this.0.anal, power.max.limit = 3000, hide.power.label = TRUE)$plot


# analyze further by adding interaction pattern to be paracrine
this.anal.immune <- AnalyzeClustersInteracts(fgenes.all, 
	pairs.ref = ref.data, 
	anno.location.ref = anno.location.mouse.ref.db,
	anno.type.ref = anno.type.mouse.ref.db,
	sub.sel.X.clusters = NULL,
	sub.sel.Y.clusters = NULL,
	sub.sel.exprs.changes = tmp.sub.sel.exprs.changes,
	sub.sel.X.Location = c("Extracellular Region"),
	sub.sel.X.Location.score.limit = c(3:5),
	sub.sel.Y.Location = c("Plasma Membrane"),
	sub.sel.Y.Location.score.limit = c(3:5),
	user.type.database = tmp.go.im.df,
	sub.sel.user.type.colname = "user.type",
	sub.sel.X.user.type = "immune",
	sub.sel.Y.user.type = "immune",
	sub.sel.user.merge.option = "union"
	)
# pics- id.1.right: shows the paracrine interaction inside these 3 clusters
GetResult.SummaryClustersInteracts(this.anal.immune, 
	show.clusters.in.x = c("CAF 1", "CAF 2", "Myeloid"), 
	show.clusters.in.y = c("CAF 1", "CAF 2", "Myeloid"), 
	power.max.limit = 580, hide.power.label = TRUE)$plot


###### deep analysis
# get the desired cluster group and their interaction relationship
test.caf1.tg <- ExtractTargetOnepairClusters(this.anal.immune, "CAF 1", "Myeloid")
# matching gene pairs to interaction action database
test.caf1.gmoc <- GenerateMapDetailOnepairClusters(test.caf1.tg, actions.mouse.ref.db)
# fetch the result
test.caf1.veinfo <- GenerateVEinfos(test.caf1.gmoc, fgenes.all, 
	direction.A.to.B = NULL, if.ignore.annos = FALSE)

# collect and grouping all gene pairs by expressiong changes and action effects
test.caf1.dgsa <- CollectHierachyOnepairClusters(test.caf1.veinfo, 
	limits.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"), 
	limits.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"))
# pics- id.2.left: show that action effects of most pairs are not well known
GetResult.SummaryOnepairClusters(test.caf1.dgsa, show.exprs.change = "Xup.Yup", 
	caption.label.x = 0.5, caption.label.y = 0.6)$plot
# pics- id.2.right: only see gene pairs with action effects being positive and negative
GetResult.SummaryOnepairClusters(test.caf1.dgsa, show.exprs.change = "Xup.Yup", 
	show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B"),
	caption.label.x = 0.5, caption.label.y = 0.6)$plot

# dicide to explore the positive relations, which is the A-->B and A<--B
test.caf1.op1.veinfo <- GenerateVEinfos(test.caf1.gmoc, fgenes.all, 
	direction.A.to.B = NULL, if.ignore.annos = FALSE)
test.caf1.op1.veinfo <- TrimVEinfos(test.caf1.op1.veinfo, 
	sel.mode.val = "activation",
	sel.action.effect.val = "positive",
	sel.mode.action.merge.option = "union")

# find special genes comparing to other cluster groups
test.caf1.spg <- FindSpecialGenesInOnepairCluster(test.caf1.op1.veinfo, this.anal.immune, 
	uq.cnt.range = 1:17)
# get the top 10 gene pairs
test.caf1.top10.spg <- GetResult.SummarySpecialGenes(test.caf1.spg, test.caf1.op1.veinfo, 
	show.uq.cnt.range = 1:17, show.uq.cnt.merged = TRUE, 
	select.genepairs.method = "logfc-sum", 
	select.by.method.decreasing = TRUE,  
	prioritize.cluster.groups = paste("CAF 1", "Myeloid", sep = kClustersSplit), 
	select.by.method.pairs.limit = 10, 
	show.topn.inside.onepair = 17)
# pics- id.3: show the top 10 gene pairs in specific expressing analysis
Tool.ShowGraph(test.caf1.top10.spg)

# get the top 60 gene pairs
test.caf1.summary.spg <- GetResult.SummarySpecialGenes(test.caf1.spg, test.caf1.op1.veinfo, 
	show.uq.cnt.range = 1:17, show.uq.cnt.merged = TRUE, 
	select.genepairs.method = "logfc-sum", 
	select.by.method.decreasing = TRUE,  
	prioritize.cluster.groups = paste("CAF 1", "Myeloid", sep = kClustersSplit), 
	select.by.method.pairs.limit = 60, 
	show.topn.inside.onepair = 17)
# extract the gene pairs
caf1.top.gene.pairs <- unique(test.caf1.summary.spg$table[, c("inter.GeneName.A", "inter.GeneName.B")])
# and further trim the VEinfos
test.caf1.topx.veinfo <- TrimVEinfos(test.caf1.op1.veinfo, 
	sel.some.gene.pairs.df = caf1.top.gene.pairs, 
	sel.some.gene.pairs.colnames = c("inter.GeneName.A", "inter.GeneName.B"))
# pics- id.4: show the top 60 gene pairs in the crosstalk graph
GetResult.PlotOnepairClusters.GeneCrosstalk(test.caf1.topx.veinfo, fgenes.all, 
	direction.A.to.B = TRUE,
	axis.order.xy = c("AlphaBet", "AlphaBet"), 
	axis.order.xy.decreasing = c(FALSE, FALSE), 
	range.to.use = list("LogFC" = c(-Inf, Inf), "p_val_adj" = c(-Inf, 40000)), 
	nodes.size.range = c(1, 6))$plot

# select 5 gene pairs to plot for
tmp.for.plot.df <- data.frame(act.gene.A = c("C3", "C3", "Cxcl12", "Cxcl12", "Csf1"),
	act.gene.B = c("C3ar1", "C5ar1", "Cxcr4", "Ccr1", "Csf1r"), stringsAsFactors = FALSE)
tmp.for.plot.select.genepairs <- paste(tmp.for.plot.df$act.gene.A, tmp.for.plot.df$act.gene.B, sep = kGenesSplit)
test.caf1.summary.spg.for.pics <- GetResult.SummarySpecialGenes(test.caf1.spg, test.caf1.op1.veinfo, 
	show.uq.cnt.range = 1:17, show.uq.cnt.merged = TRUE, 
	select.genepairs = tmp.for.plot.select.genepairs, 
	prioritize.cluster.groups = paste("CAF 1", "Myeloid", sep = kClustersSplit), 
	show.topn.inside.onepair = 17)
# pics- id.5: show the selected 5 gene pairs in specific expressing analysis
test.caf1.summary.spg.for.pics$plot

# extra plotting(not shown): show the gene pairs by their locations
GetResult.PlotOnepairClusters.CellPlot(test.caf1.topx.veinfo,
	area.extend.times = 15, 
	expand.gap.radius.list = list(ECM = 4, PM = 3, CTP = 2, OTHER = 2),
	expand.shift.degree.list = list(ECM = 90, PM = 90, CTP = 0, OTHER = 0), 
	expand.gap.degree.list = list(ECM = 180, PM = 180, CTP = 60, OTHER = 60),
	expand.center.density.list = list(ECM = 0.5, PM = 0.4, CTP = 0.25, OTHER = 0.25))$plot



# FIN
