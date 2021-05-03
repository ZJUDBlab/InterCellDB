
# This file is the workflow for analyzing scRNA-seq data on human cholangiocarcinoma
# The steps for fetching marker genes are omitted.

# load library
library(InterCellDB)

# read differentially expressed genes
tmp.markers <- read.csv("./DEGs-human-clusters.csv", stringsAsFactors = FALSE)
# select the subset of passing bonferroni correction
fgenes.all <- tmp.markers[which(tmp.markers$p_val_adj < 0.05), ]  # 43174 rows left

###### fullview analysis
# include 3 expression change status with at lease one part up-regulated
tmp.sub.sel.exprs.changes <- c("Xup.Yup", "Xup.Ydn", "Xdn.Yup")

# one more step - select the high confidence database
tmp.db.exp.only <- pairs.human.db[which(pairs.human.db$inter.Experiments.Score != 0), ]
tmp.db.rest.high.conf <- pairs.human.db[intersect(which(pairs.human.db$inter.Combined.Score >= 700),
	which(pairs.human.db$inter.Experiments.Score == 0)), ]
ref.data <- rbind(tmp.db.exp.only, tmp.db.rest.high.conf)

# analyze the cytokine-related procedure
this.anal.humanz1.v1 <- AnalyzeClustersInteracts(fgenes.all, 
	pairs.ref = ref.data,
	anno.location.ref = anno.location.human.ref.db,
	anno.type.ref = anno.type.human.ref.db,
	sub.sel.X.clusters = NULL,
	sub.sel.Y.clusters = NULL,
	sub.sel.exprs.changes = tmp.sub.sel.exprs.changes,
	sub.sel.X.Location = c("Extracellular Region"),
	sub.sel.X.Location.score.limit = c(3:5),
	sub.sel.Y.Location = c("Plasma Membrane"),
	sub.sel.Y.Location.score.limit = c(3:5),
	sub.sel.X.Type = "Cytokine"
	)
# pics- id.1.left: all interactions relating to cytokines
GetResult.SummaryClustersInteracts(this.anal.humanz1.v1)$plot
# pics- id.1.right: cancer cell associated interactions about cytokines
GetResult.SummaryClustersInteracts(this.anal.humanz1.v1, show.clusters.in.y = "Malignant")$plot


###### deep analysis
# get the desired cluster group and their interaction relationship
test.tg <- ExtractTargetOnepairClusters(this.anal.humanz1.v1, "Fibroblast", "Malignant")
# matching gene pairs to interaction action database
test.tg.gmoc <- GenerateMapDetailOnepairClusters(test.tg, actions.human.ref.db)
# fetch the result
test.tg.veinfo <- GenerateVEinfos(test.tg.gmoc, fgenes.all, 
	direction.A.to.B = NULL, if.ignore.annos = FALSE)

# collect and grouping all gene pairs by expressiong changes and action effects
test.tg.dgsa <- CollectHierachyOnepairClusters(test.tg.veinfo, 
	limits.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"), 
	limits.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"))
# pics- id.2.left: show gene pairs of all possible action effects 
GetResult.SummaryOnepairClusters(test.tg.dgsa, 
	show.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup"))$plot
# pics- id.2.right: only see gene pairs with action effects being positive and negative
GetResult.SummaryOnepairClusters(test.tg.dgsa, 
	show.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup"), 
	show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B"))$plot

# dicide to explore the positive relations, which is the A-->B and A<--B
test.op1.veinfo <- GenerateVEinfos(test.tg.gmoc, fgenes.all, 
	direction.A.to.B = TRUE, if.ignore.annos = FALSE)
test.op1.veinfo <- TrimVEinfos(test.op1.veinfo,
	sel.mode.val = "activation",
	sel.action.effect.val = "positive",
	sel.mode.action.merge.option = "union")

# find special genes comparing to other cluster groups
test.caf.tumor.spg <- FindSpecialGenesInOnepairCluster(test.op1.veinfo, this.anal.humanz1.v1, 
	uq.cnt.range = 1:8)
# get the top 10 gene pairs
test.caf.tumor.for.plot <- GetResult.SummarySpecialGenes(test.caf.tumor.spg, test.op1.veinfo, 
	show.uq.cnt.range = 1:8, show.uq.cnt.merged = TRUE, 
	select.genepairs.method = "logfc-sum",
	select.by.method.decreasing = TRUE, 
	prioritize.cluster.groups = paste("Fibroblast", "Malignant", sep = kClustersSplit), 
	select.by.method.pairs.limit = 10, 
	show.topn.inside.onepair = 8)
# pics- id.3: show the top 10 gene pairs in specific expressing analysis
test.caf.tumor.for.plot$plot

# count on cytokine and get the ranked cytokines
tmp.cytokine.kinds <- test.caf.tumor.for.plot$table[which(test.caf.tumor.for.plot$table$cluster.groups == paste("Fibroblast", "Malignant", sep = kClustersSplit)), ]
res.cytokine.kinds <- tapply(seq_len(nrow(tmp.cytokine.kinds)), tmp.cytokine.kinds[, "inter.GeneName.A"], length)

# get the top 60 gene pairs
test.caf.tumor.summary.spg <- GetResult.SummarySpecialGenes(test.caf.tumor.spg, test.op1.veinfo, 
	show.uq.cnt.range = 1:8, show.uq.cnt.merged = TRUE, 
	select.genepairs.method = "logfc-sum",
	select.by.method.decreasing = TRUE, 
	prioritize.cluster.groups = paste("Fibroblast", "Malignant", sep = kClustersSplit), 
	select.by.method.pairs.limit = 60, 
	show.topn.inside.onepair = 8)
# extract the gene pairs
caf.tumor.top.gene.pairs <- unique(test.caf.tumor.summary.spg$table[, c("inter.GeneName.A", "inter.GeneName.B")])
# and further correctly construct the VEinfos
test.caf.tumor.topx.veinfo <- GenerateVEinfos(test.tg.gmoc, fgenes.all, 
	direction.A.to.B = TRUE, if.ignore.annos = FALSE)
test.caf.tumor.topx.veinfo <- TrimVEinfos(test.caf.tumor.topx.veinfo, 
	sel.some.gene.pairs.df = caf.tumor.top.gene.pairs, 
	sel.some.gene.pairs.colnames = c("inter.GeneName.A", "inter.GeneName.B"))
# pics- id.4: plot all top 60 gene pairs by their subcellular locations and 
# get to directly see the distribution of gene pairs
GetResult.PlotOnepairClusters.CellPlot(test.caf.tumor.topx.veinfo,
	area.extend.times = 40, 
	nodes.label.size = 2.5, 
	nodes.label.colour = "white", 
	expand.gap.radius.list = list(ECM = 4, PM = 3, CTP = 2, OTHER = 2),
	expand.shift.degree.list = list(ECM = 90, PM = 90, CTP = 0, OTHER = 0), 
	expand.gap.degree.list = list(ECM = 180, PM = 180, CTP = 60, OTHER = 60),
	expand.center.density.list = list(ECM = 0.5, PM = 0.5, CTP = 0.25, OTHER = 0.25))$plot
# pics- id.5: show the top 60 gene pairs in the crosstalk graph
GetResult.PlotOnepairClusters.GeneCrosstalk(test.caf.tumor.topx.veinfo, fgenes.all, 
	direction.A.to.B = TRUE,
	axis.order.xy = c("AlphaBet", "AlphaBet"), 
	axis.order.xy.decreasing = c(FALSE, FALSE), 
	range.to.use = list("LogFC" = c(-Inf, Inf), "p_val_adj" = c(-Inf, 400000)), 
	plot.font.size.base = 8, 
	nodes.size.range = c(0.5, 5.5))$plot



# FIN
