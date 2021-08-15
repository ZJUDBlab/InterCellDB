
# This file is to compare InterCellDB with other 5 databases
library(dplyr)
library(ggplot2)
library(cowplot)

# global variables
# DB are arranged by presenting time
kfDB.namelist <- c("iTALK", "CellPhoneDB", "SingleCellSignalR", "NicheNet", "CellChatDB")
kfDB.favor.color <- c("#A1D99B", "#C6DBEF", "#74C476", "#9ECAE1", "#6BAED6")
#kfDB.favor.color <- c("#D7C0BA", "#CEC3F2", "#CDD3D4", "#96CAEA", "#8CACAF")
#kfDB.favor.color <- c("#B2DF8A", "#A6CEE3", "#33A02C","#CAB2D6","#1F78B4")

# InterCellDB are ordered BY # experimentally validated # from curated database # by prediction # All included
kfIT.namelist <- c("InterCellDB.exp", "InterCellDB.know", "InterCellDB.pred", "InterCellDB", "InterCellDB.LR")
kfIT.favor.color <- c("#C19E77", "#E6CCA9", "#E9BC80", "#8E7152",  "#D7C0BA")  # choose from ggplot2::scale_color_distiller(palette=18 <or> 8)
	# reorder IT
	kfIT.order.inds <- c(4, 1:3, 5)
	kfIT.namelist <- kfIT.namelist[kfIT.order.inds]
	kfIT.favor.color <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#FDB863")
	#kfIT.favor.color <- rev(ref.brewer.Reds[3:7])
# backup colours(dark red to yellowish red): #BC0625, #D6131F, #F03623, #F0642E, #FE9840, #FEB04B, #FFCD69, #FFE88E
# backup colours(light grey blue to deep grey blue): #1496d4, #1f90c6, #2789b8, #2f81a9, #357a9a, #3b6f88, #3a6377, #375563
# backup colour(light grey blue): #96CAEA
# backup colour(brown yelow): #D9C5A4
# backup colour(dark green brown): #8E7152; a bit light: #A6896B
ref.color <-c("#D51344", "#EE5B9A", "#F16E5A", "#FBB832","#18933E","#5C7A28","#05A3AA","#34A2DC","#696DA7","#7D5786")
ref.brewer.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
ref.redish.color <- c("#5B1D42", "#B85186", "#7F0206", "#CB1109", "#F52E27", "#C5344C")
ref.brewer.Reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
ref.brewer.OrRd <- c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000")
ref.brewer.BrBG10 <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30")
ref.select.v1.color <- c(
	"#9D1742",  # purple red
	"#D53E4E",  # fair red
	"#C5247D",  # dark pink
	"#F46D42",  # orange-ish
	"#762A83"   # fair purple
)
ref.select.v2.color <- c(
	"#FDE0EF",   # light pink
	"#E5F5D1",   # light green
	"#DFF3F8",   # light blue
	"#FEE0B6",   # light orange
	"#E7D4E7"    # light purple
)
ref.select.v3.color <- c(
	"#F1B6DA",   # light-deepen pink
	"#B8E086",   # light-deepen green
	"#C2A5CF",   # light-deepen purple
	"#DFC27D",   # light-deepen brown
	"#91C4DE"    # light-deepen blue
)
ref.select.v4.color <- c(
	"#762A83",  # fair purple --> ALL 
	"#D53E4E",  # fair red --> Exp
	"#3388BD",  # fair blue 
	"#1B7838",  # fair green
	"#C5247D"   # fair pink
) 
blues.select.v5.color <- c(
	"#CBDFE5",  # mk blue
	"#A5C9D3",  # mk deep1 blue
	"#6B91A1",  # mk deep2 blue
	"#93CBED",  # mk fair blue -- OK
	"#4A6CA4"   # mk dark blue
)


#
#
kfDB.loc.all <- c("Extracellular Region", "Plasma Membrane", 
	"Cytoplasm", "Cytosol", "Cytoskeleton", 
	"Endoplasmic Reticulum", "Golgi Apparatus", 
	"Endosome", "Lysosome", "Peroxisome",
	"Mitochondrion", "Nucleus", "Other")
kfDB.type.all <- c("Actin-binding and motor protein", "Antimicrobial", "Cytokine", "Enzyme",
	"Growth Factor", "Hormone", "Hydrolase", "Ion Channel", "Other", "Protein metabolism regulator",
	"Receptor", "Transcription regulator", "Transducer", "Transferase", "Translation regulator", "Vasoactive")

# actions
kfDB.act.mode <- c("activation", "binding", "catalysis", "expression", "inhibition", "ptmod", "reaction")
kfDB.act.effect <- c("positive", "negative", "unspecified", "undirected")

# picked GO terms
kfGO.picked.terms <- c("angiogenesis", "EMT", "complement activation", "cytokine activity")
	## GO:0001525 [angiogenesis] <BP>
	go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0001837 [epithelial to mesenchymal transition] <BP>
	go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0006955 [immune response]
	#go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.human.ref.db, go.human.ref.db, 
	#	go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0006956 [complement activation] <BP>
	go.comp.act.genes <- Tool.FindGenesFromGO("GO:0006956", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	## GO:0005125 [cytokine activity] <MF>
	go.cytokine.act.genes <- Tool.FindGenesFromGO("GO:0005125", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
kfGO.picked.list <- list(go.vessel.genes, go.emt.genes, go.comp.act.genes, go.cytokine.act.genes)
names(kfGO.picked.list) <- kfGO.picked.terms

# backup GO terms
c("cell chemotaxis")  # <BP>
# <MF>
c("cytokine activity", "GO:0005125")
c("receptor ligand activity", "GO:0048018")
c("growth factor activity", "GO:0008083")
# <CC>
c("plasma membrane protein complex", "GO:0098797")
c("extracellular matrix", "GO:0031012")
c("receptor complex", "GO:0043235")


# - picked GO terms end - 

# tool functions for main ctrl or accessory files
#
kfgp.cut <- ">"  # split for 2 genes in gene pairs
#


### Brief explanation
# [part xxx] is incremental, and never removed. Any new idea will get one new part ID.
# [deter <a-z>] is to denote finally picked ones.



##### For Describing GENE - Human
#
# 
tmp.genes.cnt.all <- list(italk.genes.res, cellphone.genes.res, scsr.genes.res, nichenet.genes.res, cellchat.human.genes.res, 
	intercelldb.exp.human.genes.res, intercelldb.know.human.genes.res, intercelldb.pred.human.genes.res, 
	intercelldb.human.genes.res, intercelldb.LR.human.genes.res)
names(tmp.genes.cnt.all) <- c(kfDB.namelist, kfIT.namelist)


## [fig2-A]: location based gene count, and core genes type distribution
	# tool function
	Ctrl.replace.locs.cytoplasm <- function(need.replace) {
		res.to.replace <- "Cytoplasm"
		if (need.replace == TRUE) {
			res.to.replace <- c(res.to.replace, c("Cytoplasm", "Cytoskeleton", "Cytosol", "Endoplasmic Reticulum", 
				"Endosome", "Endosome", "Golgi Apparatus", "Lysosome", "Mitochondrion", "Peroxisome"))
		}
		return(res.to.replace)	
	}
	#
	# Location Part
	#
	check.loc.each.gene.most.conf <- tapply(anno.location.human.ref.db$score, anno.location.human.ref.db$Gene.name, max)
	#   2    3    4    5   run at 2021.02.23
	# 106 3086 4482 9972
	# so we use the most confident locations for every gene by score selection
	check.loc.egm.df <- data.frame(gene = names(check.loc.each.gene.most.conf), loc.max.score = check.loc.each.gene.most.conf, stringsAsFactors = FALSE)
	check.loc.egm.df <- left_join(check.loc.egm.df, anno.location.human.ref.db[, c("GeneID", "Gene.name", "GO.Term.target", "score")],
		by = c("gene" = "Gene.name", "loc.max.score" = "score"))
	check.loc.egm.df <- unique(check.loc.egm.df)  # 28197 rows
	# 
	tmp.use.loc.category <- kfDB.loc.all  # change to use for different purpose
	tmp.use.loc.replace <- FALSE
	tmp.use.DB.namelist <- c(kfDB.namelist, kfIT.namelist)
	#
	tmp.check.loc.alldb <- lapply(tmp.genes.cnt.all[tmp.use.DB.namelist], 
		ref.loc.kinds = tmp.use.loc.category, gene.loc.ref.db = check.loc.egm.df, 
		need.replace = tmp.use.loc.replace, 
		function(x, ref.loc.kinds, gene.loc.ref.db, need.replace) {
			as.integer(unlist(lapply(ref.loc.kinds, use.genes = x, gene.loc.ref.db = gene.loc.ref.db, 
				function(y, use.genes, gene.loc.ref.db) {
					if (y == "Cytoplasm") {
						y <- Ctrl.replace.locs.cytoplasm(need.replace = need.replace)
					}
					this.ref.genes <- gene.loc.ref.db[which(gene.loc.ref.db$GO.Term.target %in% y), "gene"]
					tmp.ret <- length(intersect(use.genes, this.ref.genes))
					#c(tmp.ret, tmp.ret / length(use.genes))
					tmp.ret
				})))
		})
	tmp.check.loc.alldb <- as.numeric(unlist(tmp.check.loc.alldb))
	tmp.loc.anal.df <- data.frame(DBname = rep(tmp.use.DB.namelist, each = length(tmp.use.loc.category)),
		Location = rep(tmp.use.loc.category, times = length(tmp.use.DB.namelist)),
		count = as.integer(unlist(tmp.check.loc.alldb)),
		stringsAsFactors = FALSE)
	#tmp.loc.anal.df$count <- sapply(tmp.loc.anal.df$count, function(x) {if (x > 1000) x <- 1000;x}, USE.NAMES = FALSE)
	tmp.loc.anal.df$count <- sapply(tmp.loc.anal.df$count, function(x) {log2(x)}, USE.NAMES = FALSE)
	tmp.loc.anal.df$count <- sapply(tmp.loc.anal.df$count, function(x) {if (x < 5) x <- 5; x}, USE.NAMES = FALSE)
	# remove InterCellDB
	#tmp.loc.anal.df <- tmp.loc.anal.df[which(tmp.loc.anal.df$DBname != "InterCellDB"), ]
	# add factor
	tmp.loc.anal.df$Location <- factor(tmp.loc.anal.df$Location, levels = tmp.use.loc.category)
	tmp.loc.anal.df$DBname <- factor(tmp.loc.anal.df$DBname, levels = tmp.use.DB.namelist)
	tmp.loc.anal.plot <- ggplot(data = tmp.loc.anal.df) + 
			geom_raster(aes(x = DBname, y = Location, fill = count)) + 
			scale_fill_gradientn(name = "Count", colours = c("#00809D", "#EEEEEE", "#C30000"), values = c(0.0, 0.5, 1.0)) + 
			theme_cowplot(16) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	# (width = 6, height = 6.5)
	# (width = 8, height = 6.5)

#
# [Type Part]
#
	tmp.use.mf.unip.sel <- setdiff(unique(Uniprot.key.map.list[,2]), "Other")
	# use Uniprot remap to mapping those keys to get its subset of Keywords
	tmp.genes.allowed.from.loc <- check.loc.egm.df[which(check.loc.egm.df$GO.Term.target %in% c("Extracellular Region", "Plasma Membrane")), "gene"]
	tmp.use.mf.unip.genelist <- lapply(tmp.use.mf.unip.sel, 
		ref.unip.remap = Uniprot.key.map.list,
		ref.unip.anno.db = anno.type.human.ref.db, 
		ref.loc.allowed.genes = tmp.genes.allowed.from.loc,  # [here included to only see the ECR and PM genes]
		function(x, ref.unip.remap, ref.unip.anno.db, ref.loc.allowed.genes) {
			this.unip.keys <- ref.unip.remap[which(ref.unip.remap[, 2] == x), 1]
			this.unip.sel.genes <- unique(ref.unip.anno.db[which(ref.unip.anno.db$Keyword.Name %in% this.unip.keys), "Gene.name"])
			unique(intersect(this.unip.sel.genes, ref.loc.allowed.genes))
		})
	names(tmp.use.mf.unip.genelist) <- tmp.use.mf.unip.sel
	# run Uniprot Keywords related genes and DB genes intersection
	tmp.use.DB.namelist.for.unip <- c(kfDB.namelist, kfIT.namelist)
	tmp.DBgenes.unip.countlist <- lapply(seq_along(tmp.genes.cnt.all[tmp.use.DB.namelist.for.unip]),
		inside.genes.cnt.all = tmp.genes.cnt.all[tmp.use.DB.namelist.for.unip], unip.list = tmp.use.mf.unip.genelist,
		function(x, inside.genes.cnt.all, unip.list) {
			this.db.genes <- inside.genes.cnt.all[[x]]
			this.db.unip.inter <- vapply(unip.list, it.genes = this.db.genes,
				FUN.VALUE = integer(1),
				USE.NAMES = FALSE,
				FUN = function(y, it.genes) {
					length(intersect(y, it.genes))
					})
			data.frame(DBname = rep(names(inside.genes.cnt.all)[x], times = length(unip.list)),
				Unip.Keywords = names(unip.list),
				unip.cnt = this.db.unip.inter,
				stringsAsFactors = FALSE)
		})
	tmp.DBgenes.unip.countdf <- bind_rows(tmp.DBgenes.unip.countlist)
	tmp.DBgenes.unip.countdf$DBname <- factor(tmp.DBgenes.unip.countdf$DBname, levels = tmp.use.DB.namelist.for.unip)
	tmp.DBgenes.unip.countdf$log.unip.cnt <- log2(tmp.DBgenes.unip.countdf$unip.cnt)
	# check diff, and get it ordered
	tmp.unip.order.ref <- abs(as.numeric(unlist(lapply(tmp.use.mf.unip.sel, 
		ref.df = tmp.DBgenes.unip.countdf,
		function(x, ref.df) {
			tmp.use.it.df <- tmp.DBgenes.unip.countdf[which(tmp.DBgenes.unip.countdf$Unip.Keywords == x), ]
			# implicitly 5 vs 5
			mean(tmp.use.it.df[6:10, "log.unip.cnt"]) - mean(tmp.use.it.df[1:5, "log.unip.cnt"])
		}))))
	names(tmp.unip.order.ref) <- tmp.use.mf.unip.sel
	tmp.unip.order.ref <- tmp.unip.order.ref[order(tmp.unip.order.ref)]
	tmp.DBgenes.unip.countdf$Unip.Keywords <- factor(tmp.DBgenes.unip.countdf$Unip.Keywords, levels = names(tmp.unip.order.ref))
	# plot settings
	tmp.for.unip.colour <- c(kfDB.favor.color, kfIT.favor.color)
	names(tmp.for.unip.colour) <- c(kfDB.namelist, kfIT.namelist)

	#tmp.unip.plot.group.1 <- levels(tmp.DBgenes.unip.countdf$Unip.Keywords)[1:5]
	#tmp.unip.plot.group.2 <- levels(tmp.DBgenes.unip.countdf$Unip.Keywords)[6:10]
	#tmp.unip.plot.group.3 <- levels(tmp.DBgenes.unip.countdf$Unip.Keywords)[11:15]

	# plot
	#tmp.unip.plot.use.group <- tmp.unip.plot.group.1  # [NOTE] select one group once
	tmp.unip.plot.use.df <- tmp.DBgenes.unip.countdf#[which(tmp.DBgenes.unip.countdf$Unip.Keywords %in% tmp.unip.plot.use.group), ]
	tmp.unip.plot.use.df$log.unip.cnt[which(tmp.unip.plot.use.df$log.unip.cnt < 0)] <- 0
	tmp.for.unip.plot <- ggplot(data = tmp.unip.plot.use.df) + 
		geom_col(aes(x = DBname, y = log.unip.cnt, group = Unip.Keywords, fill = DBname),
			width = 0.4) + 
		facet_wrap(.~Unip.Keywords, nrow = 3) + 
		scale_fill_manual(values = tmp.for.unip.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			strip.background = element_blank(), 
			strip.text.x = element_text(size = 5),
			legend.position = "right")
	# (width = 10, height = 6.5)
	# fetch legend
	tmp.for.unip.for.leg.plot <- ggplot(data = tmp.unip.plot.use.df) + 
		geom_col(aes(x = DBname, y = log.unip.cnt, group = Unip.Keywords, fill = DBname),
			width = 0.4) + 
		facet_grid(.~Unip.Keywords) + 
		scale_fill_manual(values = tmp.for.unip.colour) + 
		scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
		#scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_blank(),  # element_text(angle = 90, hjust = 1, vjust = 0.5)
			axis.ticks.x = element_blank(), 
			strip.text.x = element_text(size = 5),
			legend.position = "right")
	# (width = 8, height = 7)





##### For Describing GENE PAIRS - Human
# pasted gene pairs
cmp.pairs.db.list <- list(remap.italk.pairs.human, remap.cellphone.pairs.human, remap.scsr.pairs.human, 
	remap.nichenet.pairs.human, remap.cellchat.pairs.human,
	ref.intercelldb.pairs.human.experiments, ref.intercelldb.pairs.human.knowledge, 
	ref.intercelldb.pairs.human.prediction, ref.intercelldb.pairs.human, 
	ref.intercelldb.pairs.human.LR)
cmp.pairs.db.list[6:10] <- cmp.pairs.db.list[6:10][kfIT.order.inds]
names(cmp.pairs.db.list) <- c(kfDB.namelist, kfIT.namelist)

# position-aligned gene pairs
pos.pairs.db.list <- list(pos.italk.pairs.human, pos.cellphone.pairs.human, pos.scsr.pairs.human, 
	pos.nichenet.pairs.human, pos.cellchat.pairs.human,
	pos.intercelldb.pairs.human.experiments, pos.intercelldb.pairs.human.knowledge,
	pos.intercelldb.pairs.human.prediction, pos.intercelldb.pairs.human,
	pos.intercelldb.pairs.human.LR)
pos.pairs.db.list[6:10] <- pos.pairs.db.list[6:10][kfIT.order.inds]
names(pos.pairs.db.list) <- c(kfDB.namelist, kfIT.namelist)


## [fig2-B]: ribbon plot for showing the differences of gene pairs in every database
# fetch all intersect-ed ligand and receptor genes (genes that exist in all database and identified as ligand or receptor)
	# receptor
	tmp.receptor.cnt.all <- list(sel.receptor.italk.res, sel.receptor.cpdb.res, sel.receptor.scsr.res, 
		sel.receptor.nichenet.res, sel.receptor.cellchat.res, sel.receptor.intercelldb.res)
	names(sel.receptor.cnt.all) <- c(kfDB.namelist, "InterCellDB")
	tmp.allagree.receptor.genes <- Reduce(intersect, tmp.receptor.cnt.all)
	# ligand
	tmp.ligand.cnt.all <- list(sel.ligand.italk.res, sel.ligand.cpdb.res, sel.ligand.scsr.res, 
		sel.ligand.nichenet.res, sel.ligand.cellchat.res, sel.ligand.intercelldb.res)
	names(tmp.ligand.cnt.all) <- c(kfDB.namelist, "InterCellDB")
	tmp.allagree.ligand.genes <- Reduce(intersect, tmp.ligand.cnt.all)
	# merged genes
	tmp.for.gp.allagree.genes <- union(tmp.allagree.receptor.genes, tmp.allagree.ligand.genes)


	# collecting data process
	prog.bar.ribbon.gp <- progress::progress_bar$new(total = length(tmp.for.gp.allagree.genes))
	prog.bar.ribbon.gp$tick(0)
	tmp.ribbon.gp.matchcnt.raw.df <- bind_rows(lapply(tmp.for.gp.allagree.genes, 
		allDB.list = pos.pairs.db.list,
		function(x, allDB.list) {
			tmp.DBnames <- names(allDB.list)
			tmp.match.res <- as.integer(unlist(lapply(allDB.list, onegene = x,
				function(y, onegene) {
					length(unique(c(which(y$gene.A == onegene), which(y$gene.B == onegene))))
				})))
			prog.bar.ribbon.gp$tick()
			data.frame(DBname = tmp.DBnames,
				Gene = rep(x, times = length(tmp.DBnames)),
				pairs.cnt = tmp.match.res,
				log.pairs.cnt = log2(tmp.match.res),
				stringsAsFactors = FALSE)
		})
	)

	# get order of genes and set each genes as one-increment-step numbers
	# ordered BY  the count of gene pairs for every gene, in InterCellDB (ALL)
	tmp.ribbon.order.ref.df <- tmp.ribbon.gp.matchcnt.raw.df[which(tmp.ribbon.gp.matchcnt.raw.df$DBname == "InterCellDB"), ]
	tmp.ribbon.order.ref.df <- tmp.ribbon.order.ref.df[order(tmp.ribbon.order.ref.df$pairs.cnt, decreasing = TRUE), ]
	# get order data.frame
	tmp.ribbon.gene.order <- tmp.ribbon.order.ref.df[, "Gene"]
	tmp.ribbon.gene.order.df <- data.frame(Gene = tmp.ribbon.gene.order, num.gene = seq_along(tmp.ribbon.gene.order), stringsAsFactors = FALSE)
	#
	tmp.ribbon.gp.matchcnt.df <- left_join(tmp.ribbon.gp.matchcnt.raw.df, tmp.ribbon.gene.order.df, by = c("Gene" = "Gene"))
	# replace log -Inf to 0
	tmp.ribbon.gp.matchcnt.df$log.pairs.cnt[which(tmp.ribbon.gp.matchcnt.df$log.pairs.cnt < 0)] <- 0
	# plot prep
	tmp.ribbon.gp.matchcnt.df$DBname <- factor(tmp.ribbon.gp.matchcnt.df$DBname, levels	= names(pos.pairs.db.list))
	tmp.ribbon.db.colour <- c(kfDB.favor.color, kfIT.favor.color)
	names(tmp.ribbon.db.colour) <- c(kfDB.namelist, kfIT.namelist)
	# top genes
	tmp.top10.ribbon.df <- tmp.ribbon.gp.matchcnt.df[intersect(which(tmp.ribbon.gp.matchcnt.df$DBname == "InterCellDB"), which(tmp.ribbon.gp.matchcnt.df$num.gene <= 10)), ]
	# plot
	tmp.ribbon.gp.plot <- ggplot(data = tmp.ribbon.gp.matchcnt.df) + 
		geom_point(aes(x = num.gene, y = log.pairs.cnt, group = DBname, colour = DBname), 
			size = 0.3, alpha = 0.5) + 
		geom_smooth(aes(x = num.gene, y = log.pairs.cnt, group = DBname, colour = DBname, fill = DBname),
			method = "loess", formula = y ~ x) + 
		geom_text(data = tmp.top10.ribbon.df, aes(x = num.gene, y = log.pairs.cnt * 1.05, label = Gene),
			size = 1) + 
		scale_colour_manual(values = tmp.ribbon.db.colour) + 
		scale_fill_manual(values = tmp.ribbon.db.colour) + 
		#geom_text(data = tmp.top10.receptor.plot.df[which(tmp.top10.receptor.plot.df$DBname == "InterCellDB"), ], aes(x = num.gene, y = pairs.log.cnt, label = Gene),
		#	nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		theme_cowplot(12)
	# (width = 10, height = 8)
	# only ribbon
	tmp.ribbon.null.gp.plot <- ggplot(data = tmp.ribbon.gp.matchcnt.df) + 
		geom_point(aes(x = num.gene, y = log.pairs.cnt, group = DBname, colour = DBname), size = 1) + 
		geom_smooth(aes(x = num.gene, y = log.pairs.cnt, group = DBname, colour = DBname, fill = DBname),
			method = "loess", formula = y ~ x) + 
		scale_colour_manual(values = tmp.ribbon.db.colour) + 
		scale_fill_manual(values = tmp.ribbon.db.colour) + 
		theme_void() + 
		theme(legend.position = "none")
	# (width = 7.5, height = 7.5)
	# only axis
	tmp.ribbon.only.axis.gp.plot <- ggplot(data = tmp.ribbon.gp.matchcnt.df) + 
		geom_smooth(aes(x = num.gene, y = log.pairs.cnt, group = DBname, colour = DBname, fill = DBname),
			method = "loess", formula = y ~ x) + 
		scale_colour_manual(values = tmp.ribbon.db.colour) + 
		scale_fill_manual(values = tmp.ribbon.db.colour) + 
		#geom_text(data = tmp.top10.receptor.plot.df[which(tmp.top10.receptor.plot.df$DBname == "InterCellDB"), ], aes(x = num.gene, y = pairs.log.cnt, label = Gene),
		#	nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		theme_cowplot(12)
	# (width = 10, height = 8)



## [fig2-C]: overlap of database in gene pairs, use 5 published as reference, 
##           check the rate how InterCellDB covers those
# fetch ref and cmp database list
	tmp.overlap.ref.gp.db.list <- cmp.pairs.db.list[kfDB.namelist]
	tmp.overlap.tocmp.gp.db.list <- cmp.pairs.db.list[kfIT.namelist]
	# co-relation matrix in list format
	tmp.gpairs.corel.list <- list()
	for (i in seq_along(tmp.overlap.ref.gp.db.list)) {  # goes from y, pick one DB-Y each time
		this.dummy.name <- names(tmp.overlap.ref.gp.db.list)[i]
		for (j in seq_along(tmp.overlap.tocmp.gp.db.list)) {
			tmp.cmp.name <- names(tmp.overlap.tocmp.gp.db.list)[j]
			tmp.inter.cnt <- length(intersect(tmp.overlap.ref.gp.db.list[[i]], tmp.overlap.tocmp.gp.db.list[[j]]))
			tmp.gpairs.corel.list <- c(tmp.gpairs.corel.list, list(data.frame(DBx = tmp.cmp.name, DBy = this.dummy.name, 
				inter.cnt = tmp.inter.cnt,
				inter.percent = tmp.inter.cnt / length(tmp.overlap.ref.gp.db.list[[i]]),
				stringsAsFactors = FALSE)))
		}
	}
	tmp.gpairs.corel.df <- bind_rows(tmp.gpairs.corel.list)	
	# plot data prep
	tmp.gpairs.corel.df$DBx <- factor(tmp.gpairs.corel.df$DBx, levels = kfIT.namelist)
	tmp.gpairs.corel.df$DBy <- factor(tmp.gpairs.corel.df$DBy, levels = kfDB.namelist)	
	tmp.gpairs.corel.df$inter.percent <- round(tmp.gpairs.corel.df$inter.percent, digits = 3)
	tmp.gpairs.corel.df$inter.percent.label <- format(tmp.gpairs.corel.df$inter.percent, scientific = FALSE)
	# plot
	tmp.db.overlap.plot <- ggplot(data = tmp.gpairs.corel.df) + 
		geom_raster(aes(x = DBx, y = DBy, fill = inter.percent)) + 
		geom_text(aes(x = DBx, y = DBy, label = inter.percent.label)) + 
		scale_fill_gradientn(name = "cover ratio", colours = c("#00809D", "#EEEEEE", "#C30000"), values = c(0.0, 0.5, 1.0)) + 
		theme_cowplot(16) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	# (width = 6.5, height = 5)






# [fig2-D]: GO terms related genes and interactions with themselves
tmp.sel.GOterms.sml <- list(
	# <BP>
	c("aging", "GO:0007568"),  # about 3k records
	#c("cell aging", "GO:0007569"), 
	c("angiogenesis", "GO:0001525"),
	c("epithelial to mesenchymal transition", "GO:0001837"),  # EMT
	c("complement activation", "GO:0006956"),
	c("cell chemotaxis", "GO:0060326"),
	c("acute inflammatory response", "GO:0002526")
	)
tmp.sel.GOterms.big <- list(
	c("cell-cell signaling", "GO:0007267"), # 28k
	c("immune response", "GO:0006955"),  # 30k
	#c("cytokine-mediated signaling pathway", "GO:0019221"),
	#c("regulation of growth", "GO:0040008"),  # 10k
	c("cell development", "GO:0048468")  # 40k
	)

# subset the GOlist
tmp.chk.genes.picked.GOlist <- tmp.sel.GOterms.sml
#
# get genes
tmp.chk.genes.picked.goid <- vapply(tmp.chk.genes.picked.GOlist, function(x) {x[2]}, FUN.VALUE = character(1))
tmp.picked.go.rel.genelist <- lapply(tmp.chk.genes.picked.goid, 
	function(x) {
		Tool.FindGenesFromGO(x, genes.human.ref.db, go.human.ref.db, 
			go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	})
names(tmp.picked.go.rel.genelist) <- vapply(tmp.chk.genes.picked.GOlist, function(x) {x[1]}, FUN.VALUE = character(1))
# check GOterm related gene pairs
prog.bar.gorel.gp <- progress::progress_bar$new(total = length(tmp.picked.go.rel.genelist))
prog.bar.gorel.gp$tick(0)
tmp.picked.go.res.df <- bind_rows(lapply(seq_along(tmp.picked.go.rel.genelist), 
	allDB.list = pos.pairs.db.list, go.rel.list = tmp.picked.go.rel.genelist,
	function(i, allDB.list, go.rel.list) {
		this.go.term <- names(go.rel.list)[i]
		this.go.rel.dbs.res <- as.integer(unlist(lapply(allDB.list, go.rel.genes = go.rel.list[[i]],
			function(y, go.rel.genes) {
				length(intersect(which(y$gene.A %in% go.rel.genes), which(y$gene.B %in% go.rel.genes)))
			})))
		prog.bar.gorel.gp$tick()
		data.frame(DBname = names(allDB.list),
			GOterm = rep(this.go.term, times = length(allDB.list)),
			pairs.cnt = this.go.rel.dbs.res,
			log.pairs.cnt = log2(this.go.rel.dbs.res),
			stringsAsFactors = FALSE)
	}))
# plot data prep
tmp.picked.go.res.df$DBname <- factor(tmp.picked.go.res.df$DBname, levels = c(kfDB.namelist, kfIT.namelist))
tmp.picked.go.res.df$GOterm <- factor(tmp.picked.go.res.df$GOterm, levels = names(tmp.picked.go.rel.genelist))

# only used in sml
#
#tmp.picked.go.res.df$log.pairs.cnt <- sapply(tmp.picked.go.res.df$log.pairs.cnt, function(x) {if (x > 12) x <- 12; x})
#


# plot
tmp.picked.go.plot <- ggplot(data = tmp.picked.go.res.df) + 
	geom_raster(aes(x = GOterm, y = DBname, fill = log.pairs.cnt)) + 
	scale_fill_gradientn(name = "log(count)", colours = c("#00809D", "#EEEEEE", "#C30000"), values = c(0.0, 0.5, 1.0)) + 
	theme_cowplot(16) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# (width = 7, height = 7) sml
# (width = 5.5, height = 7) big



