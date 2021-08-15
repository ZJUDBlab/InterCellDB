#
# This file records the workflow of mouse
# This file inherits all the function and global variables from Main-ctrl-process.R

#
# global variables [inherited]
# kDB.namelist <- c("CellPhoneDB", "NicheNet", "iTALK", "SingleCellSignalR", "CellChatDB", "InterCellDB")
# kDB.favor.color <- c("#CEC3F2", "#96CAEA", "#EF98CE", "#D9C5A4", "#8CACAF", "#BCC3C4")
#


##### For Describing GENE - Mouse
# --- 1st --- 
#   Extract plain genes for every DB containing mouse specific gene and gene pairs
tmp.mouse.genes.cnt.all <- list(cellchat.mouse.genes.res, intercelldb.mouse.genes.res)
names(tmp.mouse.genes.cnt.all) <- kDB.namelist[5:6]
# part0: directly show all intersections of all DBs
pdf(file = "../anal-peer-db-Result/GeneItself/genes-mouse-allDBs-redundancy.pdf", 
	width = 9, height = 6)
plot(euler(tmp.mouse.genes.cnt.all), fills = kDB.favor.color[5:6], 
	edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
	labels = list(fontfamily = "sans", cex = 0.8), 
	quantities = list(cex = 0.8), 
	legend = TRUE)
dev.off()

# part x: give the genes overlap percentage by referring InterCellDB
tmp.mouse.genes.overlap.list.fx <- list()
for (i in 1) {
	tmp.mouse.genes.overlap.list.fx <- c(tmp.mouse.genes.overlap.list.fx, list(Ctrl.percent.plot.add.data(tmp.mouse.genes.cnt.all[[i]], names(tmp.mouse.genes.cnt.all)[i],
		col.desc.vals = c("matched.genes", "unmatched.genes"), 
		ref.vector = list(intercelldb.mouse.genes.res),
		count.or.identity = TRUE, 
		show.unmatches = TRUE)))
}
tmp.mouse.genes.overlap.df.fx <- dplyr::bind_rows(tmp.mouse.genes.overlap.list.fx)
tmp.mouse.genes.overlap.df.fx$target.names <- factor(tmp.mouse.genes.overlap.df.fx$target.names, levels = names(tmp.mouse.genes.cnt.all)[1])
# --- plot ---
pic.mouse.genes.overlap.fx <- ggplot(data = tmp.mouse.genes.overlap.df.fx, 
	aes(x = target.names, y = result.cnt, fill = cmp.types))
pic.mouse.genes.overlap.fx <- pic.mouse.genes.overlap.fx + 
	geom_col(position = position_fill(reverse = FALSE)) + 
	scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
	#scale_fill_brewer(palette = 1, direction = -1) + 
	scale_fill_manual(values = c("#2172B5", "green")) + 
	xlab("Source") + ylab("Matches") +
	theme_cowplot(12)
dev.new(width = 3.2, height = 6)
ggsave("genes-mouse-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneItself/", device = "pdf")
dev.off()


# --- 2nd ---
# get gene function definition of each database, and to comp those
# import Func from Main-ctrl-process.R

## option 1.1: Cytokine. ALL by definition of InterCellDB
 # NOT RUN

## option 1.2: Growth Factor. ALL by definition of InterCellDB
 # NOT RUN

## option 1.3: Receptor by themselves' definition
	# part 1: percent or stack, remapping against InterCellDB
	tmp.mouse.receptor.db.all <- list(tmp.mouse.receptor.cellchat.res, tmp.mouse.receptor.intercelldb.res)
	names(tmp.mouse.receptor.db.all) <- kDB.namelist[5:6]
	tmp.mouse.receptor.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.mouse.receptor.db.all[1], tmp.mouse.receptor.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 3.2, height = 6)
	ggsave("receptor-mouse-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc-Mouse/", device = "pdf")
	dev.off()
	tmp.mouse.receptor.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.mouse.receptor.db.all[1], tmp.mouse.receptor.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 3.2, height = 6)
	ggsave("receptor-mouse-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc-Mouse/", device = "pdf")
	dev.off()
	# part 2: as overlap is ok, overlap is still needed to show
	tmp.mouse.receptor.notin.ref.all <- list(tmp.mouse.receptor.intercelldb.res, tmp.mouse.receptor.cellchat.res)
	names(tmp.mouse.receptor.notin.ref.all) <- kDB.namelist[6:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc-Mouse/receptor-mouse-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.mouse.receptor.notin.ref.all), fills = kDB.favor.color[6:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
# - done

## option 1.4: Ligand by themselves' definition
# part 1: percent or stack, remapping against InterCellDB
	tmp.mouse.ligand.db.all <- list(tmp.mouse.ligand.cellchat.res, tmp.mouse.ligand.intercelldb.res)
	names(tmp.mouse.ligand.db.all) <- kDB.namelist[5:6]
	tmp.mouse.ligand.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.mouse.ligand.db.all[1], tmp.mouse.ligand.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 3.2, height = 6)
	ggsave("ligand-mouse-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc-Mouse/", device = "pdf")
	dev.off()
	tmp.mouse.ligand.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.mouse.ligand.db.all[1], tmp.mouse.ligand.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 3.2, height = 6)
	ggsave("ligand-mouse-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc-Mouse/", device = "pdf")
	dev.off()
	# part 2: as overlap is not good, overlap is still needed to show
	tmp.mouse.ligand.notin.ref.all <- list(tmp.mouse.ligand.intercelldb.res, tmp.mouse.ligand.cellchat.res)
	names(tmp.mouse.ligand.notin.ref.all) <- kDB.namelist[6:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc-Mouse/ligand-mouse-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.mouse.ligand.notin.ref.all), fills = kDB.favor.color[6:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
# - done

## option 1.5: Integrin. ALL by definition of InterCellDB
 # NOT RUN

## option 1.6: G-protein Coupled Receptor. ALL by definition of InterCellDB
 # NOT RUN



##### For Describing PAIRS - Mouse
# ------
# ------
cmp.mouse.pairs.db.list <- list(CellChatDB = remap.cellchat.pairs.mouse)

### section 1 - collect all pairs DB and compare between each other

## way 1 - percentage matching plot
# for showing matches in 3 independent databases
tmp.mouse.gp.match.colnames <- c("match.in.exp", "match.in.know", "match.in.pred", "unmatch.all")
tmp.mouse.gp.ref.pairs <- list(ref.intercelldb.pairs.mouse.experiments, ref.intercelldb.pairs.mouse.knowledge, ref.intercelldb.pairs.mouse.prediction)
#
src.1d2.mouse.names <- c("CellChatDB")

tmp.mix2.mouse.list <- list()
for (i in seq_along(cmp.mouse.pairs.db.list)) {
	tmp.mix2.mouse.list <- c(tmp.mix2.mouse.list, 
		list(Ctrl.percent.plot.add.data(cmp.mouse.pairs.db.list[[i]], names(cmp.mouse.pairs.db.list)[i], 
			col.desc.vals = tmp.mouse.gp.match.colnames, ref.vector = tmp.mouse.gp.ref.pairs)))
}
tmp.mix2.mouse.df <- dplyr::bind_rows(tmp.mix2.mouse.list)
colnames(tmp.mix2.mouse.df) <- c("source", "match.types", "match.cnt")

# --- plot ---
	# fill pattern
	pic.1d2.mouse.prim <- ggplot(data = tmp.mix2.mouse.df, 
		aes(x = source, y = match.cnt, fill = match.types))
	pic.1d2.mouse.prim <- pic.1d2.mouse.prim + 
		geom_col(position = position_fill(reverse = FALSE)) + 
		scale_x_discrete(limits = src.1d2.mouse.names, breaks = src.1d2.mouse.names) + 
		scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
		scale_fill_brewer(palette = 1, direction = -1) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 3.2, height = 6)
	ggsave("gpairs-mouse-overlap-percent.pdf", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()

	# stack pattern
	pic.1d2.mouse.prim <- ggplot(data = tmp.mix2.mouse.df, 
		aes(x = source, y = match.cnt, fill = match.types))
	pic.1d2.mouse.prim <- pic.1d2.mouse.prim + 
		geom_col(position = position_stack(reverse = FALSE)) +  
		scale_x_discrete(limits = src.1d2.mouse.names, breaks = src.1d2.mouse.names) + 
		scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
		scale_fill_brewer(palette = 1, direction = -1) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 3.2, height = 6)
	ggsave("gpairs-mouse-overlap-stack.pdf", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()


## way 2 - use euler plot to show the overlap of all other DBs with all-DB & 3 sub-DBs of InterCellDB
# [NOTE] as other 5 DBs are mostly experimentally validated gene pairs,
# 		 so here, only compare the sub-DB[exp] of InterCellDB to them

# part w2.0: all-DB
	tmp.gpairs.overlap.all.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse))
	names(tmp.gpairs.overlap.all.list) <- kDB.namelist[5:6]
	# sub1: check overlap
	# use euler
	pdf(file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/gpairs-all-mouse-check-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.gpairs.overlap.all.list), fills = kDB.favor.color[5:6], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		#labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.gpairs.overlap.all.list, 
		filename = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/gpairs-all-mouse-check-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1, euler.d = FALSE, scaled = FALSE)


## Receptor
	tmp.receptor.for.gpairs <- list(tmp.mouse.receptor.cellchat.res, tmp.mouse.receptor.intercelldb.res)
	tmp.receptor.all.shared <- unique(Reduce(intersect, tmp.receptor.for.gpairs))

	# sub 1: use InterCellDB mouse all
	tmp.gpairs.cmp.receptor.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.receptor.run.id <- 1

	# sub 2: use InterCellDB mouse Experiments
	tmp.gpairs.cmp.receptor.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.experiments))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.receptor.run.id <- 2

	# sub 3: use InterCellDB mouse Knowledge
	tmp.gpairs.cmp.receptor.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.knowledge))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.receptor.run.id <- 3

	# sub 4: use InterCellDB mouse Prediction
	tmp.gpairs.cmp.receptor.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.prediction))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.receptor.run.id <- 4

# --- Shared Code ---
# !!!
# [!NOTE!] !!! after run db.mouse.all, sub-DBs can be extracted from it by some math calculation
# !!!
# fetching counts
	tmp.mouse.gpairs.cmp.receptor.genes.list <- lapply(tmp.gpairs.cmp.receptor.list, gp.cut = kgp.cut, 
		function(x, gp.cut) {
			Ctrl.FastOptGenePairs.GetAllRelatedGenes(x, gp.cut)
			})
	tmp.receptor.len.res.df <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.mouse.gpairs.cmp.receptor.genes.list, tmp.receptor.all.shared))
	tmp.receptor.len.res.df$pairs.log.cnt <- log(tmp.receptor.len.res.df$pairs.cnt + 1)

	# to saveRDS
	tmp.receptor.saveRDS.filename <- c("receptor-gpairs-mouse-all-len-df-log.rds",
		"receptor-gpairs-mouse-exps-len-df-log.rds",
		"receptor-gpairs-mouse-know-len-df-log.rds",
		"receptor-gpairs-mouse-pred-len-df-log.rds")
	tmp.receptor.saveRDS.filename <- paste0("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", 
		tmp.receptor.saveRDS.filename)
	saveRDS(tmp.receptor.len.res.df, tmp.receptor.saveRDS.filename[tmp.receptor.run.id])
	
# ggsave file name version 1
	tmp.receptor.ggsave.filename <- c("receptor-gpairs-mouse-all-takein.pdf",
		"receptor-gpairs-mouse-exps-takein.pdf",
		"receptor-gpairs-mouse-know-takein.pdf",
		"receptor-gpairs-mouse-pred-takein.pdf")
# template ggsave code
	dev.new(width = 32, height = 8)
	ggsave(tmp.receptor.ggsave.filename[tmp.receptor.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()

# ggsave file name version 2
	tmp.receptor.ggsave.smooth.filename <- c("receptor-gpairs-mouse-all-takein-smooth.pdf",
		"receptor-gpairs-mouse-exps-takein-smooth.pdf",
		"receptor-gpairs-mouse-know-takein-smooth.pdf",
		"receptor-gpairs-mouse-pred-takein-smooth.pdf")
# template ggsave code
	dev.new(width = 10, height = 8)
	ggsave(tmp.receptor.ggsave.smooth.filename[tmp.receptor.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()

# plot version 1: line, rank by count of InterCellDB receptor gene pairs
	tmp.receptor.getrank.only.intercelldb.v1 <- tmp.receptor.len.res.df[which(tmp.receptor.len.res.df$DBname == "InterCellDB"), ]
	tmp.receptor.getrank.only.intercelldb.v1 <- tmp.receptor.getrank.only.intercelldb.v1[order(tmp.receptor.getrank.only.intercelldb.v1$pairs.log.cnt, decreasing = TRUE), ]
	tmp.top10.receptor.intercelldb.df <- tmp.receptor.getrank.only.intercelldb.v1[1:10, ]
	# change factor
	tmp.receptor.plot.df <- tmp.receptor.len.res.df
	tmp.receptor.plot.df$Gene <- factor(tmp.receptor.plot.df$Gene, levels = unique(tmp.receptor.getrank.only.intercelldb.v1$Gene))
	rm(tmp.receptor.getrank.only.intercelldb.v1)
	tmp.plot <- ggplot(data = tmp.receptor.plot.df) + 
		geom_line(aes(x = Gene, y = pairs.log.cnt, group = DBname, colour = DBname)) + 
		geom_text(data = tmp.top10.receptor.intercelldb.df, aes(x = Gene, y = pairs.log.cnt, label = Gene),
			nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		scale_x_discrete(expand = expansion(add = c(1.2, 1))) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# - done

# plot version 2: smmoth and ribbon, rank by count of InterCellDB receptor gene pairs
	tmp.receptor.getrank.only.intercelldb.v2 <- tmp.receptor.len.res.df[which(tmp.receptor.len.res.df$DBname == "InterCellDB"), ]
	tmp.receptor.getrank.only.intercelldb.v2 <- tmp.receptor.getrank.only.intercelldb.v2[order(tmp.receptor.getrank.only.intercelldb.v2$pairs.log.cnt, decreasing = TRUE), ]
	# fetch
	tmp.receptor.plot.df <- tmp.receptor.len.res.df
	tmp.receptor.intercelldb.gene.rank.vec <- unique(tmp.receptor.getrank.only.intercelldb.v2$Gene)
	rm(tmp.receptor.getrank.only.intercelldb.v2)
	# change factor
	tmp.receptor.plot.df$Gene <- factor(tmp.receptor.plot.df$Gene, levels = tmp.receptor.intercelldb.gene.rank.vec)
	tmp.receptor.num.gene <- data.frame(num.gene = seq_along(tmp.receptor.intercelldb.gene.rank.vec), gene = tmp.receptor.intercelldb.gene.rank.vec, stringsAsFactors = FALSE)
	tmp.receptor.plot.df <- left_join(tmp.receptor.plot.df, tmp.receptor.num.gene, 
		by = c("Gene" = "gene"))
	tmp.top10.receptor.plot.df <- tmp.receptor.plot.df[which(tmp.receptor.plot.df$num.gene %in% 1:10), ]
	tmp.top10.receptor.plot.df <- tmp.top10.receptor.plot.df[order(tmp.top10.receptor.plot.df$num.gene), ]
	tmp.receptor.ribbon.colour <- tmp.receptor.ribbon.fill <- kDB.favor.color
	names(tmp.receptor.ribbon.colour) <- names(tmp.receptor.ribbon.fill) <- kDB.namelist
	tmp.plot <- ggplot(data = tmp.receptor.plot.df) + 
		geom_point(aes(x = num.gene, y = pairs.log.cnt, group = DBname, colour = DBname), size = 1) + 
		geom_smooth(aes(x = num.gene, y = pairs.log.cnt, group = DBname, colour = DBname, fill = DBname),
			method = "loess", formula = y ~ x) + 
		scale_colour_manual(values = tmp.receptor.ribbon.colour) + 
		scale_fill_manual(values = tmp.receptor.ribbon.fill, guide = guide_legend(override.aes = list(alpha = 0.4))) + 
		geom_text(data = tmp.top10.receptor.plot.df[which(tmp.top10.receptor.plot.df$DBname == "InterCellDB"), ], aes(x = num.gene, y = pairs.log.cnt, label = Gene),
			nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		theme_cowplot(12)
		#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# - done


## Ligand
	tmp.ligand.for.gpairs <- list(tmp.mouse.ligand.cellchat.res, tmp.mouse.ligand.intercelldb.res)
	tmp.ligand.all.shared <- unique(Reduce(intersect, tmp.ligand.for.gpairs))

	# sub 1: use InterCellDB mouse all
	tmp.gpairs.cmp.ligand.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.ligand.run.id <- 1

	# sub 2: use InterCellDB mouse Experiments
	tmp.gpairs.cmp.ligand.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.experiments))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.ligand.run.id <- 2

	# sub 3: use InterCellDB mouse Knowledge
	tmp.gpairs.cmp.ligand.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.knowledge))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.ligand.run.id <- 3

	# sub 4: use InterCellDB mouse Prediction
	tmp.gpairs.cmp.ligand.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse.prediction))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[5], "InterCellDB")
	tmp.ligand.run.id <- 4

# --- Shared Code ---
# !!!
# [!NOTE!] !!! after run db.mouse.all, sub-DBs can be extracted from it by some math calculation
# !!!
# fetching counts
	tmp.mouse.gpairs.cmp.ligand.genes.list <- lapply(tmp.gpairs.cmp.ligand.list, gp.cut = kgp.cut, 
		function(x, gp.cut) {
			Ctrl.FastOptGenePairs.GetAllRelatedGenes(x, gp.cut)
			})
	tmp.ligand.len.res.df <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.mouse.gpairs.cmp.ligand.genes.list, tmp.ligand.all.shared))
	tmp.ligand.len.res.df$pairs.log.cnt <- log(tmp.ligand.len.res.df$pairs.cnt + 1)

	# to saveRDS
	tmp.ligand.saveRDS.filename <- c("ligand-gpairs-mouse-all-len-df-log.rds",
		"ligand-gpairs-mouse-exps-len-df-log.rds",
		"ligand-gpairs-mouse-know-len-df-log.rds",
		"ligand-gpairs-mouse-pred-len-df-log.rds")
	tmp.ligand.saveRDS.filename <- paste0("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", 
		tmp.ligand.saveRDS.filename)
	saveRDS(tmp.ligand.len.res.df, tmp.ligand.saveRDS.filename[tmp.ligand.run.id])
	
# ggsave file name version 1
	tmp.ligand.ggsave.filename <- c("ligand-gpairs-mouse-all-takein.pdf",
		"ligand-gpairs-mouse-exps-takein.pdf",
		"ligand-gpairs-mouse-know-takein.pdf",
		"ligand-gpairs-mouse-pred-takein.pdf")
# template ggsave code
	dev.new(width = 32, height = 8)
	ggsave(tmp.ligand.ggsave.filename[tmp.ligand.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()

# ggsave file name version 2
	tmp.ligand.ggsave.smooth.filename <- c("ligand-gpairs-mouse-all-takein-smooth.pdf",
		"ligand-gpairs-mouse-exps-takein-smooth.pdf",
		"ligand-gpairs-mouse-know-takein-smooth.pdf",
		"ligand-gpairs-mouse-pred-takein-smooth.pdf")
# template ggsave code
	dev.new(width = 10, height = 8)
	ggsave(tmp.ligand.ggsave.smooth.filename[tmp.ligand.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "pdf")
	dev.off()

# plot version 1: line, rank by count of InterCellDB ligand gene pairs
	tmp.ligand.getrank.only.intercelldb.v1 <- tmp.ligand.len.res.df[which(tmp.ligand.len.res.df$DBname == "InterCellDB"), ]
	tmp.ligand.getrank.only.intercelldb.v1 <- tmp.ligand.getrank.only.intercelldb.v1[order(tmp.ligand.getrank.only.intercelldb.v1$pairs.log.cnt, decreasing = TRUE), ]
	tmp.top10.ligand.intercelldb.df <- tmp.ligand.getrank.only.intercelldb.v1[1:10, ]
	# change factor
	tmp.ligand.plot.df <- tmp.ligand.len.res.df
	tmp.ligand.plot.df$Gene <- factor(tmp.ligand.plot.df$Gene, levels = unique(tmp.ligand.getrank.only.intercelldb.v1$Gene))
	rm(tmp.ligand.getrank.only.intercelldb.v1)
	tmp.plot <- ggplot(data = tmp.ligand.plot.df) + 
		geom_line(aes(x = Gene, y = pairs.log.cnt, group = DBname, colour = DBname)) + 
		geom_text(data = tmp.top10.ligand.intercelldb.df, aes(x = Gene, y = pairs.log.cnt, label = Gene),
			nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		scale_x_discrete(expand = expansion(add = c(1.2, 1))) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# - done

# plot version 2: smmoth and ribbon, rank by count of InterCellDB ligand gene pairs
	tmp.ligand.getrank.only.intercelldb.v2 <- tmp.ligand.len.res.df[which(tmp.ligand.len.res.df$DBname == "InterCellDB"), ]
	tmp.ligand.getrank.only.intercelldb.v2 <- tmp.ligand.getrank.only.intercelldb.v2[order(tmp.ligand.getrank.only.intercelldb.v2$pairs.log.cnt, decreasing = TRUE), ]
	# fetch
	tmp.ligand.plot.df <- tmp.ligand.len.res.df
	tmp.ligand.intercelldb.gene.rank.vec <- unique(tmp.ligand.getrank.only.intercelldb.v2$Gene)
	rm(tmp.ligand.getrank.only.intercelldb.v2)
	# change factor
	tmp.ligand.plot.df$Gene <- factor(tmp.ligand.plot.df$Gene, levels = tmp.ligand.intercelldb.gene.rank.vec)
	tmp.ligand.num.gene <- data.frame(num.gene = seq_along(tmp.ligand.intercelldb.gene.rank.vec), gene = tmp.ligand.intercelldb.gene.rank.vec, stringsAsFactors = FALSE)
	tmp.ligand.plot.df <- left_join(tmp.ligand.plot.df, tmp.ligand.num.gene, 
		by = c("Gene" = "gene"))
	tmp.top10.ligand.plot.df <- tmp.ligand.plot.df[which(tmp.ligand.plot.df$num.gene %in% 1:10), ]
	tmp.top10.ligand.plot.df <- tmp.top10.ligand.plot.df[order(tmp.top10.ligand.plot.df$num.gene), ]
	tmp.ligand.ribbon.colour <- tmp.ligand.ribbon.fill <- kDB.favor.color
	names(tmp.ligand.ribbon.colour) <- names(tmp.ligand.ribbon.fill) <- kDB.namelist
	tmp.plot <- ggplot(data = tmp.ligand.plot.df) + 
		geom_point(aes(x = num.gene, y = pairs.log.cnt, group = DBname, colour = DBname), size = 1) + 
		geom_smooth(aes(x = num.gene, y = pairs.log.cnt, group = DBname, colour = DBname, fill = DBname),
			method = "loess", formula = y ~ x) + 
		scale_colour_manual(values = tmp.ligand.ribbon.colour) + 
		scale_fill_manual(values = tmp.ligand.ribbon.fill, guide = guide_legend(override.aes = list(alpha = 0.4))) + 
		geom_text(data = tmp.top10.ligand.plot.df[which(tmp.top10.ligand.plot.df$DBname == "InterCellDB"), ], aes(x = num.gene, y = pairs.log.cnt, label = Gene),
			nudge_y = seq(from=2,to=1.1,by=-0.1), size = 2) + 
		theme_cowplot(12)
		#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# - done





#### GO terms related genes
# target DB list
tmp.mouse.go.db.list <- c(cmp.mouse.pairs.db.list, list(ref.intercelldb.pairs.mouse, 
	ref.intercelldb.pairs.mouse.experiments, ref.intercelldb.pairs.mouse.knowledge, ref.intercelldb.pairs.mouse.prediction))
names(tmp.mouse.go.db.list) <- c(kDB.namelist[5], "InterCellDB.mouse.all", "InterCellDB.mouse.exps", "InterCellDB.mouse.know", "InterCellDB.mouse.pred")

tmp.mouse.go.db.related.genes.list <- lapply(tmp.mouse.go.db.list, gp.cut = kgp.cut, 
	function(x, gp.cut) {
		Ctrl.FastOptGenePairs.GetAllRelatedGenes(x, gp.cut)
		})

### version 1: plot heatmap by arrange genes along x-axis
	## GO:0006955 [immune response]
	tmp.mouse.go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.immune.response.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.immune.response.genes))
	tmp.mouse.go.immune.response.id.ref.df <- data.frame(gene = tmp.mouse.go.immune.response.genes[order(tmp.mouse.go.immune.response.genes)], 
		num.id = seq_along(tmp.mouse.go.immune.response.genes),
		stringsAsFactors = FALSE)
	tmp.mouse.go.immune.response.counts <- dplyr::left_join(tmp.mouse.go.immune.response.counts, tmp.mouse.go.immune.response.id.ref.df, by = c("Gene" = "gene"))

	## GO:0001525 [angiogenesis]
	tmp.mouse.go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.vessel.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.vessel.genes))
	tmp.mouse.go.vessel.id.ref.df <- data.frame(gene = tmp.mouse.go.vessel.genes[order(tmp.mouse.go.vessel.genes)], 
		num.id = seq_along(tmp.mouse.go.vessel.genes),  # add on
		stringsAsFactors = FALSE)
	# add on num.id
	tmp.mouse.go.vessel.id.ref.df$num.id <- tmp.mouse.go.vessel.id.ref.df$num.id + length(tmp.mouse.go.immune.response.genes)
	tmp.mouse.go.vessel.counts <- dplyr::left_join(tmp.mouse.go.vessel.counts, tmp.mouse.go.vessel.id.ref.df, by = c("Gene" = "gene"))

	## GO:0001837 [epithelial to mesenchymal transition]
	tmp.mouse.go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.emt.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.emt.genes))
	tmp.mouse.go.emt.id.ref.df <- data.frame(gene = tmp.mouse.go.emt.genes[order(tmp.mouse.go.emt.genes)], 
		num.id = seq_along(tmp.mouse.go.emt.genes),  # add on
		stringsAsFactors = FALSE)
	# add on num.id
	tmp.mouse.go.emt.id.ref.df$num.id <- tmp.mouse.go.emt.id.ref.df$num.id + length(tmp.mouse.go.immune.response.genes) + length(tmp.mouse.go.vessel.genes)
	tmp.mouse.go.emt.counts <- dplyr::left_join(tmp.mouse.go.emt.counts, tmp.mouse.go.emt.id.ref.df, by = c("Gene" = "gene"))


	## merge counts
	tmp.mouse.go.counts.merge <- Reduce(rbind, list(tmp.mouse.go.immune.response.counts, tmp.mouse.go.vessel.counts, tmp.mouse.go.emt.counts))
	# log transform
	tmp.mouse.go.counts.merge$pairs.cnt <- log(tmp.mouse.go.counts.merge$pairs.cnt + 1)
	tmp.mouse.go.counts.merge$DBname <- factor(tmp.mouse.go.counts.merge$DBname, 
		levels = names(tmp.mouse.go.db.list))

	tmp.mouse.go.len.1 <- length(tmp.mouse.go.immune.response.genes)
	tmp.mouse.go.len.2 <- length(tmp.mouse.go.vessel.genes)
	tmp.mouse.go.len.3 <- length(tmp.mouse.go.emt.genes)
	tmp.mouse.go.counts.breaks.IT <- c(1, tmp.mouse.go.len.1 / 2, tmp.mouse.go.len.1,  
		tmp.mouse.go.len.1 + tmp.mouse.go.len.2 / 2, tmp.mouse.go.len.1 + tmp.mouse.go.len.2,  
		tmp.mouse.go.len.1 + tmp.mouse.go.len.2 + tmp.mouse.go.len.3 / 2, tmp.mouse.go.len.1 + tmp.mouse.go.len.2 + tmp.mouse.go.len.3)
	tmp.mouse.go.counts.breaks.labels <- tmp.mouse.go.counts.breaks.IT
	tmp.mouse.go.counts.breaks.labels[c(2,4,6)] <- c("Immune Response", "Angiogenesis", "EMT")

	tmp.mouse.plot <- ggplot(data = tmp.mouse.go.counts.merge) + 
		geom_raster(aes(x = num.id, y = DBname, fill = pairs.cnt)) + 
		scale_x_continuous(breaks = tmp.mouse.go.counts.breaks.IT, 
			labels = tmp.mouse.go.counts.breaks.labels) + 
		theme_cowplot(12) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	dev.new(width = 12, height = 8)
	ggsave("go-mix3-gpairs-mouse.png", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "png")
	dev.off()
### - version 1 plot GO end -


### version 2: plot heatmap by arrange GO terms along x-axis
	## GO:0006955 [immune response]
	tmp.mouse.go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.immune.response.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.immune.response.genes))
	tmp.mouse.go.immune.response.mergecnt$GOterm <- "immune response"

	## GO:0001525 [angiogenesis]
	tmp.mouse.go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.vessel.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.vessel.genes))
	tmp.mouse.go.vessel.mergecnt$GOterm <- "angiogenesis"

	## GO:0001837 [epithelial to mesenchymal transition]
	tmp.mouse.go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.emt.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.emt.genes))
	tmp.mouse.go.emt.mergecnt$GOterm <- "EMT"

	## GO:0006956 [complement activation]
	tmp.mouse.go.comp.act.genes <- Tool.FindGenesFromGO("GO:0006956", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.comp.act.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.comp.act.genes))
	tmp.mouse.go.comp.act.mergecnt$GOterm <- "complement activation"

	## GO:0005125 [cytokine activity]
	tmp.mouse.go.cytokine.act.genes <- Tool.FindGenesFromGO("GO:0005125", genes.mouse.ref.db, go.mouse.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.mouse.go.cytokine.act.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.mouse.go.db.related.genes.list, tmp.mouse.go.cytokine.act.genes))
	tmp.mouse.go.cytokine.act.mergecnt$GOterm <- "cytokine activity"


	# merge result
	tmp.mouse.go.mergecnt.collect <- Reduce(rbind, list(tmp.mouse.go.immune.response.mergecnt,
		tmp.mouse.go.vessel.mergecnt, tmp.mouse.go.emt.mergecnt, tmp.mouse.go.comp.act.mergecnt, tmp.mouse.go.cytokine.act.mergecnt))
	tmp.mouse.go.mergecnt.collect$log.cnt <- log(tmp.mouse.go.mergecnt.collect$merge.pairs.cnt)
	tmp.mouse.go.mergecnt.collect$DBname <- factor(tmp.mouse.go.mergecnt.collect$DBname, levels = names(tmp.mouse.go.db.list))

	tmp.plot <- ggplot(data = tmp.mouse.go.mergecnt.collect) + 
		geom_raster(aes(x = GOterm, y = DBname, fill = log.cnt)) + 
		scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	dev.new(width = 8, height = 8)
	ggsave("go-mergecnt-mix3-gpairs-mouse.png", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB-Mouse/", device = "png")
	dev.off()
### - version 2 plot GO end -
