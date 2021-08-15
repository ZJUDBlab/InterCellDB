#
# This file records the workflow of steps (all human and some mouse steps move to other files)

# new library
library(eulerr)

# global variables
kDB.namelist <- c("CellPhoneDB", "NicheNet", "iTALK", "SingleCellSignalR", "CellChatDB", "InterCellDB")
kDB.favor.color <- c("#CEC3F2", "#96CAEA", "#EF98CE", "#D9C5A4", "#8CACAF", "#BCC3C4")

# "#CEC3F2" light purple 	[sel: CellPhoneDB]
# "#96CAEA" light blue 		[sel: NicheNet]
# "#8CACAF" light green 	[sel: CellChatDB]
# "#E3EDA0" yellow green 	[sel: ]
# "#EF98CE" young red 		[sel: iTALK]
# "#D9C5A4" soft brown 		[sel: SingleCellSignalR]
# "#BCC3C4" grey green 		[sel: InterCellDB]
# "#D94B4B" 2orange red 	[sel: merge5-genes]

# tool functions for main ctrl or accessory files
#
kgp.cut <- ">"
#
Ctrl.RemapGenePairs <- function(
	cmp.db,
	genes.ref.db,
	cmp.def.pairs.colnames  # sequential to be A~B
) {
	inter.A <- cmp.db[, cmp.def.pairs.colnames[1], drop = FALSE]
	inter.A <- cbind(inter.A, row.id = 1:nrow(cmp.db))
	inter.B <- cmp.db[, cmp.def.pairs.colnames[2], drop	= FALSE]
	inter.B <- cbind(inter.B, row.id = 1:nrow(cmp.db))

	# remapping genes
	new.inter.A <- Tool.AddUserRestrictDB(inter.A, genes.ref.db)
	new.inter.A <- new.inter.A[order(new.inter.A$row.id), ]
	colnames(new.inter.A)[1:2] <- c("align.GeneID.A", "align.GeneName.A")
	#
	new.inter.B <- Tool.AddUserRestrictDB(inter.B, genes.ref.db)
	new.inter.B <- new.inter.B[order(new.inter.B$row.id), ]
	colnames(new.inter.B)[1:2] <- c("align.GeneID.B", "align.GeneName.B")

	# 
	new.inter.it <- dplyr::left_join(new.inter.A, new.inter.B, by = c("row.id" = "row.id"))
	new.inter.it$row.id <- NULL
	#
	new.inter.it <- new.inter.it[, c(1,3,2,4)]
	# align database
	new.inter.it <- DoPartUnique(FastAlignPairs(new.inter.it, 4), 3:4)
	#
	return(new.inter.it)
}
#

# This function is default served for comparing gene pairs
# but currently extended to be use for check gene overlaps
Ctrl.percent.plot.add.data <- function(
	cmp.std.onevec,  # the values(in vector()) of the to-compared target.
	cmp.anno.name,  # the name for the to-compared target
	col.desc.vals,  # describe the types of comparing, which is corresponding to `ref.vector` 1-by-1. If show unmatches, 1 extra name will be added 
	ref.vector = list(),  
	count.or.identity = TRUE,  # TRUE use count, FALSE use identity
	show.unmatches = TRUE  # if show unmatches, then all rest ones after matching will be counted
) {
	check.len.diff <- length(col.desc.vals) - length(ref.vector)
	if (show.unmatches == TRUE && check.len.diff != 1) {
		stop("Colnames for unmatches not given or too much!")
	}
	if (show.unmatches == FALSE && check.len.diff != 0) {
		stop("Not show unmatches, but given the colname of it!")
	}
	if (count.or.identity == TRUE) {
		tmp.res <- sapply(ref.vector, USE.NAMES = FALSE, cmp.std.onevec = cmp.std.onevec, 
			function(x, cmp.std.onevec) {
				length(intersect(x, cmp.std.onevec))
			})
		if (show.unmatches == TRUE) {
			tmp.res <- c(tmp.res, length(cmp.std.onevec) - sum(tmp.res))
		}
		return(data.frame(target.names = cmp.anno.name, 
			cmp.types = col.desc.vals, 
			result.cnt = tmp.res, 
			stringsAsFactors = FALSE))
	} else {
		tmp.res <- lapply(ref.vector, cmp.std.onevec = cmp.std.onevec, 
			function(x, cmp.std.onevec) {
				intersect(cmp.std.onevec, x)
			})
		if (show.unmatches == TRUE) {
			tmp.res <- c(tmp.res, list(setdiff(cmp.std.onevec, unlist(tmp.res))))
		}
		names(tmp.res) <- col.desc.vals
		ret.res <- list(tmp.res)
		names(ret.res) <- cmp.anno.name
		return(ret.res)
	}
}
# include the plotting function inside to optimize the usage of function
Ctrl.percent.plot.2item <- function(
	tmp.some.db.all,  # the all 
	tmp.some.intercelldb.res,
	stack.or.fill = TRUE,  # stack is just stacking, fill is stacking and getting all same height
	scale.y.add = scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1)))
				# scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
) {
	# some(like receptor, ligand, etc) percent matching to InterCellDB
	tmp.some.overlap.list <- list()
	for (i in seq_along(tmp.some.db.all)) {
		tmp.some.overlap.list <- c(tmp.some.overlap.list, list(Ctrl.percent.plot.add.data(tmp.some.db.all[[i]], names(tmp.some.db.all)[i],
			col.desc.vals = c("matched.genes", "unmatched.genes"), 
			ref.vector = list(tmp.some.intercelldb.res),
			count.or.identity = TRUE, 
			show.unmatches = TRUE)))
	}
	tmp.some.overlap.df <- dplyr::bind_rows(tmp.some.overlap.list)
	tmp.some.overlap.df$target.names <- factor(tmp.some.overlap.df$target.names, levels = names(tmp.some.db.all))
	# --- plot ---
	pic.some.overlap <- ggplot(data = tmp.some.overlap.df, 
		aes(x = target.names, y = result.cnt, fill = cmp.types))
	# add item stack or fill
	if (stack.or.fill == TRUE) {
		pic.some.add.1 <- geom_col(position = position_stack(reverse = FALSE))
	} else {
		pic.some.add.1 <- geom_col(position = position_fill(reverse = FALSE))

	}
	pic.some.overlap <- pic.some.overlap + 
		pic.some.add.1 + 
		scale.y.add + 
		#scale_fill_brewer(palette = 1, direction = -1) + 
		scale_fill_manual(values = c("#2172B5", "green")) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	return(pic.some.overlap)
}

# 
# This function get count of gene pairs of all included genes
Ctrl.FastOptGenePairs.GetAllRelatedGenes <- function(
	used.cat.pairs.list,  # to explore related gene pairs in which database. The format of gene pairs is A>B, ">" is defined by kgp.cut
	gp.cut  # the char used to gap the 2 genes in one gene pair
) {
	tmp.splits <- as.character(unlist(strsplit(used.cat.pairs.list, split = gp.cut, fixed = TRUE)))
	tmp.seqs <- seq_len(length(tmp.splits) / 2)
	tmp.the.former.split <- tmp.splits[2 * tmp.seqs - 1]
	tmp.the.latter.split <- tmp.splits[2 * tmp.seqs]

	# check those A>A gene pairs and subtract those
	tmp.inds.dup <- which(tmp.the.former.split == tmp.the.latter.split)
	tmp.dup.genes <- tapply(seq_along(tmp.inds.dup), tmp.the.former.split[tmp.inds.dup], length)
	if (length(which(tmp.dup.genes > 1)) != 0) {
		stop("Error: dup genes larget than 1. Detect database not unique.")
	}

	# each gene count
	tmp.proc.res <- tapply(seq_along(tmp.splits), tmp.splits, length)

	# to substract dup
	tmp.inds.to.substact <- match(names(tmp.dup.genes), names(tmp.proc.res))
	tmp.proc.res[tmp.inds.to.substact] <- tmp.proc.res[tmp.inds.to.substact] - tmp.dup.genes
	return(tmp.proc.res)
}

# get count of participating gene pairs of a set of genes
Ctrl.FastOptGenePairs.MatchEachGene <- function(
	db.related.genes.list,
	go.some.genes
) {
	this.find.related.res <- lapply(seq_along(db.related.genes.list), 
		split.db.genes.list = db.related.genes.list, go.some.genes = go.some.genes, 
		function(x, split.db.genes.list, go.some.genes) {
			tmp.IT <- split.db.genes.list[[x]]
			tmp.inds <- match(go.some.genes, names(tmp.IT))
			names(tmp.IT) <- NULL
			tmp.df <- data.frame(DBname = names(split.db.genes.list)[x], 
				Gene = go.some.genes, 
				pairs.cnt = tmp.IT[tmp.inds],
				stringsAsFactors = FALSE)
			tmp.df[which(is.na(tmp.df$pairs.cnt)), "pairs.cnt"] <- 0
			tmp.df
			})
	return(this.find.related.res)
}

# get count of participating gene pairs of a set of genes
Ctrl.FastOptGenePairs.MergeGORelCnts <- function(
	db.related.genes.list,
	go.some.genes
) {
	this.find.related.res <- lapply(seq_along(db.related.genes.list), 
		split.db.genes.list = db.related.genes.list, go.some.genes = go.some.genes, 
		function(x, split.db.genes.list, go.some.genes) {
			tmp.IT <- split.db.genes.list[[x]]
			tmp.inds <- match(go.some.genes, names(tmp.IT))
			names(tmp.IT) <- NULL
			tmp.df <- data.frame(DBname = names(split.db.genes.list)[x], 
				merge.pairs.cnt = sum(tmp.IT[tmp.inds[which(!is.na(tmp.inds))]]),
				stringsAsFactors = FALSE)
			tmp.df
			})
	return(this.find.related.res)
}

# get overlap within items of one list
Ctrl.check.inter.overlap <- function(
	cmp.list  # names must be included
) {
	this.res.list <- list()
	for (i in seq_along(cmp.list)) {  # goes from y, pick one DB-Y each time
		this.dummy.name <- names(cmp.list)[i]
		for (j in seq_along(cmp.list)) {
			tmp.cmp.name <- names(cmp.list)[j]
			tmp.inter.cnt <- length(intersect(cmp.list[[i]], cmp.list[[j]]))
			this.res.list <- c(this.res.list, list(data.frame(DBx = tmp.cmp.name, DBy = this.dummy.name, 
				inter.cnt = tmp.inter.cnt,
				inter.percent = tmp.inter.cnt / length(cmp.list[[i]]),
				stringsAsFactors = FALSE)))
		}
	}
	bind_rows(this.res.list)
} 



##### For Describing GENE - Human
# --- 1st --- 
#   Extract plain genes for every DB
# 
# {deprecated at 2021.01.29} [NOTE] if DB has gene reference database, directly use it, or collect genes from gene pairs.
# CellChatDB, InterCellDB have gene ref DB, other 4 not.
# {added at 2021.01.29} [NOTE] count genes only when they take part in specific gene pairs
#

tmp.genes.cnt.all <- list(cellphone.genes.res, nichenet.genes.res, italk.genes.res, 
	scsr.genes.res, cellchat.human.genes.res, intercelldb.human.genes.res)
names(tmp.genes.cnt.all) <- kDB.namelist
# part0: directly show all intersections of all DBs
pdf(file = "../anal-peer-db-Result/GeneItself/genes-human-allDBs-redundancy.pdf", 
	width = 9, height = 6)
plot(euler(tmp.genes.cnt.all), fills = kDB.favor.color, 
	edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
	labels = list(fontfamily = "sans", cex = 0.8), 
	quantities = list(cex = 0.8), 
	legend = TRUE)
dev.off()

# part1.1: merge the 1:5 DB to compare between themselves
#
pdf(file = "../anal-peer-db-Result/GeneItself/genes-human-merge5-redundancy.pdf", 
	width = 9, height = 6)
plot(euler(tmp.genes.cnt.all[1:5]), quantities = TRUE, 
	fills = kDB.favor.color, edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
	labels = list(fontfamily = "sans", cex = 1), legend = TRUE)
dev.off()

# part1.2: compare the merge5 and InterCellDB
tmp.genes.cnt.last2 <- list(unique(as.character(unlist(tmp.genes.cnt.all[1:5]))), intercelldb.human.genes.res)
names(tmp.genes.cnt.last2) <- c("Merge5Genes", "InterCellDB")
pdf(file = "../anal-peer-db-Result/GeneItself/genes-human-last2-redundancy.pdf", 
	width = 9, height = 6)
plot(euler(tmp.genes.cnt.last2), quantities = TRUE, 
	fills = c("#D94B4B", kDB.favor.color[6]), edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
	labels = list(fontfamily = "sans", cex = 1), legend = TRUE)
dev.off()

# part2(form 1): give the genes lacked in InterCellDB
# directly show the genes
tmp.genes.overlap.list.f1 <- list()
for (i in 1:5) {
	tmp.inside.res <- Ctrl.percent.plot.add.data(tmp.genes.cnt.all[[i]], names(tmp.genes.cnt.all)[i],
		col.desc.vals = c("matched.genes", "unmatched.genes"), 
		ref.vector = list(intercelldb.human.genes.res),
		count.or.identity = FALSE, 
		show.unmatches = TRUE)
	tmp.genes.overlap.list.f1 <- c(tmp.genes.overlap.list.f1, 
		list(data.frame(source = names(tmp.genes.cnt.all)[i], 
			unmatched.genes = tmp.inside.res[[1]]$unmatched.genes, 
			stringsAsFactors = FALSE)))
}
tmp.genes.overlap.df.f1 <- dplyr::bind_rows(tmp.genes.overlap.list.f1)
tmp.genes.overlap.df.f1$source <- factor(tmp.genes.overlap.df.f1$source, levels = names(tmp.genes.cnt.all)[1:5])
pic.genes.overlap.f1 <- ggplot(data = tmp.genes.overlap.df.f1, 
	aes(x = source, y = unmatched.genes))
pic.genes.overlap.f1 <- pic.genes.overlap.f1 + 
	geom_point() +  # fill = "#DDEBF7"
	#scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
	#scale_fill_brewer(palette = 1, direction = -1) + 
	xlab("Source") + ylab("Unmatched count") +
	theme_cowplot(12)
dev.new(width = 7.7, height = 5.8)
ggsave("genes-human-lack-cmp-intercelldb-gene-points.pdf", path = "../anal-peer-db-Result/GeneItself", device = "pdf")
dev.off()

# part2(form 2): give the genes lacked in InterCellDB
# only the count of lacked genes
tmp.genes.overlap.list.f2 <- list()
for (i in 1:5) {
	tmp.genes.overlap.list.f2 <- c(tmp.genes.overlap.list.f2, list(Ctrl.percent.plot.add.data(tmp.genes.cnt.all[[i]], names(tmp.genes.cnt.all)[i],
		col.desc.vals = c("matched.genes", "unmatched.genes"), 
		ref.vector = list(intercelldb.human.genes.res),
		count.or.identity = TRUE, 
		show.unmatches = TRUE)))
}
tmp.genes.overlap.df.f2 <- dplyr::bind_rows(tmp.genes.overlap.list.f2)
tmp.genes.overlap.df.f2$target.names <- factor(tmp.genes.overlap.df.f2$target.names, levels = names(tmp.genes.cnt.all)[1:5])
tmp.genes.overlap.df.f2 <- tmp.genes.overlap.df.f2[which(tmp.genes.overlap.df.f2$cmp.types == "unmatched.genes"), ]
# --- plot ---
pic.genes.overlap.f2 <- ggplot(data = tmp.genes.overlap.df.f2, 
	aes(x = target.names, y = result.cnt, fill = cmp.types))
pic.genes.overlap.f2 <- pic.genes.overlap.f2 + 
	geom_col() +  # fill = "#DDEBF7"
	scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
	scale_fill_brewer(palette = 1, direction = -1) + 
	xlab("Source") + ylab("Unmatched count") +
	theme_cowplot(12)
dev.new(width = 7.2, height = 6)
ggsave("genes-human-lack-cmp-intercelldb.pdf", path = "../anal-peer-db-Result/GeneItself/", device = "pdf")
dev.off()

# part x: give the genes overlap percentage by referring InterCellDB
tmp.genes.overlap.df.fx <- dplyr::bind_rows(tmp.genes.overlap.list.f2)
tmp.genes.overlap.df.fx$target.names <- factor(tmp.genes.overlap.df.fx$target.names, levels = names(tmp.genes.cnt.all)[1:5])
# --- plot ---
pic.genes.overlap.fx <- ggplot(data = tmp.genes.overlap.df.fx, 
	aes(x = target.names, y = result.cnt, fill = cmp.types))
pic.genes.overlap.fx <- pic.genes.overlap.fx + 
	geom_col(position = position_fill(reverse = FALSE)) + 
	scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
	#scale_fill_brewer(palette = 1, direction = -1) + 
	scale_fill_manual(values = c("#2172B5", "green")) + 
	xlab("Source") + ylab("Matches") +
	theme_cowplot(12)
dev.new(width = 7.2, height = 6)
ggsave("genes-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneItself/", device = "pdf")
dev.off()


# --- 2nd ---
# (option 1): get gene function definition of each database, and to comp those
#
## option 1.1: Cytokine by themselves' definition, and if not defined, use the definition of InterCellDB
# part 1: percent or stack, remapping against InterCellDB
	tmp.cytokine.db.all <- list(tmp.cytokine.cpdb.res, tmp.cytokine.nichenet.res, tmp.cytokine.italk.res, 
		tmp.cytokine.scsr.res, tmp.cytokine.cellchat.res, tmp.cytokine.intercelldb.res)
	names(tmp.cytokine.db.all) <- kDB.namelist
	tmp.cytokine.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.cytokine.db.all[1:5], tmp.cytokine.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("cytokine-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	tmp.cytokine.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.cytokine.db.all[1:5], tmp.cytokine.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 7.2, height = 6)
	ggsave("cytokine-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# - done
# part 2: overlap is ok, still need to check if InterCellDB lost some cytokine in definition
	tmp.cytokine.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.cytokine.notin.ref.all <- c(tmp.cytokine.notin.ref.all, list(setdiff(tmp.cytokine.db.all[[i]], tmp.cytokine.intercelldb.res)))
	}
	names(tmp.cytokine.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/cytokine-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.cytokine.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn, only 2, so select 2 to show
	venn.diagram(tmp.cytokine.notin.ref.all[c(1,3)],  # CellPhoneDB and iTALK
		filename = "../anal-peer-db-Result/GeneFunc/cytokine-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done


## option 1.2: Growth Factor by themselves' definition
# part 1: percent or stack, remapping against InterCellDB
	tmp.growfactor.db.all <- list(tmp.growfactor.cpdb.res, tmp.growfactor.nichenet.res, tmp.growfactor.italk.res, 
		tmp.growfactor.scsr.res, tmp.growfactor.cellchat.res, tmp.growfactor.intercelldb.res)
	names(tmp.growfactor.db.all) <- kDB.namelist
	tmp.growfactor.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.growfactor.db.all[1:5], tmp.growfactor.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("growfactor-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	tmp.growfactor.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.growfactor.db.all[1:5], tmp.growfactor.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 7.2, height = 6)
	ggsave("growfactor-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# - done
# part 2: overlap is ok, still need to check if InterCellDB lost some growth factor in definition
	tmp.growfactor.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.growfactor.notin.ref.all <- c(tmp.growfactor.notin.ref.all, list(setdiff(tmp.growfactor.db.all[[i]], tmp.growfactor.intercelldb.res)))
	}
	names(tmp.growfactor.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/growfactor-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.growfactor.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn, only 2, so select 2 to show
	venn.diagram(tmp.growfactor.notin.ref.all[c(1,3)],  # CellPhoneDB and iTALK
		filename = "../anal-peer-db-Result/GeneFunc/growfactor-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done


## option 1.3: Receptor by themselves' definition
# part 1: percent or stack, remapping against InterCellDB
	tmp.receptor.db.all <- list(tmp.receptor.cpdb.res, tmp.receptor.nichenet.res, tmp.receptor.italk.res, 
		tmp.receptor.scsr.res, tmp.receptor.cellchat.res, tmp.receptor.intercelldb.res)
	names(tmp.receptor.db.all) <- kDB.namelist
	# receptor percent matching to InterCellDB
	tmp.receptor.overlap.list <- list()
	for (i in 1:5) {
		tmp.receptor.overlap.list <- c(tmp.receptor.overlap.list, list(Ctrl.percent.plot.add.data(tmp.receptor.db.all[[i]], names(tmp.receptor.db.all)[i],
			col.desc.vals = c("matched.genes", "unmatched.genes"), 
			ref.vector = list(tmp.receptor.intercelldb.res),
			count.or.identity = TRUE, 
			show.unmatches = TRUE)))
	}
	tmp.receptor.overlap.df <- dplyr::bind_rows(tmp.receptor.overlap.list)
	tmp.receptor.overlap.df$target.names <- factor(tmp.receptor.overlap.df$target.names, levels = names(tmp.receptor.db.all)[1:5])
	# --- plot ---
	pic.receptor.overlap <- ggplot(data = tmp.receptor.overlap.df, 
		aes(x = target.names, y = result.cnt, fill = cmp.types))
	pic.receptor.overlap <- pic.receptor.overlap + 
		geom_col(position = position_fill(reverse = FALSE)) + 
		scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
		#scale_fill_brewer(palette = 1, direction = -1) + 
		scale_fill_manual(values = c("#2172B5", "green")) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 7.2, height = 6)
	ggsave("receptor-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	# stack result
	tmp.receptor.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.receptor.db.all[1:5], tmp.receptor.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("receptor-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# - done
# part 2: overlap is ok, still need to check if InterCellDB lost some receptor in definition 
	tmp.receptor.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.receptor.notin.ref.all <- c(tmp.receptor.notin.ref.all, list(setdiff(tmp.receptor.db.all[[i]], tmp.receptor.intercelldb.res)))
	}
	names(tmp.receptor.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/receptor-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.receptor.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.receptor.notin.ref.all, 
		filename = "../anal-peer-db-Result/GeneFunc/receptor-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done


## option 1.4: Ligand by themselves' definition
# part 1: percent or stack, remapping against InterCellDB
	tmp.ligand.db.all <- list(tmp.ligand.cpdb.res, tmp.ligand.nichenet.res, tmp.ligand.italk.res, 
		tmp.ligand.scsr.res, tmp.ligand.cellchat.res, tmp.ligand.intercelldb.res)
	names(tmp.ligand.db.all) <- kDB.namelist
	# ligand percent matching to InterCellDB
	tmp.ligand.overlap.list <- list()
	for (i in 1:5) {
		tmp.ligand.overlap.list <- c(tmp.ligand.overlap.list, list(Ctrl.percent.plot.add.data(tmp.ligand.db.all[[i]], names(tmp.ligand.db.all)[i],
			col.desc.vals = c("matched.genes", "unmatched.genes"), 
			ref.vector = list(tmp.ligand.intercelldb.res),
			count.or.identity = TRUE, 
			show.unmatches = TRUE)))
	}
	tmp.ligand.overlap.df <- dplyr::bind_rows(tmp.ligand.overlap.list)
	tmp.ligand.overlap.df$target.names <- factor(tmp.ligand.overlap.df$target.names, levels = names(tmp.ligand.db.all)[1:5])
	# --- plot ---
	pic.ligand.overlap <- ggplot(data = tmp.ligand.overlap.df, 
		aes(x = target.names, y = result.cnt, fill = cmp.types))
	pic.ligand.overlap <- pic.ligand.overlap + 
		geom_col(position = position_fill(reverse = FALSE)) + 
		scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
		#scale_fill_brewer(palette = 1, direction = -1) + 
		scale_fill_manual(values = c("#2172B5", "green")) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 7.2, height = 6)
	ggsave("ligand-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	# stack result
	tmp.ligand.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.ligand.db.all[1:5], tmp.ligand.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("ligand-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# part 2: as overlap is not that good, check if InterCellDB lost some ligand in definition
	tmp.ligand.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.ligand.notin.ref.all <- c(tmp.ligand.notin.ref.all, list(setdiff(tmp.ligand.db.all[[i]], tmp.ligand.intercelldb.res)))
	}
	names(tmp.ligand.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/ligand-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.ligand.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.ligand.notin.ref.all, 
		filename = "../anal-peer-db-Result/GeneFunc/ligand-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done


## option 1.5: Integrin by themselves' definition and if not defined, use the definition of InterCellDB
# part 1: percent or stack, remapping against InterCellDB
	tmp.integrin.db.all <- list(tmp.integrin.cpdb.res, tmp.integrin.nichenet.res, tmp.integrin.italk.res, 
		tmp.integrin.scsr.res, tmp.integrin.cellchat.res, tmp.integrin.intercelldb.res)
	names(tmp.integrin.db.all) <- kDB.namelist
	tmp.integrin.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.integrin.db.all[1:5], tmp.integrin.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("integrin-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	tmp.integrin.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.integrin.db.all[1:5], tmp.integrin.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 7.2, height = 6)
	ggsave("integrin-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# - done
# part 2: as overlap is perfectly good, overlap is still needed to show
	tmp.integrin.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.integrin.notin.ref.all <- c(tmp.integrin.notin.ref.all, list(setdiff(tmp.integrin.db.all[[i]], tmp.integrin.intercelldb.res)))
	}
	names(tmp.integrin.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/integrin-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.integrin.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.integrin.notin.ref.all, 
		filename = "../anal-peer-db-Result/GeneFunc/integrin-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done


## option 1.6: G-protein Coupled Receptor. ALL by definition of InterCellDB
# part 1: percent or stack, remapping against InterCellDB
	tmp.gpcR.db.all <- list(tmp.gpcR.cpdb.res, tmp.gpcR.nichenet.res, tmp.gpcR.italk.res, 
		tmp.gpcR.scsr.res, tmp.gpcR.cellchat.res, tmp.gpcR.intercelldb.res)
	names(tmp.gpcR.db.all) <- kDB.namelist
	tmp.gpcR.percent.match.plot.stack <- Ctrl.percent.plot.2item(tmp.gpcR.db.all[1:5], tmp.gpcR.intercelldb.res,
		stack.or.fill = TRUE,  
		scale.y.add = scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
	dev.new(width = 7.2, height = 6)
	ggsave("gpcR-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
	tmp.gpcR.percent.match.plot.fill <- Ctrl.percent.plot.2item(tmp.gpcR.db.all[1:5], tmp.gpcR.intercelldb.res,
		stack.or.fill = FALSE)
	dev.new(width = 7.2, height = 6)
	ggsave("gpcR-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GeneFunc/", device = "pdf")
	dev.off()
# - done
# part 2: as overlap is perfectly good, overlap is still needed to show
	tmp.gpcR.notin.ref.all <- list()
	for (i in 1:5) {
		tmp.gpcR.notin.ref.all <- c(tmp.gpcR.notin.ref.all, list(setdiff(tmp.gpcR.db.all[[i]], tmp.gpcR.intercelldb.res)))
	}
	names(tmp.gpcR.notin.ref.all) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GeneFunc/gpcR-human-notin-intercelldb-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.gpcR.notin.ref.all), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.gpcR.notin.ref.all, 
		filename = "../anal-peer-db-Result/GeneFunc/gpcR-human-notin-intercelldb-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
# - done



# --- 3rd ---
## get gene definition count upon subcellular location and function 
#  to collect all infos from other files

# Location summary (NOT USED after 2021.01.31)
	merge.all.scloc.list <- list(splits.cellphone.scloc, splits.nichenet.scloc, splits.italk.scloc, splits.scsr.scloc, splits.cellchat.scloc, splits.intercelldb.scloc)
	names(merge.all.scloc.list) <- kDB.namelist
	merge.all.scloc.df <- data.frame(DB = kDB.namelist, 
		scloc.kinds = sapply(merge.all.scloc.list, USE.NAMES = FALSE, function(x) {
			length(x)
			}), 
		stringsAsFactors = FALSE)
	ggplot(data = merge.all.scloc.df) + geom_col(aes(x = DB, y = scloc.kinds)) + 
		scale_x_discrete(breaks = kDB.namelist, limits = kDB.namelist)
# - done
# Molecular Function summary (NOT USED after 2021.01.31)
	merge.all.mfunc.list <- list(splits.cellphone.mfunc, splits.nichenet.mfunc, splits.italk.mfunc, splits.scsr.mfunc, splits.cellchat.mfunc, splits.intercelldb.mfunc)
	names(merge.all.scloc.list) <- kDB.namelist
	merge.all.mfunc.df <- data.frame(DB = kDB.namelist, 
		mfunc.kinds = sapply(merge.all.mfunc.list, USE.NAMES = FALSE, function(x) {
			length(x)
			}), 
		stringsAsFactors = FALSE)
	ggplot(data = merge.all.mfunc.df) + geom_col(aes(x = DB, y = mfunc.kinds)) + 
		scale_x_discrete(breaks = kDB.namelist, limits = kDB.namelist)
# - done






##### For Describing PAIRS - Human
# ------
# ------
cmp.pairs.db.list <- list(CellPhoneDB = remap.cellphone.pairs.human,
	NicheNet = remap.nichenet.pairs.human,
	iTALK = remap.italk.pairs.human,
	SingleCellSignalR = remap.scsr.pairs.human,
	CellChatDB = remap.cellchat.pairs.human)


### section 1 - collect all pairs DB and compare between each other
#
## way deprecated - Venn Diagram
library(VennDiagram) 
# to avoid venn.diagram output the log file
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# [ALL Removed] For showing overlaps of gene pairs
# 


## way 1 - percentage matching plot
# for showing matches in 3 independent databases
tmp.human.gp.match.colnames <- c("match.in.exp", "match.in.know", "match.in.pred", "unmatch.all")
tmp.human.gp.ref.pairs <- list(ref.intercelldb.pairs.human.experiments, ref.intercelldb.pairs.human.knowledge, ref.intercelldb.pairs.human.prediction)
#
src.1d2.names <- c("CellPhoneDB", "NicheNet", "iTALK", "SingleCellSignalR", "CellChatDB")

tmp.mix3.list <- list()
for (i in seq_along(cmp.pairs.db.list)) {
	tmp.mix3.list <- c(tmp.mix3.list, 
		list(Ctrl.percent.plot.add.data(cmp.pairs.db.list[[i]], names(cmp.pairs.db.list)[i], 
			col.desc.vals = tmp.human.gp.match.colnames, ref.vector = tmp.human.gp.ref.pairs)))
}
tmp.mix3.df <- dplyr::bind_rows(tmp.mix3.list)
colnames(tmp.mix3.df) <- c("source", "match.types", "match.cnt")

# --- plot ---
	# fill pattern
	pic.1d2.prim <- ggplot(data = tmp.mix3.df, 
		aes(x = source, y = match.cnt, fill = match.types))
	pic.1d2.prim <- pic.1d2.prim + 
		geom_col(position = position_fill(reverse = FALSE)) + 
		#geom_hline(yintercept = 0.2, linetype = "dashed", alpha = 0.8) + 
		scale_x_discrete(limits = src.1d2.names, breaks = src.1d2.names) + 
		scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(mult = c(0, .1))) + 
		scale_fill_brewer(palette = 1, direction = -1) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 7.2, height = 6)
	ggsave("gpairs-human-overlap-percent.pdf", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
	dev.off()

	# stack pattern
	pic.1d2.prim <- ggplot(data = tmp.mix3.df, 
		aes(x = source, y = match.cnt, fill = match.types))
	pic.1d2.prim <- pic.1d2.prim + 
		geom_col(position = position_stack(reverse = FALSE)) +  
		scale_x_discrete(limits = src.1d2.names, breaks = src.1d2.names) + 
		scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
		scale_fill_brewer(palette = 1, direction = -1) + 
		xlab("Source") + ylab("Matches") +
		theme_cowplot(12)
	dev.new(width = 7.2, height = 6)
	ggsave("gpairs-human-overlap-stack.pdf", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
	dev.off()


## way 2 - use euler plot to show the overlap of all other DBs with all-DB & 3 sub-DBs of InterCellDB
# [NOTE] as other 5 DBs are mostly experimentally validated gene pairs,
# 		 so here, only compare the sub-DB[exp] of InterCellDB to them

# part w2.0: all-DB
	# sub1: notin overlap
	tmp.gpairs.notin.all.list <- list()
	for (i in 1:5) {
		tmp.gpairs.notin.all.list <- c(tmp.gpairs.notin.all.list, list(setdiff(cmp.pairs.db.list[[i]], ref.intercelldb.pairs.human)))
	}
	names(tmp.gpairs.notin.all.list) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/gpairs-all-human-notin-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.gpairs.notin.all.list), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		#labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.gpairs.notin.all.list, 
		filename = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/gpairs-all-human-notin-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
	rm(tmp.gpairs.notin.all.list)

# part w2.1: Experiments sub-DB
	# sub1: notin overlap
	tmp.gpairs.notin.exp.list <- list()
	for (i in 1:5) {
		tmp.gpairs.notin.exp.list <- c(tmp.gpairs.notin.exp.list, list(setdiff(cmp.pairs.db.list[[i]], ref.intercelldb.pairs.human.experiments)))
	}
	names(tmp.gpairs.notin.exp.list) <- kDB.namelist[1:5]
	# use euler
	pdf(file = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/gpairs-exp-human-notin-overlap.pdf", 
		width = 9, height = 6)
	plot(euler(tmp.gpairs.notin.exp.list), fills = kDB.favor.color[1:5], 
		edges = list(col = "white", lwd = 0.01, alpha = 0.1), 
		#labels = list(fontfamily = "sans", cex = 0.8), 
		quantities = list(cex = 0.8), 
		legend = TRUE)
	dev.off()
	# use venn
	venn.diagram(tmp.gpairs.notin.exp.list, 
		filename = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/gpairs-exp-human-notin-overlap-venn.png", 
		width = 3000, height = 3000, margin = 0.1)
	rm(tmp.gpairs.notin.exp.list)

# part w2.2: Knowledge sub-DB
 # NOT RUN

# part w2.3: Prediction sub-DB
 # NOT RUN



### EXTRA step, check database validity
if (FALSE) {  # WRONG see following FATAL error
	tmp.cellphone.uq2ref <- setdiff(remap.cellphone.pairs.human, ref.intercelldb.pairs.human)
	tmp.nichenet.uq2ref <- setdiff(remap.nichenet.pairs.human, ref.intercelldb.pairs.human)
	tmp.italk.uq2ref <- setdiff(remap.italk.pairs.human, ref.intercelldb.pairs.human)
	tmp.scsr.uq2ref <- setdiff(remap.scsr.pairs.human, ref.intercelldb.pairs.human)
	tmp.cellchat.uq2ref <- setdiff(remap.cellchat.pairs.human, ref.intercelldb.pairs.human)

	tmp.uq2ref.list <- list(
		CellPhoneDB = tmp.cellchat.uq2ref,     # FATAL error detected here!
		NicheNet = tmp.nichenet.uq2ref,
		iTALK = tmp.italk.uq2ref,
		SingleCellSignalR = tmp.scsr.uq2ref,   
		CellChatDB = tmp.cellchat.uq2ref)      # FATAL error detected here!

	venn.diagram(x = tmp.uq2ref.list,
		filename = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/cmp-human-other5-uq2ref.png", imagetype = "png")

	tmp.other5.uq2ref.overlap <- calculate.overlap(tmp.uq2ref.list)

	# shared level5
	tmp.other5.uq2ref.overlap$a31  # length 3

	# shared level4 -1 except NicheNet
	tmp.other5.uq2ref.overlap$a27  # length 1

	# shared level3 -1 except CellPhoneDB, CellChatDB
	tmp.other5.uq2ref.overlap$a24  # length 29
	tmp.uq2ref.level3.1 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a24, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level3 -2 except iTALK, SingleCellSignalR
	tmp.other5.uq2ref.overlap$a20  # length 23
	tmp.uq2ref.level3.2 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a20, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level2 -1 except CellPhoneDB, NicheNet, CellChatDB
	tmp.other5.uq2ref.overlap$a13  # length 77
	tmp.uq2ref.level2.1 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a13, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level2 -2 except NicheNet, iTALK, SingleCellSignalR
	tmp.other5.uq2ref.overlap$a7  # length 70
	tmp.uq2ref.level2.2 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a7, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level1 -1 is SingleCellSignalR
	tmp.other5.uq2ref.overlap$a4  # length 10
	tmp.uq2ref.level1.1 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a4, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level1 -2 is iTALK
	tmp.other5.uq2ref.overlap$a3  # length 22
	tmp.uq2ref.level1.2 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a3, kgp.cut, c("Gene.A", "Gene.B"))

	# shared level1 -3 is NicheNet
	tmp.other5.uq2ref.overlap$a2  # length 651
	tmp.uq2ref.level1.3 <- Tool.SplitToGenDataFrame(tmp.other5.uq2ref.overlap$a2, kgp.cut, c("Gene.A", "Gene.B"))
}



#### section 2 - use some specific genes to see the related pairs in every database
## In this section, all sub-DBs in InterCellDB will be fully concerned
#

# [WARNING] grep method not exactly get gene names matched right, like getting "CCL5>CCR5", if "CCL5>ABCCR5", it will be fetched back but not expectedly
Ctrl.FindRelatedGenePairs <- function(  
	onegene,
	used.cat.pairs.list,  # to explore related gene pairs in which database. The format of gene pairs is A>B, ">" is defined by kgp.cut
	gp.cut  # the char used to gap the 2 genes in one gene pair
) {
	tmp.res <- lapply(used.cat.pairs.list, onegene = onegene, gp.cut = gp.cut, 
		function(x, onegene, gp.cut) {
			tmp.conv.gp <- grep(paste0("^", as.character(onegene)[[1]], gp.cut), x, value = TRUE)  # restrict only use 1 gene
			tmp.rev.gp  <- grep(paste0(gp.cut, as.character(onegene)[[1]], "$"), x, value = TRUE)		
			length(unique(c(tmp.conv.gp, tmp.rev.gp)))
		})
	return(as.integer(unlist(tmp.res)))
}

## Cytokine
if (FALSE) {  # NOT USED yet, and after 2021.01.31
	# collect all cytokine genes, use the all shared definition
	tmp.cytokine.for.gpairs <- list(tmp.cytokine.cpdb.res, tmp.cytokine.nichenet.res, tmp.cytokine.italk.res, 
			tmp.cytokine.scsr.res, tmp.cytokine.cellchat.res, tmp.cytokine.intercelldb.res)
	tmp.cytokine.all.shared <- Reduce(intersect, tmp.cytokine.for.gpairs)

	# sub 1: use InterCellDB Human all
	tmp.gpairs.cmp.cytokine.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human))
	names(tmp.gpairs.cmp.cytokine.list) <- c(kDB.namelist[1:5], "InterCellDB.human.all")

	# sub 2: use InterCellDB Human Experiments
	tmp.gpairs.cmp.cytokine.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.experiments))
	names(tmp.gpairs.cmp.cytokine.list) <- c(kDB.namelist[1:5], "InterCellDB.human.exps")

	# sub 3: use InterCellDB Human Knowledge
	tmp.gpairs.cmp.cytokine.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.knowledge))
	names(tmp.gpairs.cmp.cytokine.list) <- c(kDB.namelist[1:5], "InterCellDB.human.know")

	# sub 4: use InterCellDB Human Prediction
	# ...

	# --- Shared code ---
	tmp.cytokine.len.calc <- sapply(tmp.cytokine.all.shared, USE.NAMES = TRUE, 
		used.cat.pairs.list = tmp.gpairs.cmp.cytokine.list, gp.cut = kgp.cut, 
		function(x, used.cat.pairs.list, gp.cut) {
			Ctrl.FindRelatedGenePairs(x, used.cat.pairs.list, gp.cut)
		})
	# log transform result
	tmp.cytokine.len.trans.log <- log(tmp.cytokine.len.calc + 1)
	rownames(tmp.cytokine.len.trans.log) <- kDB.namelist
	tmp.cytokine.len.res.df <- data.frame(DBname = rep(kDB.namelist, times = length(tmp.cytokine.all.shared)),
		Gene = rep(tmp.cytokine.all.shared, each = length(kDB.namelist)), 
		pairs.log.cnt = as.numeric(tmp.cytokine.len.trans.log),
		stringsAsFactors = FALSE
	)
	ggplot(data = tmp.cytokine.len.res.df) + 
		geom_line(aes(x = Gene, y = pairs.log.cnt, group = DBname, colour = DBname)) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	rm(tmp.gpairs.cmp.cytokine.list)
}

## Growth factor
if (FALSE) {  # NOT USED yet, and after 2021.01.31
	# collect all cytokine genes, use the all shared definition
	tmp.growfactor.all.shared <- intersect(intersect(tmp.growfactor.cpdb.res, tmp.growfactor.italk.res), intersect(tmp.growfactor.cellchat.res, tmp.growfactor.intercelldb.res))
	# 
	tmp.growfactor.len.calc <- sapply(tmp.growfactor.all.shared, USE.NAMES = TRUE, 
		used.cat.pairs.list = tmp.all.db.list, gp.cut = kgp.cut, 
		function(x, used.cat.pairs.list, gp.cut) {
			Ctrl.FindRelatedGenePairs(x, used.cat.pairs.list, gp.cut)
		})
	# log transform result
	tmp.growfactor.len.trans.log <- log(tmp.growfactor.len.calc + 1)
	rownames(tmp.growfactor.len.trans.log) <- kDB.namelist
	tmp.growfactor.len.res.df <- data.frame(DBname = rep(kDB.namelist, times = length(tmp.growfactor.all.shared)),
		Gene = rep(tmp.growfactor.all.shared, each = length(kDB.namelist)), 
		pairs.log.cnt = as.numeric(tmp.growfactor.len.trans.log),
		stringsAsFactors = FALSE
	)
	ggplot(data = tmp.growfactor.len.res.df) + 
		geom_line(aes(x = Gene, y = pairs.log.cnt, group = DBname, colour = DBname)) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}


## Receptor
	tmp.receptor.for.gpairs <- list(tmp.receptor.cpdb.res, tmp.receptor.nichenet.res, tmp.receptor.italk.res, 
			tmp.receptor.scsr.res, tmp.receptor.cellchat.res, tmp.receptor.intercelldb.res)
	tmp.receptor.all.shared <- unique(Reduce(intersect, tmp.receptor.for.gpairs))

	# sub 1: use InterCellDB Human all
	tmp.gpairs.cmp.receptor.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[1:5], "InterCellDB.human.all")
	tmp.receptor.run.id <- 1

	# sub 2: use InterCellDB Human Experiments
	tmp.gpairs.cmp.receptor.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.experiments))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[1:5], "InterCellDB.human.exps")
	tmp.receptor.run.id <- 2

	# sub 3: use InterCellDB Human Knowledge
	tmp.gpairs.cmp.receptor.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.knowledge))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[1:5], "InterCellDB.human.know")
	tmp.receptor.run.id <- 3

	# sub 4: use InterCellDB Human Prediction
	tmp.gpairs.cmp.receptor.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.prediction))
	names(tmp.gpairs.cmp.receptor.list) <- c(kDB.namelist[1:5], "InterCellDB.human.pred")
	tmp.receptor.run.id <- 4

# --- Shared Code ---
# !!!
# [!NOTE!] !!! after run db.human.all, sub-DBs can be extracted from it by some math calculation
# !!!
# fetching counts
	prog.bar.receptor.len.calc <- progress::progress_bar$new(total = length(tmp.receptor.all.shared))
	prog.bar.receptor.len.calc$tick(0)
	tmp.receptor.len.calc <- sapply(tmp.receptor.all.shared, USE.NAMES = TRUE, 
		used.cat.pairs.list = tmp.gpairs.cmp.receptor.list, gp.cut = kgp.cut, 
		function(x, used.cat.pairs.list, gp.cut) {
			prog.bar.receptor.len.calc$tick()
			Ctrl.FindRelatedGenePairs(x, used.cat.pairs.list, gp.cut)
		})
	# log transform result
	tmp.receptor.len.trans.log <- log(tmp.receptor.len.calc + 1)
	rownames(tmp.receptor.len.trans.log) <- kDB.namelist
	tmp.receptor.len.res.df <- data.frame(DBname = rep(kDB.namelist, times = length(tmp.receptor.all.shared)),
		Gene = rep(tmp.receptor.all.shared, each = length(kDB.namelist)), 
		pairs.log.cnt = as.numeric(tmp.receptor.len.trans.log),
		stringsAsFactors = FALSE
	)
	# to saveRDS
	tmp.receptor.saveRDS.filename <- c("receptor-gpairs-human-all-len-df-log.rds",
		"receptor-gpairs-human-exps-len-df-log.rds",
		"receptor-gpairs-human-know-len-df-log.rds",
		"receptor-gpairs-human-pred-len-df-log.rds")
	tmp.receptor.saveRDS.filename <- paste0("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", 
		tmp.receptor.saveRDS.filename)
	saveRDS(tmp.receptor.len.res.df, tmp.receptor.saveRDS.filename[tmp.receptor.run.id])
	
# ggsave file name version 1
	tmp.receptor.ggsave.filename <- c("receptor-gpairs-human-all-takein.pdf",
		"receptor-gpairs-human-exps-takein.pdf",
		"receptor-gpairs-human-know-takein.pdf",
		"receptor-gpairs-human-pred-takein.pdf")
# template ggsave code
	dev.new(width = 32, height = 8)
	ggsave(tmp.receptor.ggsave.filename[tmp.receptor.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
	dev.off()

# ggsave file name version 2
	tmp.receptor.ggsave.smooth.filename <- c("receptor-gpairs-human-all-takein-smooth.pdf",
		"receptor-gpairs-human-exps-takein-smooth.pdf",
		"receptor-gpairs-human-know-takein-smooth.pdf",
		"receptor-gpairs-human-pred-takein-smooth.pdf")
# template ggsave code
	dev.new(width = 10, height = 8)
	ggsave(tmp.receptor.ggsave.smooth.filename[tmp.receptor.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
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
	tmp.ligand.for.gpairs <- list(tmp.ligand.cpdb.res, tmp.ligand.nichenet.res, tmp.ligand.italk.res, 
			tmp.ligand.scsr.res, tmp.ligand.cellchat.res, tmp.ligand.intercelldb.res)
	tmp.ligand.all.shared <- unique(Reduce(intersect, tmp.ligand.for.gpairs))

	# sub 1: use InterCellDB Human all
	tmp.gpairs.cmp.ligand.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[1:5], "InterCellDB.human.all")
	tmp.ligand.run.id <- 1

	# sub 2: use InterCellDB Human Experiments
	tmp.gpairs.cmp.ligand.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.experiments))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[1:5], "InterCellDB.human.exps")
	tmp.ligand.run.id <- 2

	# sub 3: use InterCellDB Human Knowledge
	tmp.gpairs.cmp.ligand.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.knowledge))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[1:5], "InterCellDB.human.know")
	tmp.ligand.run.id <- 3

	# sub 4: use InterCellDB Human Prediction
	tmp.gpairs.cmp.ligand.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human.prediction))
	names(tmp.gpairs.cmp.ligand.list) <- c(kDB.namelist[1:5], "InterCellDB.human.pred")
	tmp.ligand.run.id <- 4


# !!! Has Code here
# [!NOTE!] !!! after run db.human.all, sub-DBs can be extracted from it by some math calculation
	tmp.123.gp.ligand.all <- readRDS("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ligand-gpairs-human-all-len-df-log.rds")
	tmp.123.gp.ligand.exps <- readRDS("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ligand-gpairs-human-exps-len-df-log.rds")
	tmp.123.gp.ligand.know <- readRDS("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ligand-gpairs-human-know-len-df-log.rds")
	   #tmp.123.gp.ligand.pred <- readRDS("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/ligand-gpairs-human-pred-len-df-log.rds")
	# to calculate to get pred values
	tmp.123.gp.ligand.pred <- tmp.123.gp.ligand.all
	tmp.inds.gp.ligand.intercelldb <- which(tmp.123.gp.ligand.pred$DBname == "InterCellDB")
	tmp.123.gp.ligand.pred.reserve <- tmp.123.gp.ligand.pred[setdiff(seq_len(nrow(tmp.123.gp.ligand.pred)), tmp.inds.gp.ligand.intercelldb), ]
	# modifies the values of InterCellDB
	tmp.123.gp.ligand.pred.to.modf <- tmp.123.gp.ligand.pred[tmp.inds.gp.ligand.intercelldb, ]
	tmp.123.gp.ligand.all.to.modf <- tmp.123.gp.ligand.all[which(tmp.123.gp.ligand.all$DBname == "InterCellDB"), ]
	tmp.123.gp.ligand.exps.to.modf <- tmp.123.gp.ligand.exps[which(tmp.123.gp.ligand.exps$DBname == "InterCellDB"), ]
	tmp.123.gp.ligand.know.to.modf <- tmp.123.gp.ligand.know[which(tmp.123.gp.ligand.know$DBname == "InterCellDB"), ]
	#
	tmp.123.gp.ligand.pred.to.modf <- tmp.123.gp.ligand.pred.to.modf[order(tmp.123.gp.ligand.pred.to.modf$Gene), ]
	tmp.123.gp.ligand.all.to.modf <- tmp.123.gp.ligand.all.to.modf[order(tmp.123.gp.ligand.all.to.modf$Gene), ]
	tmp.123.gp.ligand.exps.to.modf <- tmp.123.gp.ligand.exps.to.modf[order(tmp.123.gp.ligand.exps.to.modf$Gene), ]
	tmp.123.gp.ligand.know.to.modf <- tmp.123.gp.ligand.know.to.modf[order(tmp.123.gp.ligand.know.to.modf$Gene), ]
	# exp back to original count
	tmp.123.gp.ligand.all.to.modf[, "pairs.log.cnt"] <- exp(tmp.123.gp.ligand.all.to.modf[, "pairs.log.cnt"]) - 1
	tmp.123.gp.ligand.exps.to.modf[, "pairs.log.cnt"] <- exp(tmp.123.gp.ligand.exps.to.modf[, "pairs.log.cnt"]) - 1
	tmp.123.gp.ligand.know.to.modf[, "pairs.log.cnt"] <- exp(tmp.123.gp.ligand.know.to.modf[, "pairs.log.cnt"]) - 1
	# calculate the subtraction
	tmp.123.gp.ligand.pred.to.modf[, "pairs.log.cnt"] <- round(tmp.123.gp.ligand.all.to.modf[, "pairs.log.cnt"] - tmp.123.gp.ligand.exps.to.modf[, "pairs.log.cnt"] - tmp.123.gp.ligand.know.to.modf[, "pairs.log.cnt"])
	# log-transform the result
	tmp.123.gp.ligand.pred.to.modf[, "pairs.log.cnt"] <- log(tmp.123.gp.ligand.pred.to.modf[, "pairs.log.cnt"] + 1)
	# merge result
	tmp.123.gp.ligand.pred.to.modf.res <- rbind(tmp.123.gp.ligand.pred.to.modf, tmp.123.gp.ligand.pred.reserve)
	# get it to be len.res.df
	#  tmp.ligand.len.res.df <- tmp.123.gp.ligand.pred.to.modf.res
# !!!


# --- Shared Code ---
# fetching counts
	prog.bar.ligand.len.calc <- progress::progress_bar$new(total = length(tmp.ligand.all.shared))
	prog.bar.ligand.len.calc$tick(0)
	tmp.ligand.len.calc <- sapply(tmp.ligand.all.shared, USE.NAMES = TRUE, 
		used.cat.pairs.list = tmp.gpairs.cmp.ligand.list, gp.cut = kgp.cut, 
		function(x, used.cat.pairs.list, gp.cut) {
			prog.bar.ligand.len.calc$tick()
			Ctrl.FindRelatedGenePairs(x, used.cat.pairs.list, gp.cut)
		})
	# log transform result
	tmp.ligand.len.trans.log <- log(tmp.ligand.len.calc + 1)
	rownames(tmp.ligand.len.trans.log) <- kDB.namelist
	tmp.ligand.len.res.df <- data.frame(DBname = rep(kDB.namelist, times = length(tmp.ligand.all.shared)),
		Gene = rep(tmp.ligand.all.shared, each = length(kDB.namelist)), 
		pairs.log.cnt = as.numeric(tmp.ligand.len.trans.log),
		stringsAsFactors = FALSE
	)
	# to saveRDS
	tmp.ligand.saveRDS.filename <- c("ligand-gpairs-human-all-len-df-log.rds",
		"ligand-gpairs-human-exps-len-df-log.rds",
		"ligand-gpairs-human-know-len-df-log.rds",
		"ligand-gpairs-human-pred-len-df-log.rds")
	tmp.ligand.saveRDS.filename <- paste0("../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", 
		tmp.ligand.saveRDS.filename)
	saveRDS(tmp.ligand.len.res.df, tmp.ligand.saveRDS.filename[tmp.ligand.run.id])
	
# ggsave file name version 1
	tmp.ligand.ggsave.filename <- c("ligand-gpairs-human-all-takein.pdf",
		"ligand-gpairs-human-exps-takein.pdf",
		"ligand-gpairs-human-know-takein.pdf",
		"ligand-gpairs-human-pred-takein.pdf")
# template ggsave code
	dev.new(width = 32, height = 8)
	ggsave(tmp.ligand.ggsave.filename[tmp.ligand.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
	dev.off()

# ggsave file name version 2
	tmp.ligand.ggsave.smooth.filename <- c("ligand-gpairs-human-all-takein-smooth.pdf",
		"ligand-gpairs-human-exps-takein-smooth.pdf",
		"ligand-gpairs-human-know-takein-smooth.pdf",
		"ligand-gpairs-human-pred-takein-smooth.pdf")
# template ggsave code
	dev.new(width = 10, height = 8)
	ggsave(tmp.ligand.ggsave.smooth.filename[tmp.ligand.run.id], 
		path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "pdf")
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
tmp.go.db.list <- c(cmp.pairs.db.list, list(ref.intercelldb.pairs.human, 
	ref.intercelldb.pairs.human.experiments, ref.intercelldb.pairs.human.knowledge, ref.intercelldb.pairs.human.prediction))
names(tmp.go.db.list) <- c(kDB.namelist[1:5], "InterCellDB.human.all", "InterCellDB.human.exps", "InterCellDB.human.know", "InterCellDB.human.pred")

tmp.go.db.related.genes.list <- lapply(tmp.go.db.list, gp.cut = kgp.cut, 
	function(x, gp.cut) {
		Ctrl.FastOptGenePairs.GetAllRelatedGenes(x, gp.cut)
		})

### version 1: plot heatmap by arrange genes along x-axis
	## GO:0006955 [immune response]
	tmp.go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.immune.response.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.go.db.related.genes.list, tmp.go.immune.response.genes))
	tmp.go.immune.response.id.ref.df <- data.frame(gene = tmp.go.immune.response.genes[order(tmp.go.immune.response.genes)], 
		num.id = seq_along(tmp.go.immune.response.genes),
		stringsAsFactors = FALSE)
	tmp.go.immune.response.counts <- dplyr::left_join(tmp.go.immune.response.counts, tmp.go.immune.response.id.ref.df, by = c("Gene" = "gene"))

	## GO:0001525 [angiogenesis]
	tmp.go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.vessel.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.go.db.related.genes.list, tmp.go.vessel.genes))
	tmp.go.vessel.id.ref.df <- data.frame(gene = tmp.go.vessel.genes[order(tmp.go.vessel.genes)], 
		num.id = seq_along(tmp.go.vessel.genes),  # add on
		stringsAsFactors = FALSE)
	# add on num.id
	tmp.go.vessel.id.ref.df$num.id <- tmp.go.vessel.id.ref.df$num.id + length(tmp.go.immune.response.genes)
	tmp.go.vessel.counts <- dplyr::left_join(tmp.go.vessel.counts, tmp.go.vessel.id.ref.df, by = c("Gene" = "gene"))

	## GO:0001837 [epithelial to mesenchymal transition]
	tmp.go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.emt.counts <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MatchEachGene(tmp.go.db.related.genes.list, tmp.go.emt.genes))
	tmp.go.emt.id.ref.df <- data.frame(gene = tmp.go.emt.genes[order(tmp.go.emt.genes)], 
		num.id = seq_along(tmp.go.emt.genes),  # add on
		stringsAsFactors = FALSE)
	# add on num.id
	tmp.go.emt.id.ref.df$num.id <- tmp.go.emt.id.ref.df$num.id + length(tmp.go.immune.response.genes) + length(tmp.go.vessel.genes)
	tmp.go.emt.counts <- dplyr::left_join(tmp.go.emt.counts, tmp.go.emt.id.ref.df, by = c("Gene" = "gene"))


	## merge counts
	tmp.go.counts.merge <- Reduce(rbind, list(tmp.go.immune.response.counts, tmp.go.vessel.counts, tmp.go.emt.counts))
	# log transform
	tmp.go.counts.merge$pairs.cnt <- log(tmp.go.counts.merge$pairs.cnt + 1)
	tmp.go.counts.merge$DBname <- factor(tmp.go.counts.merge$DBname, 
		levels = names(tmp.go.db.list))

	tmp.go.len.1 <- length(tmp.go.immune.response.genes)
	tmp.go.len.2 <- length(tmp.go.vessel.genes)
	tmp.go.len.3 <- length(tmp.go.emt.genes)
	tmp.go.counts.breaks.IT <- c(1, tmp.go.len.1 / 2, tmp.go.len.1,  
		tmp.go.len.1 + tmp.go.len.2 / 2, tmp.go.len.1 + tmp.go.len.2,  
		tmp.go.len.1 + tmp.go.len.2 + tmp.go.len.3 / 2, tmp.go.len.1 + tmp.go.len.2 + tmp.go.len.3)
	tmp.go.counts.breaks.labels <- tmp.go.counts.breaks.IT
	tmp.go.counts.breaks.labels[c(2,4,6)] <- c("Immune Response", "Angiogenesis", "EMT")

	tmp.plot <- ggplot(data = tmp.go.counts.merge) + 
		geom_raster(aes(x = num.id, y = DBname, fill = pairs.cnt)) + 
		scale_x_continuous(breaks = tmp.go.counts.breaks.IT, 
			labels = tmp.go.counts.breaks.labels) + 
		theme_cowplot(12) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	dev.new(width = 12, height = 8)
	ggsave("go-mix3-gpairs-human.png", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "png")
	dev.off()
### - version 1 plot GO end -


### version 2: plot heatmap by arrange GO terms along x-axis
	## GO:0006955 [immune response]
	tmp.go.immune.response.genes <- Tool.FindGenesFromGO("GO:0006955", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.immune.response.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.go.db.related.genes.list, tmp.go.immune.response.genes))
	tmp.go.immune.response.mergecnt$GOterm <- "immune response"

	## GO:0001525 [angiogenesis]
	tmp.go.vessel.genes <- Tool.FindGenesFromGO("GO:0001525", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.vessel.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.go.db.related.genes.list, tmp.go.vessel.genes))
	tmp.go.vessel.mergecnt$GOterm <- "angiogenesis"

	## GO:0001837 [epithelial to mesenchymal transition]
	tmp.go.emt.genes <- Tool.FindGenesFromGO("GO:0001837", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.emt.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.go.db.related.genes.list, tmp.go.emt.genes))
	tmp.go.emt.mergecnt$GOterm <- "EMT"

	## GO:0006956 [complement activation]
	tmp.go.comp.act.genes <- Tool.FindGenesFromGO("GO:0006956", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.comp.act.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.go.db.related.genes.list, tmp.go.comp.act.genes))
	tmp.go.comp.act.mergecnt$GOterm <- "complement activation"

	## GO:0005125 [cytokine activity]
	tmp.go.cytokine.act.genes <- Tool.FindGenesFromGO("GO:0005125", genes.human.ref.db, go.human.ref.db, 
		go.use.relative = TRUE, go.relative.option = "offspring")[[1]]
	tmp.go.cytokine.act.mergecnt <- dplyr::bind_rows(Ctrl.FastOptGenePairs.MergeGORelCnts(tmp.go.db.related.genes.list, tmp.go.cytokine.act.genes))
	tmp.go.cytokine.act.mergecnt$GOterm <- "cytokine activity"


	# merge result
	tmp.go.mergecnt.collect <- Reduce(rbind, list(tmp.go.immune.response.mergecnt,
		tmp.go.vessel.mergecnt, tmp.go.emt.mergecnt, tmp.go.comp.act.mergecnt, tmp.go.cytokine.act.mergecnt))
	tmp.go.mergecnt.collect$log.cnt <- log(tmp.go.mergecnt.collect$merge.pairs.cnt)
	tmp.go.mergecnt.collect$DBname <- factor(tmp.go.mergecnt.collect$DBname, levels = names(tmp.go.db.list))

	tmp.plot <- ggplot(data = tmp.go.mergecnt.collect) + 
		geom_raster(aes(x = GOterm, y = DBname, fill = log.cnt)) + 
		scale_fill_gradient(low = "yellow", high = "red") + 
		theme_cowplot(16) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	dev.new(width = 8, height = 8)
	ggsave("go-mergecnt-mix3-gpairs-human.png", path = "../anal-peer-db-Result/GenePairs-ref-to-InterCellDB/", device = "png")
	dev.off()
### - version 2 plot GO end -
