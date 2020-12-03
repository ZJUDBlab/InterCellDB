


#' Find special genes in one pair of interacting clusters
#'
#' @description
#' This function is used to find special genes in one pair of interacting clusters. Genes are special
#' if it passes some limitations when comparing to other pairs of interacting clusters. 
#'
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param target.gene.pairs.df [TODO]
#' @param involved.clusters.pair [TODO]
#' @param to.cmp.clusters.pairs [TODO]
#' @param to.cmp.clusters.pairs.sel.strategy [TODO]
#' @param merge.cnt.quantile [TODO]
#' @param merge.cnt.decreasing [TODO]
#' @param merge.cnt.topn.shown [TODO]
#' @param plot.fill [TODO]
#' @param sep.inside.gene.pair [TODO]
#'
#' @details
#' [TODO]
#'
#'
#' @return List. [TODO]
#'
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#'
#' 
#' @export
#'
FindSpecialGenesInOnepairCluster <- function(
	interact.pairs.acted,
	target.gene.pairs.df,  # get from VEinfos 
	involved.clusters.pair, # = paste("cluster.C", "cluster.D", sep = kClustersSplit),
	to.cmp.clusters.pairs = character(),
	to.cmp.clusters.pairs.sel.strategy = "inter-cross",
	merge.cnt.quantile = 0.5,  # median
	merge.cnt.decreasing = TRUE,
	merge.cnt.topn.shown = 20,
	plot.fill = c("red", "green"),
	sep.inside.gene.pair = "->"
) {
	this.gp.sep <- sep.inside.gene.pair
	# calculate the multiply of fold change of interacting genes
	inside.c.short.interacts <- function(
		tmp.pairs,
		tmp.sep
	) {
		paste(tmp.pairs[, "inter.GeneName.A"], tmp.pairs[, "inter.GeneName.B"], sep = tmp.sep)
	}

	# all pairs
	all.pairs.names <- interact.pairs.acted$name.allpairs
	all.pairs.interacts <- interact.pairs.acted$data.allpairs
	# this pair
	this.pair.name <- involved.clusters.pair
	if ((this.pair.name %in% all.pairs.names) == FALSE) {
		stop("Given name of cluster pairs NOT exist!")
	}
	this.clusters.separate <- strsplit(this.pair.name, split = kClustersSplit, fixed = TRUE)[[1]]
	this.pair.interacts <- inside.c.short.interacts(target.gene.pairs.df, this.gp.sep)
	this.gene.pairs <- data.frame(gp.name = this.pair.interacts, gp.belongs = this.pair.name, stringsAsFactors = FALSE)

	# to compare pairs
	if (length(to.cmp.clusters.pairs) == 0) {  # only if no list is given then use the pre-defined strategy
		if (to.cmp.clusters.pairs.sel.strategy == "inter-cross") {
			other.pairs.names.C <- setdiff(grep(paste0("^", this.clusters.separate[1]), all.pairs.names, value = TRUE), this.pair.name)
			other.pairs.names.D <- setdiff(grep(paste0(this.clusters.separate[2], "$"), all.pairs.names, value = TRUE), this.pair.name)
			to.cmp.clusters.pairs <- unique(c(other.pairs.names.C, other.pairs.names.D))
		} else {
			if (to.cmp.clusters.pairs.sel.strategy == "all-other") {
				to.cmp.clusters.pairs <- setdiff(all.pairs.names, this.pair.name)
			} else {
				stop("Undefined pre-defined strategy used: ", to.cmp.clusters.pairs.sel.strategy)
			}
		}
	}  # else use the directly specified clusters

	# generate gene pairs to compare
	other.gene.pairs.df.list <- lapply(to.cmp.clusters.pairs, 
		all.pairs.interacts = all.pairs.interacts, tmp.sep = this.gp.sep, 
		function(x, all.pairs.interacts, tmp.sep) {
			tmp.pairs.df <- all.pairs.interacts[[x]]
			tmp.pairs.interacts <- inside.c.short.interacts(tmp.pairs.df, tmp.sep)
			if (length(tmp.pairs.interacts) > 0) {
				res.df <- data.frame(gp.name = tmp.pairs.interacts, gp.belongs = x, stringsAsFactors = FALSE)
			} else {
				res.df <- data.frame(gp.name = "", gp.belongs = "", stringsAsFactors = FALSE)
				res.df <- res.df[-1, ]
			}
			res.df
		})
	other.gene.pairs.packed <- bind_rows(other.gene.pairs.df.list)
	
	# get diff pairs counts
	tmp.all.gp.packed <- rbind(this.gene.pairs, other.gene.pairs.packed)
	tmp.diff.cnts <- tapply(tmp.all.gp.packed$gp.belongs, tmp.all.gp.packed$gp.name, length)
	tmp.diff.indexs <- tapply(1:nrow(tmp.all.gp.packed), tmp.all.gp.packed$gp.name, function(x) { x	})

	tmp.diff.res.mat <- sapply(this.pair.interacts, all.gp = tmp.all.gp.packed, 
		diff.cnts = tmp.diff.cnts, diff.indexs = tmp.diff.indexs, 
		function(x, all.gp, diff.cnts, diff.indexs) {
			tmp.i <- which(names(diff.cnts) == x)  # must be the same as it is in diff.indexs
			tmp.cnt <- diff.cnts[tmp.i]
			tmp.indexs <- diff.indexs[[tmp.i]]
			tmp.involved.cps <- all.gp$gp.belongs[tmp.indexs]
			tmp.involved.cps <- tmp.involved.cps[order(tmp.involved.cps)]
			tmp.takein.cps <- paste0(tmp.involved.cps, collapse = " ")
			c(x, tmp.takein.cps, as.character(tmp.cnt))
		})
	tmp.diff.res.mat <- t(tmp.diff.res.mat)
	this.diff.res.df <- data.frame(gp.name = tmp.diff.res.mat[, 1], gp.takein.cp = tmp.diff.res.mat[, 2], 
		gp.uq.cnt = as.integer(tmp.diff.res.mat[, 3]), stringsAsFactors = FALSE)

	# plot preparation
	# plot the median or quantile of the participated genes
	tmp.gp.splits <- as.character(unlist(strsplit(this.diff.res.df$gp.name, split = this.gp.sep, fixed = TRUE)))
	tmp.gp.cnts.match <- rep(this.diff.res.df$gp.uq.cnt, each = 2)
	# sender part
	tmp.seq.sender <- 1:(length(tmp.gp.splits) / 2) * 2 - 1
	tmp.sender.gpq.df <- data.frame(gp.gene.name = tmp.gp.splits[tmp.seq.sender], 
		gp.gene.cnt = tmp.gp.cnts.match[tmp.seq.sender], stringsAsFactors = FALSE)
	this.sender.gpq.collect <- tapply(X = tmp.sender.gpq.df$gp.gene.cnt, INDEX = tmp.sender.gpq.df$gp.gene.name, 
		cnt.base = tmp.sender.gpq.df$gp.gene.cnt, quantile.it = merge.cnt.quantile, 
		FUN = function(x, cnt.base, quantile.it) {
			quantile(cnt.base[x], probs = quantile.it)
		})
	this.sender.gpq.collect <- this.sender.gpq.collect[order(this.sender.gpq.collect, decreasing = merge.cnt.decreasing)]
	# select top n
	tmp.sender.topn <- ifelse(length(this.sender.gpq.collect) > merge.cnt.topn.shown, merge.cnt.topn.shown, length(this.sender.gpq.collect))
	this.sender.gpq.collect <- this.sender.gpq.collect[1:tmp.sender.topn]
	this.sender.gpq.res <- data.frame(gp.gene.name = names(this.sender.gpq.collect), gp.cnt.ql = this.sender.gpq.collect, stringsAsFactors = FALSE)

	# receiver part
	tmp.seq.receiver <- 1:(length(tmp.gp.splits) / 2) * 2
	tmp.receiver.gpq.df <- data.frame(gp.gene.name = tmp.gp.splits[tmp.seq.receiver], 
		gp.gene.cnt = tmp.gp.cnts.match[tmp.seq.receiver], stringsAsFactors = FALSE)
	this.receiver.gpq.collect <- tapply(X = tmp.receiver.gpq.df$gp.gene.cnt, INDEX = tmp.receiver.gpq.df$gp.gene.name, 
		cnt.base = tmp.receiver.gpq.df$gp.gene.cnt, quantile.it = merge.cnt.quantile, 
		FUN = function(x, cnt.base, quantile.it) {
			quantile(cnt.base[x], probs = quantile.it)
		})
	this.receiver.gpq.collect <- this.receiver.gpq.collect[order(this.receiver.gpq.collect, decreasing = merge.cnt.decreasing)]
	# select top n
	tmp.receiver.topn <- ifelse(length(this.receiver.gpq.collect) > merge.cnt.topn.shown, merge.cnt.topn.shown, length(this.receiver.gpq.collect))
	this.receiver.gpq.collect <- this.receiver.gpq.collect[1:tmp.receiver.topn]
	this.receiver.gpq.res <- data.frame(gp.gene.name = names(this.receiver.gpq.collect), gp.cnt.ql = this.receiver.gpq.collect, stringsAsFactors = FALSE)

	# draw plots
	this.plot.sender <- ggplot(this.sender.gpq.res, aes(x = gp.gene.name, y = gp.cnt.ql))
	this.plot.sender <- this.plot.sender + 
			geom_col(fill = plot.fill[1], colour = "black") + 
			scale_x_discrete(breaks = this.sender.gpq.res$gp.gene.name, limits = this.sender.gpq.res$gp.gene.name) + 
			coord_flip() + 
			labs(title = this.clusters.separate[1], y = "Gene Name", x = "Count Quantile")


	this.plot.receiver <- ggplot(this.receiver.gpq.res, aes(x = gp.gene.name, y = gp.cnt.ql))
	this.plot.receiver <- this.plot.receiver + 
			geom_col(fill = plot.fill[2], colour = "black") + 
			scale_x_discrete(breaks = this.receiver.gpq.res$gp.gene.name, limits = this.receiver.gpq.res$gp.gene.name) + 
			coord_flip() + 
			labs(title = this.clusters.separate[2], y = "Gene Name", x = "Count Quantile")
	#
	this.m.plot <- plot_grid(this.plot.sender, this.plot.receiver, ncol = 2, align = "vh")

	# before return, generate res table
	tmp.gpname.splits.df <- Tool.SplitToGenDataFrame(this.diff.res.df[, "gp.name"], to.split.by = this.gp.sep, 
		res.colnames = c("inter.GeneName.A", "inter.GeneName.B"))
	this.diff.res.df <- cbind(tmp.gpname.splits.df, this.diff.res.df[, c("gp.name", "gp.takein.cp", "gp.uq.cnt")])
	
	return(list(plot = this.m.plot, table = list(special.table = this.diff.res.df)))
}

