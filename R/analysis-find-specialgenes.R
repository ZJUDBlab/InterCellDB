
#' Find special genes in one pair of interacting clusters
#'
#' @description
#' This function is used to find special genes in one pair of interacting clusters. Genes are special
#' if it passes some limitations when comparing to other pairs of interacting clusters. 
#'
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param clusters.onepair.select List. Return value of \code{\link{ExtractTargetOnepairClusters}}.
#' @param merge.confidence.on.diff Numeric. Range (0,1) is available. One gene pair is special when it is unique  
#' all or part of other pairs of interacting clusters. This param gives the percentage number for the "part of" unique.
#' @param merge.confidence.on.shared Numeric. Range (0,1) is available. When one gene pair is not unqiue to single 
#' interacting cluster, it could be evaluated as special when the overall expression level changes is different to 
#' all or part of other pairs of interacting clusters. This param gives the percentage number for the "part of" different.
#' @param twist.fold.change.mul Numeric. It defines the lower and upper bound of acceptable difference of expression level.
#' Gene pairs that are out of this bound will be seen as special gene pairs.
#' @param top.ignored.genes.applier Character. It is used to remove some uncared genes from analysis, and is applied on the 
#' former one in a interacting cluster.
#' @param top.ignored.genes.receiver Character. Like \code{top.ignored.genes.applier}, but it is for the latter one 
#' in a interacting cluster.
#' @param top.n.score.positive Numeric. It specifies the count of top positive scores.
#' @param top.n.score.negative Numeric. It specifies the count of top negative scores.
#' @param option.calc.score Integer. Defining the method use in calculating score. The default settings is: 
#' use \code{sum(all)}. In other cases, If it is 1, use \code{sum(abs(all))}.
#'
#' @details
#' If the pair of interacting clusters is C -> D, then the C will be called applier, 
#' and D will be called reciever.
#'
#' To be noted, as this function is based on interact.pairs.acted, if limits have been
#' put upon the clusters in x-axis or y-axis, the interacting pairs compared will be limited 
#' corresponding to the limits put upon clusters.
#'
#' @return List.
#' \itemize{
#'   \item plot: plot top ranked special genes.
#'   \item tables: a list of 2 data.frames.
#'         \itemize{
#'           \item part.C: it records all significant special genes and their scores for the "applier" part.
#'           \item part.D: it records all significant special genes and their scores for the "reciever" part.
#'		   }
#'   \item genes.top.on.score: a list of 2 list of genes.
#'         \itemize{
#'           \item part.C: it gives the specified number of top ranked special genes for the "applier" part.
#'           \item part.D: it gives the specified number of top ranked special genes for the "reciever" part.
#'         }
#'   \item special.pairs.df: this table records all special gene pairs for this pair of interacting cluters.
#' }
#' 
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
	clusters.onepair.select,
	merge.confidence.on.diff = 0.8,
	merge.confidence.on.shared = 0.8,
	twist.fold.change.mul = 1,
	top.ignored.genes.applier = character(),
	top.ignored.genes.receiver = character(),
	top.n.score.positive = 10,
	top.n.score.negative = 10,
	option.calc.score = 0,
	...
) {
	# restrict the range of params
	merge.confidence.on.diff <- ifelse(merge.confidence.on.diff < 0, 0, merge.confidence.on.diff)
	merge.confidence.on.shared <- ifelse(merge.confidence.on.shared < 0, 0, merge.confidence.on.shared)
	twist.fold.change.mul <- ifelse(twist.fold.change.mul > 0, twist.fold.change.mul, 0)

	# calculate the multiply of fold change of interacting genes
	inside.c.short.interacts <- function(
		tmp.pairs,
		collect.flip = FALSE
	) {
		if (nrow(tmp.pairs) == 0) {
			return(data.frame(interacts.name = character(), interacts.mul.fc = numeric(), ex.id = numeric()))
		}
		tmp.new.pairs.names <- paste(tmp.pairs[, "inter.GeneName.A"], tmp.pairs[, "inter.GeneName.B"], sep = ">")
		tmp.new.eval <- tmp.pairs[, "inter.LogFC.A"] + tmp.pairs[, "inter.LogFC.B"]
		tmp.exid.flip <- ifelse(collect.flip == TRUE, -1, 1)
		tmp.new.exid <- tmp.exid.flip * (1:nrow(tmp.pairs))
		tmp.res.df <- data.frame(interacts.name = tmp.new.pairs.names, interacts.mul.fc = tmp.new.eval, ex.id = tmp.new.exid, stringsAsFactors = FALSE)
		return(tmp.res.df)
	}

	# all pairs
	all.pairs.names <- interact.pairs.acted$name.allpairs
	all.pairs.interacts <- interact.pairs.acted$data.allpairs
	# this pair
	this.clusters <- clusters.onepair.select$clusters.name
	this.cluster.C <- this.clusters[1]
	this.cluster.D <- this.clusters[2]
	this.pair.name <- paste0(this.clusters, collapse = kClustersSplit)
	this.pair.interacts <- clusters.onepair.select$bt.pairs
	this.pair.tgs <- inside.c.short.interacts(this.pair.interacts, collect.flip = TRUE)

	# get other pairs
	other.pairs.names.C <- setdiff(grep(paste0("^", this.cluster.C), all.pairs.names, value = TRUE), this.pair.name)
	other.pairs.names.D <- setdiff(grep(paste0(this.cluster.D, "$"), all.pairs.names, value = TRUE), this.pair.name)

	## tool function
	inside.score.on.selected.interacts <- function(
		target.pairs,
		optionx
	) {
		tmp.C.df <- data.frame(genes = character(), interacts.cnt = numeric(), interact.score = numeric())
		tmp.D.df <- tmp.C.df
		# function for tapply
		tiny.calc <- function(x, optionx) {
			if (optionx == 1) {
				sum(abs(x))
			} else {
				sum(x)
			}
		}
		if (nrow(target.pairs) != 0) {
			# calculate the score of applier cells
			tmp.genes.C <- tapply(1:nrow(target.pairs), target.pairs$inter.GeneName.A, length)
			tmp.C.df <- data.frame(genes = names(tmp.genes.C), interacts.cnt = tmp.genes.C, stringsAsFactors = FALSE)
			tmp.C.geneA.sums <- tapply(target.pairs[, "inter.LogFC.B"], target.pairs$inter.GeneName.A, optionx = optionx, tiny.calc)
			tmp.C.geneA.itself <- tapply(target.pairs[, "inter.LogFC.A"], target.pairs$inter.GeneName.A, mean)
			tmp.C.df$interact.score <- tmp.C.geneA.sums * abs(tmp.C.geneA.itself)
			tmp.C.df$interact.avg.score <- tmp.C.df$interact.score / tmp.C.df$interacts.cnt
			tmp.C.df$interact.marker <- ifelse(tmp.C.geneA.itself > 0, "UPreg", "DNreg")  # marker indicate the A if it is upregulated or not
			# calculate the score of reciever cells
			tmp.genes.D <- tapply(1:nrow(target.pairs), target.pairs$inter.GeneName.B, length)
			tmp.D.df <- data.frame(genes = names(tmp.genes.D), interacts.cnt = tmp.genes.D, stringsAsFactors = FALSE)
			tmp.D.geneB.sums <- tapply(target.pairs[, "inter.LogFC.A"], target.pairs$inter.GeneName.B, optionx = optionx, tiny.calc)
			tmp.D.geneB.itself <- tapply(target.pairs[, "inter.LogFC.B"], target.pairs$inter.GeneName.B, mean)
			tmp.D.df$interact.score <- tmp.D.geneB.sums * abs(tmp.D.geneB.itself)
			tmp.D.df$interact.avg.score <- tmp.D.df$interact.score / tmp.D.df$interacts.cnt
			tmp.D.df$interact.marker <- ifelse(tmp.D.geneB.itself > 0, "UPreg", "DNreg")  # marker indicate the B if it is upregulated or not
		}
		return(list(part.C = tmp.C.df, part.D = tmp.D.df))
	}

	# compare this to ohter under 1 v-line & 1 h-line
	twist.fold.change.mul.up <- twist.fold.change.mul + 1
	twist.fold.change.mul.it <- c(1 / twist.fold.change.mul.up, twist.fold.change.mul.up)

	# compare with C ~ x & x ~ D
	inside.collect.pass.twist.one.pair.each <- function(
		all.pairs.interacts,
		other.pairs.names,
		this.pair.tgs,
		twist.fold.change.mul.it
	) {
		p.x.collect.diff.pairs <- character()
		p.x.collect.shared.pairs <- character()
		p.x.collect.shared.orig.pairs <- character()
		for (p.x in other.pairs.names) {
			p.x.pairs <- all.pairs.interacts[[p.x]]
			p.x.pairs.tgs <- inside.c.short.interacts(p.x.pairs)
			# get diff pairs
			tmp.p.x.diff.inds <- match(setdiff(this.pair.tgs[, "interacts.name"], p.x.pairs.tgs[, "interacts.name"]), this.pair.tgs[, "interacts.name"])
			# get shared pairs
			tmp.p.x.shared.inds <- integer()
			if (nrow(this.pair.tgs) > 0) {
				tmp.p.x.shared.inds <- setdiff(1:nrow(this.pair.tgs), tmp.p.x.diff.inds)
			}
			tmp.shared.pairs <- this.pair.tgs[tmp.p.x.shared.inds, ]
			colnames(tmp.shared.pairs)[which(colnames(tmp.shared.pairs) == "interacts.mul.fc")] <- "this.mul.fc"
			tmp.shared.pairs <- left_join(tmp.shared.pairs[, c("interacts.name", "this.mul.fc")], p.x.pairs.tgs[, c("interacts.name", "interacts.mul.fc")], by = c("interacts.name"))
			# compare all the shared ones' twist
			tmp.shared.eval <- exp(tmp.shared.pairs[, "this.mul.fc"] - tmp.shared.pairs[, "interacts.mul.fc"])
			tmp.shared.eval.within.twist <- (tmp.shared.eval > twist.fold.change.mul.it[1]) & (tmp.shared.eval < twist.fold.change.mul.it[2])
			p.x.collect.diff.pairs <- c(p.x.collect.diff.pairs, this.pair.tgs[tmp.p.x.diff.inds, "interacts.name"])
			p.x.collect.shared.orig.pairs <- c(p.x.collect.shared.orig.pairs, this.pair.tgs[tmp.p.x.shared.inds, "interacts.name"])
			p.x.collect.shared.pairs <- c(p.x.collect.shared.pairs, tmp.shared.pairs[which(tmp.shared.eval.within.twist == FALSE), "interacts.name"])
		}
		return(list(diff = p.x.collect.diff.pairs, shared.orig = p.x.collect.shared.orig.pairs, shared = p.x.collect.shared.pairs))
	}
	#
	p.C.collect.pairs <- inside.collect.pass.twist.one.pair.each(all.pairs.interacts, other.pairs.names.C, this.pair.tgs, twist.fold.change.mul.it)
	p.D.collect.pairs <- inside.collect.pass.twist.one.pair.each(all.pairs.interacts, other.pairs.names.D, this.pair.tgs, twist.fold.change.mul.it)
	#
	p.C.diff.tgs.pnames <- character()
	p.D.diff.tgs.pnames <- character()
	if (length(p.C.collect.pairs$diff) > 0) {
		# for diff pairs, it is special if it only appears in restricted number of interacting clusters
		p.C.collect.diff.dups <- tapply(1:length(p.C.collect.pairs$diff), p.C.collect.pairs$diff, length)
		# collect names of gene pairs in diff pairs
		p.C.diff.tgs.pnames <- names(p.C.collect.diff.dups[which(p.C.collect.diff.dups >= floor(merge.confidence.on.diff * length(other.pairs.names.C)))])
	} 
	if (length(p.D.collect.pairs$diff) > 0) {
		p.D.collect.diff.dups <- tapply(1:length(p.D.collect.pairs$diff), p.D.collect.pairs$diff, length)
		p.D.diff.tgs.pnames <- names(p.D.collect.diff.dups[which(p.D.collect.diff.dups >= floor(merge.confidence.on.diff * length(other.pairs.names.D)))])
	}
	#
	p.C.shared.tgs.pnames <- character()
	p.D.shared.tgs.pnames <- character()
	if (length(p.C.collect.pairs$shared) > 0) {
		# for shared pairs, it is special if it is differently expressed against some percentage of interacting clusters that it appears
		p.C.collect.shared.dups <- tapply(1:length(p.C.collect.pairs$shared), p.C.collect.pairs$shared, length)
		p.C.collect.shared.dups <- p.C.collect.shared.dups[order(names(p.C.collect.shared.dups))]
		# get the shared orig collection
		p.C.collect.shared.orig <- tapply(1:length(p.C.collect.pairs$shared.orig), p.C.collect.pairs$shared.orig, length)
		p.C.collect.shared.orig <- p.C.collect.shared.orig[which(names(p.C.collect.shared.orig) %in% names(p.C.collect.shared.dups))]
		p.C.collect.shared.orig <- p.C.collect.shared.orig[order(names(p.C.collect.shared.orig))]
		# collect names of gene pairs in shared pairs
		p.C.shared.tgs.pnames <- names(p.C.collect.shared.dups[which(p.C.collect.shared.dups >= floor(merge.confidence.on.shared * p.C.collect.shared.orig))])
	}
	if (length(p.D.collect.pairs$shared) > 0) {
		p.D.collect.shared.dups <- tapply(1:length(p.D.collect.pairs$shared), p.D.collect.pairs$shared, length)
		p.D.collect.shared.dups <- p.D.collect.shared.dups[order(names(p.D.collect.shared.dups))]
		p.D.collect.shared.orig <- tapply(1:length(p.D.collect.pairs$shared.orig), p.D.collect.pairs$shared.orig, length)
		p.D.collect.shared.orig <- p.D.collect.shared.orig[which(names(p.D.collect.shared.orig) %in% names(p.D.collect.shared.dups))]
		p.D.collect.shared.orig <- p.D.collect.shared.orig[order(names(p.D.collect.shared.orig))]
		p.D.shared.tgs.pnames <- names(p.D.collect.shared.dups[which(p.D.collect.shared.dups >= floor(merge.confidence.on.shared * p.D.collect.shared.orig))])
	}
	# merge result
	p.C.tgs.pnames <- c(p.C.diff.tgs.pnames, p.C.shared.tgs.pnames)
	p.D.tgs.pnames <- c(p.D.diff.tgs.pnames, p.D.shared.tgs.pnames)
	# get the special pairs, seperately
	this.C.res.df <- this.pair.interacts[match(p.C.tgs.pnames, this.pair.tgs[, "interacts.name"]), ]
	this.D.res.df <- this.pair.interacts[match(p.D.tgs.pnames, this.pair.tgs[, "interacts.name"]), ]
	this.merge.res.df <- DoPartUnique(rbind(this.C.res.df, this.D.res.df), c(1,2))

	# calculate the scores of special pairs
	this.merge.res.tl <- inside.score.on.selected.interacts(this.merge.res.df, option.calc.score)
	this.C.spdf <- this.merge.res.tl$part.C
	this.D.spdf <- this.merge.res.tl$part.D

	# ignore some genes in ranking plot
	this.C.spdf <- this.C.spdf[which(this.C.spdf$genes %in% (setdiff(this.C.spdf$genes, top.ignored.genes.applier))), ]
	this.D.spdf <- this.D.spdf[which(this.D.spdf$genes %in% (setdiff(this.D.spdf$genes, top.ignored.genes.receiver))), ]

	#
	inside.topn.score.sort <- function(df, target.col, n.top.sel, decreasing) {
		#
		df.sub.1 <- df[which(df[, target.col] > 0), ]
		df.sub.1 <- df.sub.1[order(df.sub.1[, target.col], decreasing = decreasing), ]
		if (nrow(df.sub.1) > n.top.sel[1]) {
			df.sub.1 <- df.sub.1[1:(n.top.sel[1]), ]
		}
		#
		df.sub.2 <- df[which(df[, target.col] <= 0), ]
		df.sub.2 <- df.sub.2[order(df.sub.2[, target.col], decreasing = !decreasing), ]
		if (nrow(df.sub.2) > n.top.sel[2]) {
			df.sub.2 <- df.sub.2[1:(n.top.sel[2]), ]
		}
		#
		ret.df <- df.sub.1
		if (nrow(df.sub.1) != 0 && nrow(df.sub.2) != 0) {
			ret.df <- rbind(df.sub.2, df.sub.1[nrow(df.sub.1):1, ])
		} else {
			if (nrow(df.sub.2) != 0) {
				ret.df <- df.sub.2[nrow(df.sub.2):1, ]
			} else {
				if (nrow(df.sub.1) != 0) {
					ret.df <- df.sub.1[nrow(df.sub.1):1, ]
				}
			}
		}
		ret.df
	}

	if (nrow(this.C.spdf) == 0 || nrow(this.D.spdf) == 0) {
		cat("Quit with no speical genes, in interact [", this.cluster.C, " and ", this.cluster.D, "]!\n", 
			"Re-check if too strict restrictions were given!\n")
		return(NULL)
	}

	### draw plots
	## score plots (consider > 0 and < 0)
	this.C.spdf.score <- inside.topn.score.sort(this.C.spdf, "interact.score", c(top.n.score.positive, top.n.score.negative), TRUE)
	this.D.spdf.score <- inside.topn.score.sort(this.D.spdf, "interact.score", c(top.n.score.positive, top.n.score.negative), TRUE)
	
	this.plot.C.score <- ggplot(this.C.spdf.score, aes(x = genes, y = interact.score))
	this.plot.C.score <- this.plot.C.score + 
			geom_col(aes(fill = interact.marker), colour = "black") + 
			scale_fill_manual(name = "UpDn", values = c("grey", "red"), breaks = c("DNreg", "UPreg")) + 
			scale_x_discrete(breaks = this.C.spdf.score$genes, limits = this.C.spdf.score$genes) + 
			coord_flip() + 
			labs(title = this.cluster.C)

	this.plot.D.score <- ggplot(this.D.spdf.score, aes(x = genes, y = interact.score))
	this.plot.D.score <- this.plot.D.score + 
			geom_col(aes(fill = interact.marker), colour = "black") + 
			scale_fill_manual(name = "UpDn", values = c("lightgrey", "green"), breaks = c("DNreg", "UPreg")) + 
			scale_x_discrete(breaks = this.D.spdf.score$genes, limits = this.D.spdf.score$genes) + 
			coord_flip() + 
			labs(title = this.cluster.D)
	#
	this.final.4plots <- plot_grid(this.plot.C.score, this.plot.D.score, ncol = 2, align = "vh")

	#
	list(plot = this.final.4plots, tables = list(part.C = this.C.spdf, part.D = this.D.spdf),
		genes.top.on.score = list(part.C = this.C.spdf.score$genes, part.D = this.D.spdf.score$genes), 
		special.pairs.df = this.merge.res.df)
}





#' Find common changed genes in one cluster
#'
#' @description
#' TODO
#'
#' @param interact.pairs.acted TODO
#' @param cluster.name TODO
#' @param cluster.applier.or.receiver TODO
#' @param common.change.confidence TODO
#' @param cmp.extend.positive TODO
#' @param cmp.extend.negative TODO
#' @param ... [TODO] a lot of param in \code{FindSpecialGenesInOnepairCluster} will be passed
#'
#'
#' @details
#' TODO
#'
#'
#'
#'
#' @return TODO
#'
#'
#'
#' @export
#'
FindCommonChangedGenesInOnepairCluster <- function(
	interact.pairs.acted,
	cluster.name,
	cluster.applier.or.receiver = TRUE,
	common.change.confidence = 0.8,
	cmp.extend.positive = 0.25,
	cmp.extend.negative = 0.25,
	...
) {
	#
	cluster.to.cmp <- cluster.name
	all.clusters <- unique(as.character(unlist(interact.pairs.acted$list.clusters)))
	#
	other.clusters <- setdiff(all.clusters, cluster.to.cmp)
	other.imp.genes.df <- lapply(other.clusters, cluster = cluster.to.cmp, 
		cluster.position = cluster.applier.or.receiver, full.ref = interact.pairs.acted,
		over0.percent = cmp.extend.positive, under0.percent = cmp.extend.negative,
		...,
		function(x, cluster, cluster.position, full.ref, over0.percent, under0.percent, ...) {
			if (cluster.position == TRUE) {
				this.tg <- ExtractTargetOnepairClusters(full.ref, cluster, x)
			} else {
				this.tg <- ExtractTargetOnepairClusters(full.ref, x, cluster)
			}
			#
			this.res <- FindSpecialGenesInOnepairCluster(full.ref, this.tg, ...)
			if (is.null(this.res)) {
				return(list(change.genes = character()))
			} else {
				over0.percent <- ifelse(over0.percent >= 1, 1, ifelse(over0.percent < 0, 0, over0.percent))
				under0.percent <- ifelse(under0.percent >= 1, 1, ifelse(under0.percent < 0, 0, under0.percent))
				
				if (cluster.position == TRUE) {
					tb.change <- this.res$tables$part.C
				} else {
					tb.change <- this.res$tables$part.D
				}
				tb.change <- tb.change[order(tb.change$interact.score, decreasing = TRUE), ]
				tb.up.change <- floor(length(which(tb.change$interact.score >= 0)) * over0.percent)
				tb.dn.change <- floor(length(which(tb.change$interact.score < 0)) * under0.percent)
				list.up.change <- if (tb.up.change < 1) {0:0} else {1:tb.up.change}  # 0:0 to emphasize the result
				list.dn.change <- if (tb.dn.change < 1) {0:0} else {(nrow(tb.change) - tb.dn.change + 1):nrow(tb.change)}
				cut.list <- c(list.up.change, list.dn.change)
				tmp.genes.change <- tb.change[cut.list, "genes"]
				return(list(change.genes = tmp.genes.change))
			}
		})

	tmp.collect <- list(genes = character(), clusters = character())
	for (i in 1:length(other.imp.genes.df)) {
		tmp.collect$genes <- c(tmp.collect$genes, other.imp.genes.df[[i]]$change.genes)
		tmp.collect$clusters <- c(tmp.collect$clusters, rep(other.clusters[i], times = length(other.imp.genes.df[[i]]$change.genes)))
	}

	tmp.res <- tapply(tmp.collect$clusters, tmp.collect$genes, length)
	tmp.res <- tmp.res[order(tmp.res, decreasing = TRUE)]
	tmp.change.genes <- names(tmp.res[which(tmp.res >= common.change.confidence * length(other.clusters))])
	#
	return(tmp.change.genes)
}


