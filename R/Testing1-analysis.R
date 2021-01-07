
#' Find special genes in one pair of interacting clusters
#'
#' @description
#' This function is used to find special genes in one pair of interacting clusters. Genes are special
#' if it passes some limitations when comparing to other pairs of interacting clusters. 
#'
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param clusters.onepair.select List. Return value of \code{\link{ExtractTargetOnepairClusters}}.
#' @param evaluation.method Character. This defines the method used upon evaluation of ranks of gene pairs, and 
#' it supports 2 methods currently, "special-score", and "special-index", see details for the differences between these 2 methods.
#' @param merge.confidence.on.diff Numeric. Range (0,1) is available. One gene pair is special when it is unique  
#' all or part of other pairs of interacting clusters. This param gives the percentage number for the "part of" unique.
#' @param merge.confidence.on.shared Numeric. Range (0,1) is available. When one gene pair is not unqiue to single 
#' interacting cluster, it could be evaluated as special when the overall expression level changes is different to 
#' all or part of other pairs of interacting clusters. This param gives the percentage number for the "part of" different.
#' @param twist.fold.change.mul Numeric. It defines the lower and upper bound of acceptable difference of expression level.
#' Gene pairs that are out of this bound will be seen as special gene pairs.
#' @param top.ignored.genes.sender Character. It is used to remove some uncared genes from analysis, and is applied on the 
#' former one in a interacting cluster.
#' @param top.ignored.genes.receiver Character. Like \code{top.ignored.genes.sender}, but it is for the latter one 
#' in a interacting cluster.
#' @param top.n.col.positive Numeric. It specifies the count of top positive evaluation results.
#' @param top.n.col.negative Numeric. It specifies the count of top negative evaluation results. To be noted, If the 
#' \code{evaluation.method} is "special-index", then no negative results will be found.
#' @param option.calc.score Integer. Defining the method use in calculating score. The default settings is: 
#' use \code{sum(all)}. In other cases, If it is 1, use \code{sum(abs(all))}.
#'
#' @details
#' If the pair of interacting clusters is C -> D, then the C will be called sender, 
#' and D will be called receiver.
#'
#' To be noted, as this function is based on \code{interact.pairs.acted}, if limits have been
#' put upon the clusters in x-axis(sender) or y-axis(receiver), the interacting pairs compared will be limited 
#' corresponding to the limits put upon clusters.
#'
#' \bold{evaluation.method}: "special-score" & "special-index". To emphasize, whatever method is used, the result is depending on interacting clusters, 
#' which means the special genes and their ranks will be different if different limits are applied upon interacting clusters, let alone different clusters.
#'
#' \itemize{
#'   \item "special-score": The logFCs(log value of the fold change of genes) is important in this method. 
#'   When calculates, every gene will get all gene pairs that it participates merged to get one single value, called Effect Value.
#'   The Effect Value is calculating in the formula like the following:
#'       If the gene is A, and it is the DEG of the cell cluster Gx, the Effect Value of it will be denoted as E(Gx, A). The gene pairs
#'       that involve A will be denoted as p1 ... pn, and the interacting counterpart gene will be denoted as Bi(i = 1...n). Then, 
#'
#'       E(Gx, A) = exp(logFC(A)) * SUM(logFC(Bi))(i=1...n).
#'   
#'   \item "special-index": This method doesn't consider the logFCs, and only takes the frequency of the genes into account. 
#'   For example, if the 2 interacting clusters are G3 and G5, and gene pairs between these 2 clusters are denoted as p1...pn, 
#'   and if gene A is in some pairs pi...pj, and the occurrence frequency of A is the count of pi...pj. 
#' }
#'
#' @return List.
#' \itemize{
#'   \item plot: plot top ranked special genes.
#'   \item tables: a list of 2 data.frames.
#'         \itemize{
#'           \item part.C: it records all significant special genes and their scores for the "sender" part.
#'           \item part.D: it records all significant special genes and their scores for the "reciever" part.
#'		   }
#'   \item genes.top.on.col: a list of 2 list of genes.
#'         \itemize{
#'           \item part.C: it gives the specified number of top ranked special genes for the "sender" part.
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
# NOT export
Testing.FindSpecialGenesInOnepairCluster <- function(
	interact.pairs.acted,
	clusters.onepair.select,
	evaluation.method = "special-score",
	merge.confidence.on.diff = 0.8,
	merge.confidence.on.shared = 0.8,
	twist.fold.change.mul = 1,
	top.ignored.genes.sender = character(),
	top.ignored.genes.receiver = character(),
	top.n.col.positive = 10,
	top.n.col.negative = 10,
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
			# calculate the score of sender cells
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

	inside.uq.cnt.on.selected.interacts <- function(
		target.pairs
	) {
		tmp.C.df <- data.frame(genes = character(), interacts.cnt = numeric(), interacts.uq.cnt = numeric())
		tmp.D.df <- tmp.C.df
		# function for tapply
		if (nrow(target.pairs) != 0) {
			# calculate the score of sender cells
			tmp.genes.C <- tapply(1:nrow(target.pairs), target.pairs$inter.GeneName.A, length)
			tmp.C.df <- data.frame(genes = names(tmp.genes.C), interacts.cnt = tmp.genes.C, stringsAsFactors = FALSE)
			tmp.C.geneA.itself <- tapply(target.pairs[, "inter.LogFC.A"], target.pairs$inter.GeneName.A, mean)
			tmp.C.df$interacts.uq.cnt <- tmp.C.df$interacts.cnt  * exp(tmp.C.geneA.itself)
			tmp.C.df$interact.marker <- ifelse(tmp.C.geneA.itself > 0, "UPreg", "DNreg")  # marker indicate the A if it is upregulated or not
			# calculate the score of reciever cells
			tmp.genes.D <- tapply(1:nrow(target.pairs), target.pairs$inter.GeneName.B, length)
			tmp.D.df <- data.frame(genes = names(tmp.genes.D), interacts.cnt = tmp.genes.D, stringsAsFactors = FALSE)
			tmp.D.geneB.itself <- tapply(target.pairs[, "inter.LogFC.B"], target.pairs$inter.GeneName.B, mean)
			tmp.D.df$interacts.uq.cnt <- tmp.D.df$interacts.cnt * exp(tmp.D.geneB.itself)
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
			tmp.p.x.shared.inds <- setdiff(seq_len(nrow(this.pair.tgs)), tmp.p.x.diff.inds)
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
	# different strategy
	this.C.spdf <- this.D.spdf <- data.frame()
	if (evaluation.method == "special-score") {
		# get speical genes in diff ones
		p.C.diff.tgs.pnames <- p.D.diff.tgs.pnames <- character()

		# for diff pairs, it is special if it only appears in restricted number of interacting clusters
		p.C.collect.diff.dups <- tapply(seq_along(p.C.collect.pairs$diff), p.C.collect.pairs$diff, length)
		# collect names of gene pairs in diff pairs
		p.C.diff.tgs.pnames <- names(p.C.collect.diff.dups[which(p.C.collect.diff.dups >= floor(merge.confidence.on.diff * length(other.pairs.names.C)))])
		p.D.collect.diff.dups <- tapply(seq_along(p.D.collect.pairs$diff), p.D.collect.pairs$diff, length)
		p.D.diff.tgs.pnames <- names(p.D.collect.diff.dups[which(p.D.collect.diff.dups >= floor(merge.confidence.on.diff * length(other.pairs.names.D)))])
		
		# get speical genes in shared ones
		p.C.shared.tgs.pnames <- p.D.shared.tgs.pnames <- character()
		# for shared pairs, it is special if it is differently expressed against some percentage of interacting clusters that it appears
		p.C.collect.shared.dups <- tapply(seq_along(p.C.collect.pairs$shared), p.C.collect.pairs$shared, length)
		p.C.collect.shared.dups <- p.C.collect.shared.dups[order(names(p.C.collect.shared.dups))]
		# get the shared orig collection
		p.C.collect.shared.orig <- tapply(seq_along(p.C.collect.pairs$shared.orig), p.C.collect.pairs$shared.orig, length)
		p.C.collect.shared.orig <- p.C.collect.shared.orig[which(names(p.C.collect.shared.orig) %in% names(p.C.collect.shared.dups))]
		p.C.collect.shared.orig <- p.C.collect.shared.orig[order(names(p.C.collect.shared.orig))]
		# collect names of gene pairs in shared pairs
		p.C.shared.tgs.pnames <- names(p.C.collect.shared.dups[which(p.C.collect.shared.dups >= floor(merge.confidence.on.shared * p.C.collect.shared.orig))])
		p.D.collect.shared.dups <- tapply(seq_along(p.D.collect.pairs$shared), p.D.collect.pairs$shared, length)
		p.D.collect.shared.dups <- p.D.collect.shared.dups[order(names(p.D.collect.shared.dups))]
		p.D.collect.shared.orig <- tapply(seq_along(p.D.collect.pairs$shared.orig), p.D.collect.pairs$shared.orig, length)
		p.D.collect.shared.orig <- p.D.collect.shared.orig[which(names(p.D.collect.shared.orig) %in% names(p.D.collect.shared.dups))]
		p.D.collect.shared.orig <- p.D.collect.shared.orig[order(names(p.D.collect.shared.orig))]
		p.D.shared.tgs.pnames <- names(p.D.collect.shared.dups[which(p.D.collect.shared.dups >= floor(merge.confidence.on.shared * p.D.collect.shared.orig))])
		
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
	} else {
		if (evaluation.method == "special-index") {
			# get diff index (how diff it is among those interacting pairs)
				# for diff pairs, it is special if it only appears in restricted number of interacting clusters
			p.C.collect.diff.dups <- tapply(seq_along(p.C.collect.pairs$diff), p.C.collect.pairs$diff, length)
			p.D.collect.diff.dups <- tapply(seq_along(p.D.collect.pairs$diff), p.D.collect.pairs$diff, length)
			
			# get speical genes in shared ones
			# for shared pairs, it is special if it is differently expressed against some percentage of interacting clusters that it appears
			p.C.collect.shared.dups <- tapply(seq_along(p.C.collect.pairs$shared), p.C.collect.pairs$shared, length)
			p.C.collect.shared.dups <- p.C.collect.shared.dups[order(names(p.C.collect.shared.dups))]
			p.D.collect.shared.dups <- tapply(seq_along(p.D.collect.pairs$shared), p.D.collect.pairs$shared, length)
			p.D.collect.shared.dups <- p.D.collect.shared.dups[order(names(p.D.collect.shared.dups))]
			
			# merge result
			p.CD.merge.dups <- c(p.C.collect.diff.dups, p.D.collect.diff.dups, p.C.collect.shared.dups, p.D.collect.shared.dups)
			p.CD.tgs <- tapply(p.CD.merge.dups, names(p.CD.merge.dups), sum)
			if (merge.confidence.on.diff == merge.confidence.on.shared) {
				tmp.merge.conf <- merge.confidence.on.diff
			} else {
				tmp.merge.conf <- (merge.confidence.on.diff + merge.confidence.on.shared) / 2
				warning("Under strategy: ", evaluation.method, ". With 2 different *.merge.confidence.*, automatically use the average!")
			}
			p.passed.CD.tgs.pnames <- names(p.CD.tgs[which(p.CD.tgs >= floor(tmp.merge.conf * (length(other.pairs.names.C) + length(other.pairs.names.D))))])
			# get the special pairs
			this.merge.res.df <- this.pair.interacts[match(p.passed.CD.tgs.pnames, this.pair.tgs[, "interacts.name"]), ]

			# calculate the cnt of special pairs in the gene-special manner
			this.merge.res.tl <- inside.uq.cnt.on.selected.interacts(this.merge.res.df)
			this.C.spdf <- this.merge.res.tl$part.C
			this.D.spdf <- this.merge.res.tl$part.D
		} else {
			stop("Error: Undefined strategy! ", evaluation.method)
		}
	}
	
	# ignore some genes in ranking plot
	this.C.spdf <- this.C.spdf[which(this.C.spdf$genes %in% (setdiff(this.C.spdf$genes, top.ignored.genes.sender))), ]
	this.D.spdf <- this.D.spdf[which(this.D.spdf$genes %in% (setdiff(this.D.spdf$genes, top.ignored.genes.receiver))), ]

	if (nrow(this.C.spdf) == 0 || nrow(this.D.spdf) == 0) {
		cat("Quit with no speical genes, in interact [", this.cluster.C, " and ", this.cluster.D, "]!\n", 
			"Re-check if too strict restrictions were given!\n")
		return(NULL)
	}

	#
	inside.topn.col.sort <- function(df, target.col, n.top.sel, decreasing) {
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

	### draw plots
	if (evaluation.method == "special-score") {
		## score plots (consider > 0 and < 0)
		tmp.sym.Y <- sym("interact.score")
		this.C.spdf.col <- inside.topn.col.sort(this.C.spdf, "interact.score", c(top.n.col.positive, top.n.col.negative), TRUE)
		this.D.spdf.col <- inside.topn.col.sort(this.D.spdf, "interact.score", c(top.n.col.positive, top.n.col.negative), TRUE)
	} else {
		## uq.cnt plots
		tmp.sym.Y <- sym("interacts.uq.cnt")
		this.C.spdf.col <- inside.topn.col.sort(this.C.spdf, "interacts.uq.cnt", c(top.n.col.positive, top.n.col.negative), TRUE)
		this.D.spdf.col <- inside.topn.col.sort(this.D.spdf, "interacts.uq.cnt", c(top.n.col.positive, top.n.col.negative), TRUE)
	}
	this.plot.C.col <- ggplot(this.C.spdf.col, aes(x = genes, y = !!tmp.sym.Y))
	this.plot.C.col <- this.plot.C.col + 
			geom_col(aes(fill = interact.marker), colour = "black") + 
			scale_fill_manual(name = "UpDn", values = c("grey", "red"), breaks = c("DNreg", "UPreg")) + 
			scale_x_discrete(breaks = this.C.spdf.col$genes, limits = this.C.spdf.col$genes) + 
			coord_flip() + 
			labs(title = this.cluster.C)

	this.plot.D.col <- ggplot(this.D.spdf.col, aes(x = genes, y = !!tmp.sym.Y))
	this.plot.D.col <- this.plot.D.col + 
			geom_col(aes(fill = interact.marker), colour = "black") + 
			scale_fill_manual(name = "UpDn", values = c("lightgrey", "green"), breaks = c("DNreg", "UPreg")) + 
			scale_x_discrete(breaks = this.D.spdf.col$genes, limits = this.D.spdf.col$genes) + 
			coord_flip() + 
			labs(title = this.cluster.D)
	#
	this.final.4plots <- plot_grid(this.plot.C.col, this.plot.D.col, ncol = 2, align = "vh")

	#
	list(plot = this.final.4plots, tables = list(part.C = this.C.spdf, part.D = this.D.spdf),
		genes.top.on.col = list(part.C = this.C.spdf.col$genes, part.D = this.D.spdf.col$genes), 
		special.pairs.df = this.merge.res.df)
}





#' Find common changed genes in one cluster
#'
#' @description
#' This function finds the common changed genes in one cluster. In other words, 
#' the common changed genes are collected as the intersection of special genes of one cluster in all its 
#' interacting clusters.
#'
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param cluster.name Character. The name of the cluster to be explored.
#' @param cluster.sender.or.receiver Character. It sets the position of the cluster in interacting clusters. 
#' The cluster to be explored is either sender or receiver. For gene pair A-B, the cluster is set as sender, then 
#' all genes in the position A of the gene pairs that the cluster participates will be collected and calculated to get the result.
#' @param common.change.confidence Numeric. It specifies how common the changed genes are. One gene is regarded as common changed 
#' gene if it are the top special genes among several interacting clusters associated with the target cluster specified in the parameter \code{cluster.name}.
#' @param cmp.extend.positive Numeric. It sets the collected range of special genes with positive evaluation results while evaluating by function 
#' \code{FindSpecialGenesInOnepairCluster}. Please refer to that function on the detailed meaning of positive evaluation result.
#' @param cmp.extend.negative Numeric. It sets the collected range of special genes with negative evaluation results while evaluating by function 
#' \code{FindSpecialGenesInOnepairCluster}. Please refer to that function on the detailed meaning of negative evaluation result.
#' @param ... A set of parameters used in \code{FindSpecialGenesInOnepairCluster} will be passed.
#'
#'
#' @details
#' The common changed genes are generated from the special genes of some pairs of interacting clusters. 
#' One gene is regarded as common changed gene if it are the top special genes among several pairs of interacting clusters. 
#'
#'
#'
#' @return Character.
#' Names of common changed genes. 
#'
#'
#'
# NOT export
Testing.FindCommonChangedGenesInOneCluster <- function(
	interact.pairs.acted,
	cluster.name,
	cluster.sender.or.receiver = TRUE,
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
		cluster.position = cluster.sender.or.receiver, full.ref = interact.pairs.acted,
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





#' Draw Circos plot for multiple interaction pairs
#'
#' @description
#' This function draws circos plot for interaction pairs of several target clusters, and returns
#' the detailed result tables about these.
#'
#' @inheritParams Inside.DummyInteractPairsActed
#' @inheritParams Inside.DummyFgenes
#' @inheritParams Inside.DummyActionsRefDB
#' @param return.table.vals Logic. If TRUE, the table of interaction pairs will be given in return value, use \code{<-} to save the result.
#' @param select.interacts Character. Interacts names, e.g. \code{paste0("ClusterC", kClustersSplit, "ClusterD")}, and 
#' the direction will be from \code{"ClusterC"} to \code{"ClusterD"}.
#' @param clusters.select.auto Character. If parameter \code{select.interacts} is not specified, it gives all related clusters,
#' and all permutation of interacts between these clusters will be used.
#' @param is.directional Logic. It is passed to \code{GenerateVEinfos}. If TRUE, it only uses single direction 
#' of clusters' interaction, otherwise use the bi-directional.
#' @param if.ignore.annos Logic. It is passed to \code{GenerateVEinfos}. If TRUE, genes with different locations or types documented will
#' be treated as the same, and only one row information will be reserved.
#' @param sel.mode.val Character. If set NULL, it uses all values in global variables \code{InterCellDB::kpred.mode}, or
#' please specify detailed and accurate values in subset of \code{InterCellDB::kpred.mode}.
#' @param sel.action.effect.val Character. If set NULL, it uses all values in global variables \code{InterCellDB::kpred.action.effect}, or
#' please specify detailed and accurate values in subset of \code{InterCellDB::kpred.action.effect}.
#' @param show.legend Logic. Show legend(TRUE) or not(FALSE).
#' @param cluster.label.cex Numeric. The size of letters when labeling cluster names.
#' @param cluster.colour Character. Predefined some colours in harmony, generated by colorbrew2, 
#' see \url{http://colorbrewer2.org} for details.
#' @param link.cor.colour Character. Colour of links, whose length can be 1 or the same as another parameter \code{pred.mode}.
#' @param link.cor.border.colour Character. Colour of border of links.
#' @param link.cor.arrow.drawn Logic. It decides if the arrows are drawn upon the links. If TRUE, draw the arrows, otherwise not.
#' @param link.cor.border.colour Character. Colour of link borders, whose value length can be 1 or the same as another parameter \code{pred.mode}.
#' @param link.cor.arrow.type Character. Arrow type of links, whose length must be same as another parameter \code{pred.action.effect}.
#' @param link.cor.arrow.length Numeric. Arrow length of links, whose length must be same as another parameter \code{pred.action.effect}.
#'
#' @details
#' This function draws circos plot to show multiple clusters invovled interaction pairs. The final circos plot
#' contains 2 tracks and 1 inner links gragh. (If you are not familiar with those terms, please see \url{} for details.)
#' \itemize{
#'   \item {base track: } {It represents all clusters that are drawn in different colours. 
#'         The colours are specified by parameter \code{cluster.colour}}
#'   \item {expression level track: } {It gives all genes' expression levels by values of "LogFC".
#'         The up-regulated ones are marked in yellowish, and the down-regulated ones are marked in green.}
#'   \item {links: } {It draws links pair-by-pair, and shows mode and action type of every interaction pair.}
#' }
#'
#'
#'
#' @import circlize
#' @importFrom graphics par plot.new
#' @importFrom gridBase gridOMI
#' @import grid
#' @importFrom ComplexHeatmap Legend packLegend
#'
#'
# NOT export
Testing.PlotInteractsInMultipleClusters <- function(
  interact.pairs.acted,
  fgenes.remapped.all,
  actions.ref.db,
  return.table.vals = FALSE,
  select.interacts = NULL,
  clusters.select.auto = NULL,
  is.directional = TRUE,
  if.ignore.annos = TRUE,
  sel.mode.val = NULL,
  sel.action.effect.val = NULL,
  show.legend = FALSE,
  cluster.label.cex = 0.6,
  cluster.colour = c("#fb8072", "#80b1d3", "#fdb462", "#ffffb3", "#bebada", "#b3de69", "#fccde5", "#8dd3c7", "#d9d9d9"),
  link.cor.colour = c("grey"),
  link.cor.border.colour = c("grey"),
  link.cor.arrow.drawn = c(TRUE, TRUE, TRUE, FALSE),
  link.cor.arrow.type = c("triangle", "ellipse", "curved", "circle"),
  link.cor.arrow.length = c(0.4, 0.3, 0.3, 0.2),
  ...
) {
  # pre-process set
  pred.mode <- kpred.mode
  pred.action.effect <- kpred.action.effect
  # pre-process check
  if ((length(link.cor.colour) != 1) && (length(link.cor.colour) < length(pred.mode))) {
    warning("Given insufficent mode-specific colour, use unified colour-'grey' instead.!")
    link.cor.colour <- c("grey")
  }
  # interact pairs selection and alignment of given parameter
  this.sel.interacts.names <- character()
  if (!is.null(select.interacts) && !(length(select.interacts) == 0)) {
    # check if these are valid interact.pairs' names
    logic.ifinnames <- select.interacts %in% interact.pairs.acted$name.allpairs
    if (sum(logic.ifinnames) == length(select.interacts)) {
      this.sel.interacts.names <- select.interacts
    } else {
      stop("Given interacts' names not found: ", paste0(select.interacts[which(logic.ifinnames == FALSE)], collapse = ", "), ".")
    }
  } else {
    if (!is.null(clusters.select.auto) && !((length(clusters.select.auto) == 0))) {
      # check if these are valid clusters' names
      logic.ifclusters <- clusters.select.auto %in% as.character(unique(unlist(interact.pairs.acted$list.clusters)))
      if (sum(logic.ifclusters) == length(logic.ifclusters)) {
        # generating interacts' names
        for (i.cluster in clusters.select.auto) {
          for (j.cluster in clusters.select.auto) {
            if (i.cluster == j.cluster) {  # remove self-interactions
              next
            }
            this.sel.interacts.names <- append(this.sel.interacts.names, paste0(i.cluster, kClustersSplit, j.cluster))
          }
        }
      } else {
        stop("Given clusters' names not found: ", paste0(clusters.select.auto[which(logic.ifclusters == FALSE)], collapse = ", "), ".")
      }
    } else {  # use all interact pairs
      print("Given params are not available, automatically use all clusters to this process!")
      this.sel.interacts.names <- interact.pairs.acted$name.allpairs
    }
  }
  # get data of all interaction pairs as specified
  this.all.interacts <- list()  # get all pairs specified by interacts.names
  this.upd.interacts.names <- character()  # handle with the two same clusters interacts
  if (is.null(sel.mode.val)) {
    sel.mode.val <- pred.mode
  }
  if (is.null(sel.action.effect.val)) {
    sel.action.effect.val <- pred.action.effect
  }
  for (one.interact in this.sel.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.AB.1p <- ExtractTargetOnepairClusters(interact.pairs.acted, tmp.rcluster.A, tmp.rcluster.B)
    tmp.gmoc <- GenerateMapDetailOnepairClusters(tmp.AB.1p, actions.ref.db)
    tmp.VEinfos <- GenerateVEinfos(tmp.gmoc, fgenes.remapped.all, is.directional = is.directional, if.ignore.annos = if.ignore.annos)
    tmp.VEinfos <- TrimVEinfos(tmp.VEinfos, sel.mode.val = sel.mode.val, sel.action.effect.val = sel.action.effect.val)
    if (tmp.rcluster.A == tmp.rcluster.B) {
      tmp.rcluster.B <- paste0(tmp.rcluster.B, ".mirror")
    }
    tmp.new.interact.name <- paste0(tmp.rcluster.A, kClustersSplit, tmp.rcluster.B)
    this.upd.interacts.names <- append(this.upd.interacts.names, tmp.new.interact.name)
    this.all.interacts[[tmp.new.interact.name]] <- tmp.VEinfos
  }  # after GenerateVEinfos, .mirror names will appear
  # get unique clusters related
  this.related.clusters <- unique(unlist(strsplit(this.upd.interacts.names, split = kClustersSplit, fixed = TRUE)))
  ### extract all genes from related clusters
  # create data structure
  this.genes.clusters.based <- list()
  for (i.rcluster in this.related.clusters) {
    this.genes.clusters.based[[i.rcluster]] <- NULL  # data.frame
  }
  # data preparation for circos plot
  for (one.interact in this.upd.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.interact <- this.all.interacts[[one.interact]]
    tmp.vertex.infos <- tmp.interact$vertices.infos
    tmp.vinfos.A <- tmp.vertex.infos[which(tmp.vertex.infos[, "ClusterName"] == tmp.rcluster.A), ]
    this.genes.clusters.based[[tmp.rcluster.A]] <- rbind(this.genes.clusters.based[[tmp.rcluster.A]], tmp.vinfos.A[, c("GeneName", "Location", "LogFC")])
    tmp.vinfos.B <- tmp.vertex.infos[which(tmp.vertex.infos[, "ClusterName"] == tmp.rcluster.B), ]
    this.genes.clusters.based[[tmp.rcluster.B]] <- rbind(this.genes.clusters.based[[tmp.rcluster.B]], tmp.vinfos.B[, c("GeneName", "Location", "LogFC")])
  }
  # doing (1) unique, (2) sorting, (3) give x-axis values
  for (i.r.cluster in this.related.clusters) {
    this.genes.clusters.based[[i.r.cluster]] <- DoPartUnique(this.genes.clusters.based[[i.r.cluster]], c(1,2))
    # here, use sort in alphabet for "Location" to arrange these genes
    tmp.st <- this.genes.clusters.based[[i.r.cluster]]
    this.genes.clusters.based[[i.r.cluster]] <- tmp.st[order(tmp.st[, "Location"]), ]
    tmp.st <- this.genes.clusters.based[[i.r.cluster]]  # give x-axis values
    rownames(tmp.st) <- NULL  # this gives the default rownames, from 1 to nrow(*)
    tmp.st[, "x.val"] <- as.integer(rownames(tmp.st))
    this.genes.clusters.based[[i.r.cluster]] <- tmp.st
  }
  # get packed values
  this.genes.clusters.packed <- NULL  # data.frame
  for (i.r.cluster in this.related.clusters) {
    this.genes.clusters.packed <- rbind(this.genes.clusters.packed, cbind(this.genes.clusters.based[[i.r.cluster]], ClusterName = i.r.cluster))
  }
  #
  ### --- circos plot ---
  # grid the plot
  plot.new()
  kunit.plot = unit(1, "npc")
  if (show.legend) {
    pushViewport(viewport(x = unit(0, "npc"), y = unit(0.5, "npc"), 
      width = 0.8 * kunit.plot, height = kunit.plot, just = c("left", "center")))
  } else {
    pushViewport(viewport(x = unit(0, "npc"), y = unit(0.5, "npc"), 
      width = kunit.plot, height = kunit.plot, just = c("left", "center")))
  }
  # grid the region
  par(omi = gridOMI())
  par(new = TRUE)
  # extra clear to avoid unexpected program errors
  circos.clear()
  ## draw the plot
  circos.par("track.height" = 0.1, "gap.degree" = 3)
  circ.factors <- factor(this.genes.clusters.packed[, "ClusterName"])
  circ.xlim <- NULL
  for (i.rcluster in levels(circ.factors)) {
    tmp.this.xval <- this.genes.clusters.based[[i.rcluster]][, "x.val"]
    circ.xlim <- rbind(circ.xlim, c(min(tmp.this.xval) - 1, max(tmp.this.xval) + 1))
  }
  # initialize the circos
  circos.initialize(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], xlim = circ.xlim)
  # add base track
  kcolset.basetrack <- cluster.colour
  circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], ylim = c(0, 1),
    track.height = 0.05,
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), CELL_META$sector.index)
      circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
        col = kcolset.basetrack[(CELL_META$sector.numeric.index-1) %% length(kcolset.basetrack)+1])
  })
  # add gene name track
  circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], ylim = c(0, 32),
    track.height = 0.20, cell.padding = c(0, 0, 0, 0),
    bg.col = NA, bg.border = NA,
    panel.fun = function(x, y) {
      labels.genes <- this.genes.clusters.packed[which(this.genes.clusters.packed[, "ClusterName"] == CELL_META$sector.index), "GeneName"]
      circos.text(x, CELL_META$ycenter, labels.genes, 
        facing = "clockwise", niceFacing = TRUE,
        cex = cluster.label.cex)
  })
  # add or logfc
  kcolor.mval <- "#dcdcdc"
  if (min(this.genes.clusters.packed[, "LogFC"] > 0)) {
    circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], y = this.genes.clusters.packed[, "LogFC"],
      panel.fun = function(x, y) {
        circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = kcolor.mval)
        circos.rect(x - 0.5, CELL_META$cell.ylim[1], x + 0.5, y, col = "#efe91a")
    })
  } else {
    circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], y = this.genes.clusters.packed[, "LogFC"],
      panel.fun = function(x, y) {
        circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = kcolor.mval)
        circos.rect(x - 0.5, 0, x + 0.5, y, col = sapply(y, function(y) {if (y < 0) "green" else "#efe91a"}))
    })
  }
  ## draw links for interaction pairs
  thisf.used.mode <- NULL  # collect all used mode
  thisf.used.action.effect <- NULL  # collect all used action effects
  for (one.interact in this.upd.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.this.data <- this.all.interacts[[one.interact]]
    tmp.edges <- tmp.this.data$edges.infos
    tmp.vertices <- tmp.this.data$vertices.infos
    # get all links pairs with c(x1.val, x2.val) & (mode, action.effect)
    tmp.links.xvals <- apply(tmp.edges, MARGIN = 1, 
      vertices = tmp.vertices, gcb.ref = this.genes.clusters.based,
      FUN = function(x, vertices, gcb.ref) {
        # for x1
        uid.x1 <- as.integer(x["from"])  # get UID in this interaction pair
        ind.mv.x1 <- which(vertices[, "UID"] == uid.x1)  # get correspond row in vertices.infos
        tmp.cluster.x1 <- vertices[ind.mv.x1, "ClusterName"]
        tgcb.x1 <- gcb.ref[[tmp.cluster.x1]]  # get this cluster infos in gcb.ref
        ind.mgcb.x1 <- which(tgcb.x1[, "GeneName"] == vertices[ind.mv.x1, "GeneName"])  # find the same gene in gcb.ref[[this.cluster]]
        ind.mgcb.loc.x1 <- which(tgcb.x1[, "Location"] == vertices[ind.mv.x1, "Location"])
        val.x1 <- tgcb.x1[intersect(ind.mgcb.x1, ind.mgcb.loc.x1), "x.val"]  # get x.val in unique 1 row
        # for x2
        uid.x2 <- as.integer(x["to"])
        ind.mv.x2 <- which(vertices[, "UID"] == uid.x2)
        tmp.cluster.x2 <- vertices[ind.mv.x2, "ClusterName"]
        tgcb.x2 <- gcb.ref[[tmp.cluster.x2]]
        ind.mgcb.x2 <- which(tgcb.x2[, "GeneName"] == vertices[ind.mv.x2, "GeneName"])
        ind.mgcb.loc.x2 <- which(tgcb.x2[, "Location"] == vertices[ind.mv.x2, "Location"])
        val.x2 <- tgcb.x2[intersect(ind.mgcb.x2, ind.mgcb.loc.x2), "x.val"]  # get x.val in unique 1 row
        # other
        save.mode <- x["mode"]
        save.action.effect <- x["action.effect"]
        # return
        as.character(c(tmp.cluster.x1, val.x1, tmp.cluster.x2, val.x2, save.mode, save.action.effect))
    })
    tmp.links.xvals <- t(tmp.links.xvals)
    if (is.null(nrow(tmp.links.xvals)) || (nrow(tmp.links.xvals) == 0)) {  # go to next loop
      next
    }
    for (i in 1:nrow(tmp.links.xvals)) {
      # map different mode to different colour
      tmp.link.colour <- link.cor.colour  # assume using unified colour
      tmp.link.border <- link.cor.border.colour  # assume using unified colour
      if (length(link.cor.colour) > 1) {
        tmp.link.colour <- link.cor.colour[which(tmp.links.xvals[i, 5] == pred.mode)]
        tmp.link.border <- link.cor.border.colour[which(tmp.links.xvals[i, 5] == pred.mode)]
        thisf.used.mode <- c(thisf.used.mode, tmp.links.xvals[i, 5])
      }
      # arrow type
      thisf.used.action.effect <- c(thisf.used.action.effect, tmp.links.xvals[i, 6])
      ind.tmp.link.type <- which(tmp.links.xvals[i, 6] == pred.action.effect)
      tmp.this.arr.drawn <- link.cor.arrow.drawn[ind.tmp.link.type]
      tmp.this.arr.type <- link.cor.arrow.type[ind.tmp.link.type]
      tmp.this.arr.length <- link.cor.arrow.length[ind.tmp.link.type]
      if (tmp.this.arr.drawn == TRUE) {
        circos.link(tmp.links.xvals[i, 1], c(as.numeric(tmp.links.xvals[i, 2]) - 0.2, as.numeric(tmp.links.xvals[i, 2]) + 0.2), tmp.links.xvals[i, 3], as.numeric(tmp.links.xvals[i, 4]),
            col = tmp.link.colour, border = tmp.link.border, 
            directional = 1, arr.length = tmp.this.arr.length, arr.type = tmp.this.arr.type)
      } else {
        circos.link(tmp.links.xvals[i, 1], c(as.numeric(tmp.links.xvals[i, 2]) - 0.2, as.numeric(tmp.links.xvals[i, 2]) + 0.2), 
            tmp.links.xvals[i, 3], c(as.numeric(tmp.links.xvals[i, 4]) - 0.2, as.numeric(tmp.links.xvals[i, 4]) + 0.2),
            col = tmp.link.colour, border = tmp.link.border, 
            directional = 0, arr.length = tmp.this.arr.length, arr.type = tmp.this.arr.type)
      }
    }
  }
  circos.clear()
  # grid work, move to parent grid frame
  upViewport()
  ## draw the legend
  # legend exprs track
  leg.logfc.max <- max(this.genes.clusters.packed[, "LogFC"])
  leg.logfc.min <- min(this.genes.clusters.packed[, "LogFC"])
  exprs.p.at <- c("down-regulation", "up-regulation")
  exprs.p.col <- c("green", "#efe91a")
  if (leg.logfc.max * leg.logfc.min < 0) {
    # keep default values
  } else {
    if (leg.logfc.min > 0) {
      exprs.p.at <- c("up-regulation")
      exprs.p.col <- c("#efe91a")
    } else {
      exprs.p.at <- c("down-regulation")
      exprs.p.col <- c("green")
    }
  }
  lgd.exprs <- Legend(at = exprs.p.at, type = "points", 
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(col = exprs.p.col),
            title_position = "topleft", title = paste0("Track2", ":LogFC"))
  # legend arrow type
  thisf.used.action.effect <- levels(factor(thisf.used.action.effect))
  inds.lgd.arr.type <- match(thisf.used.action.effect, pred.action.effect)
  thisf.used.arr.type <- link.cor.arrow.type[inds.lgd.arr.type]
  lgd.links.arr.type <- Legend(at = c(toupper(thisf.used.arr.type), thisf.used.action.effect), ncol = 2,
            grid_width = unit(3, "mm"),
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(fill = c(rep("#ffffff", times = length(thisf.used.arr.type)), rep("#eeeeee", times = length(thisf.used.action.effect)))),
            title_position = "topleft", title = "Links-ActionType")
  # packed legend if no link colour is used,
  lgd.packed.list <- packLegend(lgd.exprs, lgd.links.arr.type)
  # conditionally, legend links colour
  if (length(link.cor.colour) > 1) {  # legend for link colour is required
    thisf.used.mode <- levels(factor(thisf.used.mode))
    inds.lgd.links.colour <- match(thisf.used.mode, pred.mode)
    thisf.used.links.colour <- link.cor.colour[inds.lgd.links.colour]
    lgd.links.colour <- Legend(at = thisf.used.mode, type = "lines",
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(col = thisf.used.links.colour, lwd = 2),
            title_position = "topleft", title = paste0("Links-Mode"))
    lgd.packed.list <- packLegend(lgd.exprs, lgd.links.arr.type, lgd.links.colour)
  }
  pushViewport(viewport(x = 0.8 * kunit.plot, y = unit(0.5, "npc"), 
    width = 0.12 * kunit.plot, height = kunit.plot, just = c("left", "center")))
  if (show.legend) {
    grid.draw(lgd.packed.list)
  }
  # grid work
  upViewport()
  # end plotting

  # get result tables
  if (return.table.vals) {
    this.result.table.list <- list()
    for (i in 1:length(this.all.interacts)) {
      tmp.interact.name <- this.upd.interacts.names[i]
      tmp.it <- this.all.interacts[[i]]
      tmp.vertices <- tmp.it$vertices.infos
      tmp.edges <- tmp.it$edges.infos
      inds.tmp.from.match <- match(tmp.edges$from, tmp.vertices$UID)
      inds.tmp.to.match <- match(tmp.edges$to, tmp.vertices$UID)
      tmp.edges[, c("from.ClusterName", "from.GeneName")] <- tmp.vertices[inds.tmp.from.match, c("ClusterName", "GeneName")]
      tmp.edges[, c("to.ClusterName", "to.GeneName")] <- tmp.vertices[inds.tmp.to.match, c("ClusterName", "GeneName")]
      tmp.edges <- tmp.edges[, c("from", "from.ClusterName", "from.GeneName", "to", "to.ClusterName", "to.GeneName", "mode", "action.effect")]
      tmp.vertices <- tmp.vertices[, c("UID", "ClusterName", "GeneName", "Location", "LogFC")]
      tmp.list <- list(tmp.vertices, tmp.edges)
      names(tmp.list) <- c(paste0(tmp.interact.name, "-vertices"), paste0(tmp.interact.name, "-edges"))
      this.result.table.list <- c(this.result.table.list, tmp.list)
    }
    #end#
    return(list(table = this.result.table.list))
  } else {
    return(NULL)  # no return value
  }
}
