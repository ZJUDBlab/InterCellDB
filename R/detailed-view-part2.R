


#' Find special genes in one pair of interacting clusters
#'
#' @description
#' This function is used to find special genes in one pair of interacting clusters. Genes are special
#' if it passes some limitations when comparing to other pairs of interacting clusters. 
#'
#' @inheritParams Inside.DummyVEinfos
#' @param interact.pairs.acted List. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param to.cmp.clusters.pairs [TODO]
#' @param to.cmp.clusters.pairs.sel.strategy [TODO]
#' @param uq.cnt.range [TODO]
#' @param formula.to.use.onLogFC [TODO]
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
	VEinfos,
	interact.pairs.acted,
	to.cmp.clusters.pairs = character(),  # should be a subset of interact.pairs.acted$name.allpairs
	to.cmp.clusters.pairs.sel.strategy = "inter-cross",
	uq.cnt.range = c(1:2),  # [TODO] GetResult.* futher select once more. 100 clusters as desired default maximum
	formula.to.use.onLogFC = Tool.formula.onLogFC.default
) {
	# pre-check
	if (class(formula.to.use.onLogFC) != "function") {
		stop("Please provide usable function in parameter formula.to.use.onLogFC!")
	}
	# pre-params
	this.gp.sep <- kGenesSplit
	involved.clusters.pair <- paste(VEinfos$cluster.name.A, VEinfos$cluster.name.B, sep = kClustersSplit)
	target.gene.pairs.df <- Tool.GenStdGenePairs.from.VEinfos(VEinfos)

	# calculate the multiply of fold change of interacting genes
	inside.c.short.interacts <- function(
		tmp.pairs,
		tmp.sep
	) {
		paste(tmp.pairs[, "inter.GeneName.A"], tmp.pairs[, "inter.GeneName.B"], sep = tmp.sep)
	}
	# calculate the logfc
	inside.calc.upon.gp.logfc <- function(
		tmp.pairs,
		tmp.formula,
		tmp.colnames = c("inter.LogFC.A", "inter.LogFC.B")
	) {
		tmp.formula(tmp.pairs[, tmp.colnames[1]], tmp.pairs[, tmp.colnames[2]])
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
	# this.pair.logfc.mul <- inside.calc.upon.gp.logfc(target.gene.pairs.df, formula.to.use.onLogFC)
	this.gene.pairs <- data.frame(gp.name = this.pair.interacts, 
		gp.belongs = rep(this.pair.name, times = length(this.pair.interacts)), 
		gp.logfc.A = target.gene.pairs.df[, "inter.LogFC.A"], 
		gp.logfc.B = target.gene.pairs.df[, "inter.LogFC.B"], 
		stringsAsFactors = FALSE)

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

	# set the uq.cnt range
	tmp.all.uq.clusters <- unique(as.character(unlist(strsplit(to.cmp.clusters.pairs, split = kClustersSplit, fixed = TRUE))))
	uq.cnt.range <- uq.cnt.range[which(uq.cnt.range %in% seq_along(tmp.all.uq.clusters))]
	if (length(uq.cnt.range) == 0) {
		stop("Please specify available uq.cnt.range to continue the analysis!\nThe allowed values are ", 
			paste0(seq_along(tmp.all.uq.clusters), collapse = ", "), ".")
	}

	# generate gene pairs to compare
	  # only those pairs occurs in target cluster group will be further analyzed
	other.gene.pairs.df.list <- lapply(to.cmp.clusters.pairs, 
		all.pairs.interacts = all.pairs.interacts, tg.pairs = this.pair.interacts,
		tmp.sep = this.gp.sep, formula.to.use.onLogFC = formula.to.use.onLogFC, 
		function(x, all.pairs.interacts, tg.pairs, tmp.sep, formula.to.use.onLogFC) {
			tmp.pairs.df <- all.pairs.interacts[[x]]
			tmp.pairs.interacts <- inside.c.short.interacts(tmp.pairs.df, tmp.sep)
			tmp.inds.tg <- which(tmp.pairs.interacts %in% tg.pairs)
			res.df <- data.frame(gp.name = tmp.pairs.interacts[tmp.inds.tg], 
				gp.belongs = rep(x, times = length(tmp.inds.tg)), 
				gp.logfc.A = tmp.pairs.df[tmp.inds.tg, "inter.LogFC.A"], 
				gp.logfc.B = tmp.pairs.df[tmp.inds.tg, "inter.LogFC.B"], 
				stringsAsFactors = FALSE)
			res.df
		})
	other.gene.pairs.packed <- bind_rows(other.gene.pairs.df.list)
	
	# get unique counts result
	this.all.gp.packed <- rbind(this.gene.pairs, other.gene.pairs.packed)
	this.uq.cnts <- tapply(this.all.gp.packed$gp.belongs, this.all.gp.packed$gp.name, length)
	this.uq.indices <- tapply(seq_len(nrow(this.all.gp.packed)), this.all.gp.packed$gp.name, function(x) { x })
	# select 

	## select diff.cnt subset
	# uq.cnt
	tmp.inds.uq.valid <- this.uq.cnts %in% uq.cnt.range
	this.uq.cnts <- this.uq.cnts[tmp.inds.uq.valid]
	this.uq.indices <- this.uq.indices[tmp.inds.uq.valid]
	# [TODO] further subset options

	# collect by gene pairs each
	tmp.diff.res.mat <- vector(mode = "list", length = length(this.uq.indices))
	prog.bar.use.plot.collect <- progress::progress_bar$new(total = length(this.uq.indices))
	prog.bar.use.plot.collect$tick(0)
	tmp.diff.res.mat <- lapply(names(this.uq.indices), all.gp = this.all.gp.packed, 
	diff.cnts = this.uq.cnts, diff.indices = this.uq.indices, formula.to.use = formula.to.use.onLogFC, 
		function(x, all.gp, diff.cnts, diff.indices, formula.to.use) {
			tmp.i <- which(names(diff.cnts) == x)  # must be the same as it is in diff.indices
			tmp.cnt <- diff.cnts[tmp.i]
			names(tmp.cnt) <- NULL  # remove the name
			tmp.indices <- diff.indices[[tmp.i]]
			tmp.properties <- all.gp[tmp.indices, c("gp.belongs", "gp.logfc.A", "gp.logfc.B")]
			tmp.properties$gp.logfc.calc <- inside.calc.upon.gp.logfc(tmp.properties, formula.to.use, tmp.colnames = c("gp.logfc.A", "gp.logfc.B"))
			tmp.properties <- tmp.properties[order(tmp.properties$gp.logfc.calc, decreasing = TRUE), ]
			prog.bar.use.plot.collect$tick()  # tick
			# return
			list(uq.cnt = tmp.cnt, uq.details = tmp.properties)
		})
	names(tmp.diff.res.mat) <- names(this.uq.indices)

	# collect gene pairs in one data.frame
	tmp.diff.df.list <- vector(mode = "list", length = length(this.uq.indices))
	prog.bar.for.res.collect <- progress::progress_bar$new(total = length(this.uq.indices))
	prog.bar.for.res.collect$tick(0)
	tmp.diff.df.list <- lapply(names(this.uq.indices), all.gp = this.all.gp.packed, diff.indices = this.uq.indices, 
		function(x, all.gp, diff.indices) {
			tmp.i <- which(names(diff.indices) == x)  # must be the same as it is in diff.indices
			tmp.indices <- diff.indices[[tmp.i]]
			tmp.genepair <- names(diff.indices[tmp.i])  # get the gene pair name
			tmp.gene.part <- strsplit(tmp.genepair, split = kGenesSplit, fixed = TRUE)[[1]]
			tmp.properties <- all.gp[tmp.indices, c("gp.belongs", "gp.logfc.A", "gp.logfc.B")]
			tmp.clusters.prop.df <- Tool.SplitToGenDataFrame(tmp.properties[, "gp.belongs"], to.split.by = kClustersSplit, 
				res.colnames = c("Cluster.G1", "Cluster.G2"))
			tmp.res.df <- cbind(data.frame(gp.name = rep(tmp.genepair, times = nrow(tmp.properties)), 
					inter.GeneName.A = rep(tmp.gene.part[1], times = nrow(tmp.properties)),
					inter.GeneName.B = rep(tmp.gene.part[2], times = nrow(tmp.properties)),
					stringsAsFactors = FALSE), 
				tmp.properties$gp.belongs, tmp.clusters.prop.df, tmp.properties[, c("gp.logfc.A", "gp.logfc.B")], stringsAsFactors = FALSE)
			prog.bar.for.res.collect$tick()
			# result
			tmp.res.df
		})
	tmp.diff.df.res <- bind_rows(tmp.diff.df.list)

	return(list(res = tmp.diff.df.res, for.plot.use = tmp.diff.res.mat))
}





# GetResult.SummarySpecialGenes Select Gene Pairs Method: random
Inside.select.genepairs.method.random.default <- function(
	this.spgenes, 
	select.by.method.pairs.limit,
	...
) {
	tmp.inds.sel <- sample(seq_along(this.spgenes), select.by.method.pairs.limit)
	select.genepairs <- names(this.spgenes)[tmp.inds.sel]
	return(select.genepairs)
}
Inside.select.genepairs.method.random.IT <- Inside.select.genepairs.method.random.default

# GetResult.SummarySpecialGenes Select Gene Pairs Method: LogFC sum(decreasing or increasing)
Inside.select.genepairs.method.logfc.sum.default <- function(
	this.spgenes, 
	VEinfos, 
	select.by.method.pairs.limit, 
	select.by.method.decreasing = TRUE, 
	...
) {
	tmp.cluster.A <- VEinfos$cluster.name.A
	tmp.cluster.B <- VEinfos$cluster.name.B
	tmp.gp.name <- paste(tmp.cluster.A, tmp.cluster.B, sep = kClustersSplit)
	tmp.tg.calc <- unlist(lapply(seq_along(this.spgenes), 
		spgenes = this.spgenes, tg.gp.name = tmp.gp.name, 
		function(x, spgenes, tg.gp.name) {
		tmp.df <- spgenes[[x]]$uq.details
		tmp.df[which(tmp.df$gp.belongs == tg.gp.name), "gp.logfc.calc"]
		}))
	# max -> min or min -> max
	tmp.inds.sel <- order(tmp.tg.calc, decreasing = select.by.method.decreasing)
	length(tmp.inds.sel) <- select.by.method.pairs.limit
	tmp.inds.sel <- tmp.inds.sel[which(!is.na(tmp.inds.sel))]
	select.genepairs <- names(this.spgenes)[tmp.inds.sel]
	return(select.genepairs)
}
Inside.select.genepairs.method.logfc.sum.IT <- Inside.select.genepairs.method.logfc.sum.default

Inside.select.genepairs.method.diff.logfc.sum.default <- function(
	this.spgenes, 
	VEinfos, 
	select.by.method.pairs.limit, 
	select.by.method.decreasing = TRUE, 
	...
) {
	# [TODO]
}
Inside.select.genepairs.method.diff.logfc.sum.IT <- Inside.select.genepairs.method.diff.logfc.sum.default



#' Summary special genes in cluster groups
#'
#' @description
#' [TODO]
#'
#' @param onepair.spgenes [TODO]
#' @inheritParams Inside.DummyVEinfos
#' @param show.uq.cnt.range [TODO]
#' @param show.uq.cnt.merged [TODO]
#' @param select.genepairs [TODO]
#' @param select.genepairs.method [TODO]
#' @param select.by.method.pairs.limit [TODO]
#' @param show.topn.inside.onepair [TODO]
#' @param prioritize.cluster.groups [TODO]
#' @param plot.font.size.base [TODO]
#' @param facet.scales [TODO]
#' @param facet.space [TODO]
#' @param facet.text.x [TODO]
#' @param facet.background [TODO]
#' @param bar.colour [TODO]
#' @param bar.width [TODO]
#' @param axis.text.x.pattern [TODO]
#'
#' @details
#' [TODO]
#'
#'
#' @return List. [TODO]
#'
#'
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import cowplot
#'
#' 
#' @export
#'
GetResult.SummarySpecialGenes <- function(
	onepair.spgenes,
	VEinfos, 
	show.uq.cnt.range = c(1), 
	show.uq.cnt.merged = TRUE,  # merged shows different uq.cnt in one plot, or in several plots
	select.genepairs = character(), 
	select.genepairs.method = "random",  # logfc-sum, diff-logfc-sum, are other 2 options pre-defined
	select.by.method.pairs.limit = 10, 
	select.by.method.decreasing = TRUE,     # used in logfc-sum, diff-logfc-sum
	select.by.method.diff.option = "1-mean",  # only used in diff-logfc-sum, mean-mean, 1-1, 1-mean, mean-1
	show.topn.inside.onepair = 2,
	prioritize.cluster.groups = character(),  # names put here will be showed in left and order as it is in here 
	# plotting param
	plot.font.size.base = 12, 
	facet.scales = "free_x", 
	facet.space  = "free_x", 
	facet.text.x = element_text(size = 8, colour = "black"), 
	facet.background = element_rect(fill = "lightgrey", colour = "white"), 
	bar.colour = character(), 
	bar.width = 0.8, 
	axis.text.x.pattern = element_text(angle = 90, vjust = 0.5, hjust = 1),
	...
) {
	this.spgenes <- onepair.spgenes$for.plot.use
	## pre-check
	this.property.valid.uq.cnt <- unique(as.integer(unlist(lapply(this.spgenes, function(x) { x$uq.cnt }))))
	# check valid uq.cnt
	tmp.inds.valid.uq <- which(show.uq.cnt.range %in% this.property.valid.uq.cnt)
	if (length(tmp.inds.valid.uq) != length(show.uq.cnt.range)) {
		warning("Select not-existed uq.cnt: ", paste0(show.uq.cnt.range[setdiff(seq_along(show.uq.cnt.range), tmp.inds.valid.uq)], collapse = ", "), 
			", which will be ignored!")
	}
	show.uq.cnt.range <- show.uq.cnt.range[tmp.inds.valid.uq]
	# check topn
	if (show.topn.inside.onepair < 0) {
		stop("Top N must be larger than 0! Check parameter show.uq.cnt.range!")
	}
	# check cluster group given
	this.property.valid.cluster.group <- unique(as.character(unlist(lapply(this.spgenes, function(x) { x$uq.details$gp.belongs }))))
	tmp.inds.valid.cluster.group <- which(prioritize.cluster.groups %in% this.property.valid.cluster.group)
	if (length(tmp.inds.valid.cluster.group) != length(prioritize.cluster.groups)) {
		warning("Given cluster group order has some items not existed: ",
			paste0(prioritize.cluster.groups[setdiff(seq_along(prioritize.cluster.groups), tmp.inds.valid.cluster.group)]),
			", which will be automatically removed!")
	}

	# check select gene pairs
	if (length(select.genepairs) == 0) {  # then use method
		# check maximum pairs limit
		if (select.by.method.pairs.limit > length(this.spgenes)) {
			warning("Maximum gene pairs are:", length(this.spgenes), ", and given limit is automatically shrinked to that value.")
			select.by.method.pairs.limit <- length(this.spgenes)
		}
		# select by methods
		select.genepairs <- switch(select.genepairs.method, 
			"random" = Inside.select.genepairs.method.random.IT(this.spgenes, select.by.method.pairs.limit, ...), 
			"logfc-sum" = Inside.select.genepairs.method.logfc.sum.IT(this.spgenes, VEinfos, select.by.method.pairs.limit, select.by.method.decreasing, ...), 
			stop("Undefined method for selecting gene pairs. Please re-check the given param: select.genepairs.method!!")
		)
	} else {  # check validity of those select gene pairs
		tmp.inds.vd.gp <- which(select.genepairs %in% names(this.spgenes))
		if (length(tmp.inds.vd.gp) != length(select.genepairs)) {
			warning("Selected gene pairs not exist: ", 
				paste0(select.genepairs[setdiff(seq_along(select.genepairs), tmp.inds.vd.gp)], collapse = ", "), 
				", which will be removed from following analysis!")
		}
		select.genepairs <- select.genepairs[tmp.inds.vd.gp]
	}
	# unique on gene pairs
	select.genepairs <- unique(select.genepairs)

	# template function for dealing every valid uq.cnt
	inside.collect.uq.cnt.each <- function(
		all.spgenes,
		uq.cnt.it,  # target uq.cnt
		topn.it
	) {
		tmp.inds <- unlist(lapply(all.spgenes, uq.cnt.it = uq.cnt.it, function(x, uq.cnt.it) {
			x$uq.cnt == uq.cnt.it
		}))
		tmp.spgenes <- all.spgenes[tmp.inds]
		tmp.spgenes <- lapply(tmp.spgenes, topn = topn.it, function(x, topn) {
			tmp.selrows <- ifelse(nrow(x$uq.details) > topn, topn, nrow(x$uq.details))
			list(uq.cnt = x$uq.cnt, uq.details = x$uq.details[seq_len(tmp.selrows), c("gp.belongs", "gp.logfc.calc")])
			})
		# remove those with NO valid uq.cnt 
		tmp.check.0row <- unlist(lapply(tmp.spgenes, function(x) { nrow(x$uq.details) }))
		tmp.spgenes <- tmp.spgenes[which(tmp.check.0row > 0)]
		return(tmp.spgenes)
	}

	inside.trans.coords.uq <- function(
		std.spgenes,
		show.genepairs.order,  # ordered as it is in select.genepairs
		prioritize.cluster.groups,
		std.width.col = 2,  # may be export as param, so as the gap
		std.width.gap = 3
	) {
		#[NOTE]# x-axis coords go from 0 -> +Inf, y-axis use the original value
		tmp.order.df <- data.frame(orig.ind = seq_along(std.spgenes), 
			new.ind = match(names(std.spgenes), show.genepairs.order), 
			stringsAsFactors = FALSE)
		tmp.order.df <- tmp.order.df[order(tmp.order.df$new.ind, decreasing = FALSE), ]
		std.spgenes <- std.spgenes[tmp.order.df$orig.ind]

		# calculate the gap - every col width 2, gap width 3
		tmp.df.list <- lapply(seq_along(std.spgenes), std.spgenes = std.spgenes, 
			std.width.col = std.width.col, std.width.gap = std.width.gap, 
			prioritize.cluster.groups = prioritize.cluster.groups, 
			function(x, std.spgenes, std.width.col, std.width.gap, prioritize.cluster.groups) {
				tmp.begin <- x * std.width.gap + (x - 1) * std.width.col
				tmp.n.items <- std.spgenes[[x]]$uq.cnt
				tmp.coords.seq <- (seq_len(tmp.n.items) - 1) * std.width.col + tmp.begin
				tmp.ref.rows <- nrow(std.spgenes[[x]]$uq.details)
				# get cluster group order. Matched ones will be put in front all other.
				tmp.belongs <- std.spgenes[[x]]$uq.details$gp.belongs
				tmp.inds.prior <- which(tmp.belongs %in% prioritize.cluster.groups)
				tmp.inds.inferior <- setdiff(seq_along(tmp.belongs), tmp.inds.prior)
				# result
				data.frame(uq.name = rep(names(std.spgenes)[x], times = tmp.ref.rows), 
					uq.label = tmp.belongs[c(tmp.inds.prior, tmp.inds.inferior)], 
					uq.cnt = rep(tmp.n.items, times = tmp.ref.rows), 
					uq.x.axis = tmp.coords.seq[seq_len(tmp.ref.rows)], 
					uq.y.axis = std.spgenes[[x]]$uq.details$gp.logfc.calc[c(tmp.inds.prior, tmp.inds.inferior)],
					stringsAsFactors = FALSE)
				})
		tmp.df.res <- bind_rows(tmp.df.list)
		return(tmp.df.res)
	}

	## process: reconstruct data format of every uq.cnt to fit plotting
	if (show.uq.cnt.merged == TRUE) {
		plot.data.uq <- list()
		for (iuq in show.uq.cnt.range) {
			tmp.spgenes <- inside.collect.uq.cnt.each(this.spgenes, iuq, show.topn.inside.onepair)
			tmp.spgenes <- tmp.spgenes[which(names(tmp.spgenes) %in% select.genepairs)]
			plot.data.uq <- c(plot.data.uq, tmp.spgenes)
		}
		plot.data.uq.df <- inside.trans.coords.uq(plot.data.uq, select.genepairs, prioritize.cluster.groups)
	} else {
		plot.data.uq.notm.list <- list()
		tmp.uq.cnt.valid.list <- integer()
		for (iuq in show.uq.cnt.range) {
			tmp.spgenes <- inside.collect.uq.cnt.each(this.spgenes, iuq, show.topn.inside.onepair)
			tmp.spgenes <- tmp.spgenes[which(names(tmp.spgenes) %in% select.genepairs)]
			if (length(tmp.spgenes) == 0) {
				next  # not added to the result list
			}
			tmp.spgenes.trans <- inside.trans.coords.uq(tmp.spgenes, select.genepairs, prioritize.cluster.groups)
			plot.data.uq.notm.list <- c(plot.data.uq.notm.list, list(tmp.spgenes.trans))
			tmp.uq.cnt.valid.list <- c(tmp.uq.cnt.valid.list, iuq)
		}
		names(plot.data.uq.notm.list) <- paste("uq.cnt=", as.character(tmp.uq.cnt.valid.list), sep = ".")
	}

	inside.uq.single.plot <- function(
		plot.data,
		prioritize.cluster.groups,
		bar.colour
	) {
		# generate template 20 colours to fit most circumstances
		colours.group <- brewer.pal.info[brewer.pal.info$category == 'seq', ]
		colour.sel.author.prefer <- unlist(mapply(brewer.pal, colours.group$maxcolors, rownames(colours.group)))
		colour.sel.author.prefer <- colour.sel.author.prefer[c(9:3, 39:45, 23:26, 51, 53)]  # selected by author's well
		# colour align with the prioritize.cluster.groups
		tiny.cg.prior <- function(x) {
			tmp.inds.prior <- which(x %in% prioritize.cluster.groups)
			tmp.inds.inferior <- setdiff(seq_along(x), tmp.inds.prior)
			x[c(tmp.inds.prior, tmp.inds.inferior)]
		}
		colour.sel.cor.breaks <- tiny.cg.prior(levels(factor(plot.data$uq.label)))
		if (is.null(bar.colour) || length(bar.colour) == 0) {
			colour.sel.it <- colour.sel.author.prefer
		} else {  # use user defined colours
			colour.sel.it <- bar.colour
			if (length(bar.colour) < length(colour.sel.cor.breaks)) {
				warning("Given kinds of colours are not enough to cover all cluster groups! The program automatically fills some colours.")
				for (try.i in 1:100) {  # try several times to fill needed number of colours
					colour.sel.it <- c(colour.sel.it, colour.sel.author.prefer)
					if (length(colour.sel.it) >= length(colour.sel.cor.breaks)) {
						break
					}
					if (try.i == 100) {
						stop("Cannot give as many colours as needed, program failed!")
					}
				}
			}
		}
		colour.sel.it <- colour.sel.it[seq_along(colour.sel.cor.breaks)]
		names(colour.sel.it) <- colour.sel.cor.breaks
		# before plot, get facet levels correctly arranged
		plot.data$uq.name <- factor(plot.data$uq.name)
		levels(plot.data$uq.name) <- unique(plot.data$uq.name)
		# the plot
		this.plot <- ggplot(plot.data, aes(x = uq.label, y = uq.y.axis))
		this.plot <- this.plot + 
			geom_col(aes(fill = uq.label), width = bar.width) + 
			facet_grid(cols = vars(uq.name), scales = facet.scales, space = facet.space) + 
			scale_x_discrete(breaks = plot.data$uq.label, 
				limits = tiny.cg.prior,   # function to change the limit
				labels = plot.data$uq.label) + 
			scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # keep the default 5% in y top
			scale_fill_manual(name = "Cluster Group", 
				values = colour.sel.it, 
				breaks = colour.sel.cor.breaks) + 
			labs(x = "Cluster Groups", y = "LogFC Calc")
		this.plot <- this.plot + 
			theme_half_open(font_size = plot.font.size.base) + 
			theme(axis.text.x = axis.text.x.pattern,
				strip.text.x = facet.text.x,
				strip.background = facet.background)
		return(this.plot)
	}

	inside.gen.ret.table <- function(
		ret.data
	) {
		ret.res <- ret.data[, c("uq.name", "uq.label", "uq.cnt", "uq.y.axis")]
		colnames(ret.res) <- c("gene.pairs", "cluster.groups", "uq.cnt", "gp.logfc.calc")
		ret1.genepairs.df <- Tool.SplitToGenDataFrame(ret.res[, "gene.pairs"], 
			to.split.by = kGenesSplit, 
			res.colnames = c("inter.GeneName.A", "inter.GeneName.B"))
		ret2.clustergroup.df <- Tool.SplitToGenDataFrame(ret.res[, "cluster.groups"],
			to.split.by = kClustersSplit, 
			res.colnames = c("Cluster.G1", "Cluster.G2"))
		ret.res <- cbind(ret.res[, "gene.pairs", drop = FALSE], 
			ret1.genepairs.df, ret.res[, "cluster.groups", drop = FALSE], 
			ret2.clustergroup.df, ret.res[, c("uq.cnt", "gp.logfc.calc")],
			stringsAsFactors = FALSE)
		return(ret.res)
	}

	## process: draw plots
	if (show.uq.cnt.merged == TRUE) {
		this.plot.mg <- inside.uq.single.plot(plot.data.uq.df, prioritize.cluster.groups, bar.colour)
		return(list(plot = this.plot.mg, table = inside.gen.ret.table(plot.data.uq.df)))
	} else {
		this.notm.plot.list <- list()
		this.notm.ret.table.list <- list()
		for (i.item in names(plot.data.uq.notm.list)) {
			this.notm.plot.list <- c(this.notm.plot.list, list(inside.uq.single.plot(plot.data.uq.notm.list[[i.item]], prioritize.cluster.groups, bar.colour)))
			this.notm.ret.table.list <- c(this.notm.ret.table.list, list(inside.gen.ret.table(plot.data.uq.notm.list[[i.item]])))
		}
		this.notm.plot <- plot_grid(plotlist = this.notm.plot.list, ncol = 1)

		return(list(plot = this.notm.plot, table = this.notm.ret.table.list))
	}
}


