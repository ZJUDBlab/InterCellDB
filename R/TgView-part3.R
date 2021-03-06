

#' Analyze Gene Pairs by Specificity
#'
#' @description
#' This function is used to evaluate specificity for gene pairs from one interaction. In other word, 
#' it explores every gene pair with all its occuring interactions. The one interaction is controlled by 
#' \code{\link{FetchInterOI}}.
#'
#' @inheritParams InsideObjectInterCell
#' @param to.cmp.cluster.groups It gives the interactions (cluster groups) to be compared. 
#'  \code{\link{ListClusterGroups}} helps to give the names of interactions.
#' @param extended.search Decide whether to add more gene pairs (those from compared cluster groups) into analysis.
#' @param uq.cnt.options Integer. It defines the allowed co-occurence count of one gene pair, and count goes larger impling 
#' lower specificity. The recommended value is \code{1:(length(to.cmp.cluster.groups) + 1)}.
#'
#' @return A \code{InterCell} object.
#'
#' @import dplyr
#' 
#' @export
#'
AnalyzeInterSpecificity <- function(
	object, 
	to.cmp.cluster.groups = character(),  # should be a subset of interact.pairs.acted$name.allpairs
	extended.search = FALSE, 
	uq.cnt.options = c(1:100)
) {
	interact.pairs.acted <- getFullViewResult(object)
	VEinfos <- getTgVEInfo(object)
	Itinfos <- getTgItInfo(object)
	kGenesSplit <- getGenePairSplit(object)
	kClustersSplit <- getClusterSplit(object)
	musthave.colnames <- object@misc$musthave.colnames

	# check those cluster groups to compare
	not.valid.cluster.groups <- setdiff(to.cmp.cluster.groups, interact.pairs.acted$name.allpairs)
	if (length(not.valid.cluster.groups) > 0) {
		warning("Given invalid cluster groups: ", paste0(not.valid.cluster.groups, collapse = ", ", ". "))
	}
	to.cmp.cluster.groups <- unique(intersect(to.cmp.cluster.groups, interact.pairs.acted$name.allpairs))
	# check 0 for to.cmp.cluster.groups is put below

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
	# this pair (handle .mirror situation)
	this.p.clusters <- getOrigClusterNameTgVEInfo(object)
	this.pair.name <- paste(this.p.clusters$cluster.name.A, this.p.clusters$cluster.name.B, sep = kClustersSplit)
	if ((this.pair.name %in% all.pairs.names) == FALSE) {
		stop("Given name of cluster pairs NOT exist!")
	}
	target.gene.pairs.df <- DataPrep.GenStdGenePairs.from.VEinfos(VEinfos, musthave.colnames)
	# add power and confidence here
	tmp.ref.colnames <- c(paste("inter", rep(c("GeneName"), each = 2), c("A", "B"), sep = "."), "inter.Power", "inter.Confidence")
	target.gene.pairs.df <- left_join(target.gene.pairs.df, Itinfos$bt.pairs[, tmp.ref.colnames], 
		by = c("inter.GeneName.A", "inter.GeneName.B"))
	
	this.clusters.separate <- strsplit(this.pair.name, split = kClustersSplit, fixed = TRUE)[[1]]
	this.pair.interacts <- inside.c.short.interacts(target.gene.pairs.df, kGenesSplit)
	this.gene.pairs <- data.frame(gp.name = this.pair.interacts, 
		gp.belongs = rep(this.pair.name, times = length(this.pair.interacts)), 
		gp.power = target.gene.pairs.df[, "inter.Power"], 
		gp.confidence = target.gene.pairs.df[, "inter.Confidence"],
		stringsAsFactors = FALSE)
	# further check if comparison of itself exist, which may cause error in downstream analysis
	if (this.pair.name %in% to.cmp.cluster.groups == TRUE) {
		warning("It cannot compare with itself. Automatically remove ", this.pair.name, " from to-comparing cluster groups.")
		to.cmp.cluster.groups <- setdiff(to.cmp.cluster.groups, this.pair.name)	
	}

	if (length(to.cmp.cluster.groups) == 0) {
		stop("No cluster groups are ready to compare. Please check parameter `to.cmp.cluster.groups`. ",
			 "By the way, the name of the interaction fetched by `FetchInterOI` will be removed from `to.cmp.cluster.groups` automatically if presented. ")
	}

	# set the uq.cnt range
	tmp.all.uq.clusters <- unique(as.character(unlist(strsplit(c(this.pair.name, to.cmp.cluster.groups), split = kClustersSplit, fixed = TRUE))))
	uq.cnt.options <- uq.cnt.options[which(uq.cnt.options %in% seq_along(tmp.all.uq.clusters))]
	if (length(uq.cnt.options) == 0) {
		stop("Please specify available uq.cnt.options to continue the analysis!\nThe allowed values are ", 
			paste0(seq_along(tmp.all.uq.clusters), collapse = ", "), ".")
	}

	# generate gene pairs to compare
	this.tg.pairs <- this.pair.interacts
	if (extended.search == TRUE) {
		tmp.added.tg <- as.character(unlist(lapply(to.cmp.cluster.groups, 
			all.pairs.interacts = all.pairs.interacts, tg.pairs = this.tg.pairs, tmp.sep = kGenesSplit,
			function(x, all.pairs.interacts, tg.pairs, tmp.sep) {
				inside.c.short.interacts(all.pairs.interacts[[x]], tmp.sep)
			})))
		this.tg.pairs <- unique(c(this.tg.pairs, tmp.added.tg))
	}  # only those pairs occurs in target cluster group will be further analyzed

	# fetch uq in cmp group
	other.gene.pairs.df.list <- lapply(to.cmp.cluster.groups, 
		all.pairs.interacts = all.pairs.interacts, tg.pairs = this.tg.pairs, tmp.sep = kGenesSplit, 
		function(x, all.pairs.interacts, tg.pairs, tmp.sep) {
			tmp.pairs.df <- all.pairs.interacts[[x]]
			tmp.pairs.interacts <- inside.c.short.interacts(tmp.pairs.df, tmp.sep)
			tmp.inds.tg <- which(tmp.pairs.interacts %in% tg.pairs)
			res.df <- data.frame(gp.name = tmp.pairs.interacts[tmp.inds.tg], 
				gp.belongs = rep(x, times = length(tmp.inds.tg)), 
				gp.power = tmp.pairs.df[tmp.inds.tg, "inter.Power"], 
				gp.confidence = tmp.pairs.df[tmp.inds.tg, "inter.Confidence"],
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
	tmp.inds.uq.valid <- this.uq.cnts %in% uq.cnt.options
	this.uq.cnts <- this.uq.cnts[tmp.inds.uq.valid]
	this.uq.indices <- this.uq.indices[tmp.inds.uq.valid]

	# collect by gene pairs each
	tmp.diff.res.mat <- vector(mode = "list", length = length(this.uq.indices))
	prog.bar.use.plot.collect <- progress::progress_bar$new(total = length(this.uq.indices))
	prog.bar.use.plot.collect$tick(0)
	tmp.diff.res.mat <- lapply(names(this.uq.indices), all.gp = this.all.gp.packed, 
	diff.cnts = this.uq.cnts, diff.indices = this.uq.indices, 
		function(x, all.gp, diff.cnts, diff.indices) {
			tmp.i <- which(names(diff.cnts) == x)  # must be the same as it is in diff.indices
			tmp.cnt <- diff.cnts[tmp.i]
			names(tmp.cnt) <- NULL  # remove the name
			tmp.indices <- diff.indices[[tmp.i]]
			tmp.properties <- all.gp[tmp.indices, c("gp.belongs", "gp.power", "gp.confidence")]
			tmp.properties <- tmp.properties[order(tmp.properties$gp.power, decreasing = TRUE), ]
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
			tmp.properties <- all.gp[tmp.indices, c("gp.belongs", "gp.power", "gp.confidence")]
			tmp.clusters.prop.df <- Tool.SplitToGenDataFrame(tmp.properties[, "gp.belongs"], to.split.by = kClustersSplit, 
				res.colnames = c("Cluster.X", "Cluster.Y"))
			tmp.res.df <- cbind(data.frame(gp.name = rep(tmp.genepair, times = nrow(tmp.properties)), 
					inter.GeneName.A = rep(tmp.gene.part[1], times = nrow(tmp.properties)),
					inter.GeneName.B = rep(tmp.gene.part[2], times = nrow(tmp.properties)),
					stringsAsFactors = FALSE), 
				tmp.properties[, "gp.belongs", drop = FALSE], tmp.clusters.prop.df, tmp.properties[, c("gp.power", "gp.confidence")], stringsAsFactors = FALSE)
			prog.bar.for.res.collect$tick()
			# result
			tmp.res.df
		})
	tmp.diff.df.res <- bind_rows(tmp.diff.df.list)

	object <- setTgSpGenes(object, list(res = tmp.diff.df.res, for.plot.use = tmp.diff.res.mat))
	return(object)
}





# GetResult.SummarySpecialGenes Select Gene Pairs Method: random
Inside.sel.gene.pairs.method.random.default <- function(
	this.spgenes, 
	sel.by.method.count,
	...
) {
	tmp.inds.sel <- sample(seq_along(this.spgenes), sel.by.method.count)
	sel.gene.pairs <- names(this.spgenes)[tmp.inds.sel]
	return(sel.gene.pairs)
}
Inside.sel.gene.pairs.method.random.IT <- Inside.sel.gene.pairs.method.random.default

# GetResult.SummarySpecialGenes Select Gene Pairs Method: By Power(decreasing or increasing)
Inside.sel.gene.pairs.method.power.default <- function(
	this.spgenes, 
	VEinfos, 
	sel.by.method.count, 
	sel.by.method.decreasing = TRUE,
	kClustersSplit,
	...
) {
	tmp.cluster.A <- VEinfos$cluster.name.A
	tmp.cluster.B <- VEinfos$cluster.name.B
	tmp.gp.name <- paste(tmp.cluster.A, tmp.cluster.B, sep = kClustersSplit)
	tmp.tg.calc <- unlist(lapply(seq_along(this.spgenes), 
		spgenes = this.spgenes, tg.gp.name = tmp.gp.name, 
		function(x, spgenes, tg.gp.name) {
			tmp.df <- spgenes[[x]]$uq.details
			tmp.df[which(tmp.df$gp.belongs == tg.gp.name), "gp.power"]
		}))
	# max -> min or min -> max
	tmp.inds.sel <- order(tmp.tg.calc, decreasing = sel.by.method.decreasing)
	length(tmp.inds.sel) <- sel.by.method.count
	tmp.inds.sel <- tmp.inds.sel[which(!is.na(tmp.inds.sel))]
	sel.gene.pairs <- names(this.spgenes)[tmp.inds.sel]
	return(sel.gene.pairs)
}
Inside.sel.gene.pairs.method.power.IT <- Inside.sel.gene.pairs.method.power.default

#' Get Result for Specificity Analysis
#'
#' @description
#' This function is to summary gene pair specificity by showing all gene pairs with their 
#' co-occurring interactions. Gene pairs are evaluated on both power and confidence like it in \code{\link{GetResultTgCrosstalk}}.
#' 
#' @inheritParams InsideObjectInterCell
#' @param sel.uq.cnt.options It defines the range of co-occurence count to be shown in result. If no option is given,
#'  it will use all available co-occurence count.
#' @param sel.gene.pairs Directly specify the desired gene pairs. It should be given in standard table that is generated 
#'  by \code{\link{FormatCustomGenePairs}}. The order of gene pairs will be implicitly used for result graph and table.
#'  In addition, it's strictly aligned to clusters, which is illustrated in \code{\link{AnalyzeInterInFullView}}.
#' @param sel.gene.pairs.method Options are: "random", "by-power" or not to use (by set NULL). It works only when no specific gene pairs are given in 
#'  parameter \code{sel.gene.pairs}. Method "random" will randomly pick up some gene pairs. Method "by-power" will calculate the power of 
#'  the 2 genes in every gene pairs, and order them by the calculated values. 
#'  The order of gene pairs generated by method will be implicitly used for result graph and table.
#' @param sel.by.method.count It defines the maximum number of gene pairs to be fetched by any method.
#' @param sel.by.method.decreasing It works for method "by-power". If TRUE, result will be ordered in decreasing way, otherwise the increasing direction.
#' @param prioritize.cluster.groups It defines the most concerning cluster groups, and those cluster groups given in this parameter, will be 
#'  finally plotted from the most left-side to right. 
#' @param plot.name.X.to.Y If set FALSE, switch the position of 2 involving clusters in the original names of interaction. If set TRUE, keep still.
#' @param plot.uq.cnt.merged If set TRUE, then gene pairs of different count of co-occurence will be merged to show in the result, or gene pairs 
#'  will be grouped by their count of co-occurence and be plotted accordingly.
#' @param bound.to.use It specifies the user specified bound for evaluation params. The values out of bound will be coerced to 
#'  either lower bound or upper bound. Default no bound is set, i.e. [-Inf, +Inf]
#' @param func.to.use The function used to further transform the values, e.g. \code{log1p}. 
#' @param grid.plot.ncol It gives the number of columns for plotting grid when \code{plot.uq.cnt.merged = FALSE}.
#' @param barplot.or.dotplot If TRUE, use 'barplot', or, use 'dotplot'.
#' @param plot.font.size.base It defines the base font size for all texts. 
#' @param axis.text.x.pattern It defines the axis text style in x-axis. 
#' @param bar.facet.scales It controls the scales that facet uses, and gets 4 options as defined by \pkg{ggplot2}: "fixed", "free", "free_x", "free_y".
#' @param bar.facet.space It controls the space allocating strategy that facet uses, and gets 4 options as defined by \pkg{ggplot2}: "fixed", "free", "free_x", "free_y".
#' @param bar.facet.text.x It defines the style of facet text on the top horizontal position. 
#' @param bar.facet.background It defines the background style of labels among facets.
#' @param bar.colour It gives colors that plotting bars get to use. If no specific colour is given, then the 
#'  built-in 20 kinds of colours will be used.
#' @param bar.width It defines the bar width.
#' @param dot.colour.seq It specifies the colour sequence for dots.
#' @param dot.colour.value.seq It is along with the param \code{dot.colour.seq}, and changes the colour expansion.
#' @param dot.size.range It specifies the size range of the dots
#' @param dot.y.order.in.alphabet If set TRUE, force to order gene pairs in alphabet,  
#'  which ignores the implicit order by \code{sel.gene.pairs} or \code{sel.gene.pairs.method}.
#' @param ... Other parameter that can be passed to select by method functions.
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import cowplot
#' 
#' @export
#'
GetResultTgSpecificity <- function(
	object,
	sel.uq.cnt.options = integer(), 
	sel.gene.pairs = NULL, 
	sel.gene.pairs.method = NULL,
	sel.by.method.count = 10, 
	sel.by.method.decreasing = TRUE, 
	prioritize.cluster.groups = character(),
	plot.name.X.to.Y = TRUE,
	plot.uq.cnt.merged = TRUE,
	bound.to.use = list("Power" = c(-Inf, +Inf), "Confidence" = c(-Inf, +Inf)), 
	func.to.use = list("Power" = function(x) {x}, "Confidence" = function(x) {1/(1+x)}),
	grid.plot.ncol = 1, 
	barplot.or.dotplot = FALSE,
	plot.font.size.base = 12, 
	axis.text.x.pattern = element_text(angle = 90, vjust = 0.5, hjust = 1),
	bar.facet.scales = "free_x", 
	bar.facet.space  = "free_x", 
	bar.facet.text.x = element_text(size = 8, colour = "black"), 
	bar.facet.background = element_rect(fill = "lightgrey", colour = "white"), 
	bar.colour = character(),
	bar.width = 0.8, 
	dot.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
	dot.colour.value.seq = c(0.0, 0.5, 1.0),
	dot.size.range = c(2, 8), 
	dot.y.order.in.alphabet = TRUE,
	...
) {
	kGenesSplit <- getGenePairSplit(object)
	kClustersSplit <- getClusterSplit(object)
	VEinfos <- getTgVEInfo(object)
	onepair.spgenes <- getTgSpGenes(object)
	this.spgenes <- onepair.spgenes$for.plot.use
	#
	show.topn.inside.onepair <- +Inf

	## pre-check
	this.property.valid.uq.cnt <- unique(as.integer(unlist(lapply(this.spgenes, function(x) { x$uq.cnt }))))
	# check valid uq.cnt
	if (length(sel.uq.cnt.options) == 0) {  # if length uq.cnt = 0 then set it directly
		sel.uq.cnt.options <- this.property.valid.uq.cnt
	} else {
		# check valid
		tmp.inds.valid.uq <- which(sel.uq.cnt.options %in% this.property.valid.uq.cnt)
		if (length(tmp.inds.valid.uq) != length(sel.uq.cnt.options)) {
			warning("Select not-existed uq.cnt: ", paste0(sel.uq.cnt.options[setdiff(seq_along(sel.uq.cnt.options), tmp.inds.valid.uq)], collapse = ", "), 
				", which will be ignored!")
		}
		sel.uq.cnt.options <- sel.uq.cnt.options[tmp.inds.valid.uq]
	}

	# check and properly set cluster group given
	this.property.valid.cluster.group <- unique(as.character(unlist(lapply(this.spgenes, function(x) { x$uq.details$gp.belongs }))))
	# if not set prioritize.cluster.groups, use the original order of cluster groups
	if (length(prioritize.cluster.groups) == 0) {
		prioritize.cluster.groups <- this.property.valid.cluster.group
	}
	tmp.inds.valid.cluster.group <- which(prioritize.cluster.groups %in% this.property.valid.cluster.group)
	if (length(tmp.inds.valid.cluster.group) != length(prioritize.cluster.groups)) {
		warning("Given cluster group order has some items not existed: ",
			paste0(prioritize.cluster.groups[setdiff(seq_along(prioritize.cluster.groups), tmp.inds.valid.cluster.group)], collapse = " "),
			", which will be automatically removed!")
	}
	prioritize.cluster.groups <- prioritize.cluster.groups[tmp.inds.valid.cluster.group]


	# check dot plot parameters
	allowed.data.modf.items <- c("Power", "Confidence")
	given.bound.to.use <- CheckParamStd(names(bound.to.use), allowed.data.modf.items, "bound", stop.on.zero = FALSE)
	tmp.not.in.bound <- setdiff(c("Power", "Confidence"), given.bound.to.use)
	if (length(tmp.not.in.bound) != 0) {
		for (tmp.i in tmp.not.in.bound) {
			bound.to.use <- c(bound.to.use, list(c(-Inf, Inf)))
			names(bound.to.use)[length(bound.to.use)] <- tmp.i
		}
		warning("Auto add bounds on: ", paste0(tmp.not.in.bound, collapse = ", "), ".")
	}
	#
	if (length(dot.colour.seq) != length(dot.colour.value.seq)) {
		stop("Colors and their gradient values should be of same length! Check parameter `dot.colour.seq` and `dot.colour.value.seq`.")
	}

	# check func to modify data
	tmp.valid.funcnames <- CheckParamStd(names(func.to.use), allowed.data.modf.items, "Transform functions", stop.on.zero = FALSE)
  func.to.use <- func.to.use[tmp.valid.funcnames]
  tmp.not.in.func <- setdiff(allowed.data.modf.items, tmp.valid.funcnames)
  if (length(tmp.not.in.func) != 0) {
    for (tmp.j in tmp.not.in.func) {
      func.to.use <- c(func.to.use, list(function(x) {x}))
      names(func.to.use)[length(func.to.use)] <- tmp.j
    }
    warning("Auto add transform function on: ", paste0(tmp.not.in.func, collapse = ", "), ".")
  }

	# fetch selected gene pairs
	inside.fetch.sel.genepairs <- function(
		input.spgenes,
		VEinfos,
		sel.gene.pairs,
		sel.gene.pairs.method,
		sel.by.method.count,
		sel.by.method.decreasing,
		...
	) {
		if (length(sel.gene.pairs) == 0) {  
			if (is.null(sel.gene.pairs.method) || length(sel.gene.pairs.method) == 0) {
				sel.gene.pairs <- names(input.spgenes)  # use all 
			} else {  # then use method
				# check maximum pairs limit
				if (sel.by.method.count > length(input.spgenes)) {
					warning("Maximum gene pairs are:", length(input.spgenes), ", and given limit is automatically shrinked to that value.")
					sel.by.method.count <- length(input.spgenes)
				}
				# select by methods
				sel.gene.pairs <- switch(sel.gene.pairs.method, 
					"random" = Inside.sel.gene.pairs.method.random.IT(input.spgenes, sel.by.method.count, ...), 
					"by-power" = Inside.sel.gene.pairs.method.power.IT(input.spgenes, VEinfos, sel.by.method.count, sel.by.method.decreasing, kClustersSplit, ...), 
					stop("Undefined method for selecting gene pairs. Please re-check the given param: sel.gene.pairs.method!!")
				)
			}
		} else {  # check validity of those select gene pairs
			## check format of given `sel.gene.pairs`
			# check data.frame
			if (!is.data.frame(sel.gene.pairs)) {
				stop("User-selected gene pairs should be given in table (data.frame format)!")
			}
			# check if it is standardized data.frame
			std.colnames.1t4 <- c("inter.GeneID.A", "inter.GeneID.B", "inter.GeneName.A", "inter.GeneName.B")
			if (!identical(colnames(sel.gene.pairs)[1:4], std.colnames.1t4)) {
				stop("Non-standardized table of gene pairs are given, please use function `FormatCustomGenePairs` first to get standardized one!")
			}

			# transform data.frame to character
			sel.gene.pairs <- paste(sel.gene.pairs[, "inter.GeneName.A"], sel.gene.pairs[, "inter.GeneName.B"], sep = kGenesSplit)

			tmp.inds.vd.gp <- which(sel.gene.pairs %in% names(input.spgenes))
			if (length(tmp.inds.vd.gp) != length(sel.gene.pairs)) {
				warning("Selected gene pairs not exist: ", 
					paste0(sel.gene.pairs[setdiff(seq_along(sel.gene.pairs), tmp.inds.vd.gp)], collapse = ", "), 
					", which will be removed from following analysis!")
			}
			sel.gene.pairs <- sel.gene.pairs[tmp.inds.vd.gp]
		}
		# unique on gene pairs
		sel.gene.pairs <- unique(sel.gene.pairs)
		return(sel.gene.pairs)
	}
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
			list(uq.cnt = x$uq.cnt, uq.details = x$uq.details[seq_len(tmp.selrows), c("gp.belongs", "gp.power", "gp.confidence")])
			})
		# remove those with NO valid uq.cnt 
		tmp.check.0row <- unlist(lapply(tmp.spgenes, function(x) { nrow(x$uq.details) }))
		tmp.spgenes <- tmp.spgenes[which(tmp.check.0row > 0)]
		return(tmp.spgenes)
	}

	inside.trans.coords.uq <- function(
		std.spgenes,
		show.genepairs.order,  # ordered as it is in sel.gene.pairs
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
					uq.y.axis = std.spgenes[[x]]$uq.details$gp.power[c(tmp.inds.prior, tmp.inds.inferior)],
					uq.z.axis = std.spgenes[[x]]$uq.details$gp.confidence[c(tmp.inds.prior, tmp.inds.inferior)],
					stringsAsFactors = FALSE)
				})
		tmp.df.res <- bind_rows(tmp.df.list)
		if (nrow(tmp.df.res) == 0) {
			stop("No available gene pairs in current settings. Consider enlarge the range of `sel.uq.cnt.options`, like set 1:100 or larger.")
		}
		return(tmp.df.res)
	}

	## process: reconstruct data format of every uq.cnt to fit plotting
	if (plot.uq.cnt.merged == TRUE) {
		plot.data.uq <- list()
		for (iuq in sel.uq.cnt.options) {
			tmp.spgenes <- inside.collect.uq.cnt.each(this.spgenes, iuq, show.topn.inside.onepair)
			plot.data.uq <- c(plot.data.uq, tmp.spgenes)
		}
		tmp.sel.gpairs <- inside.fetch.sel.genepairs(plot.data.uq, VEinfos, sel.gene.pairs, sel.gene.pairs.method, sel.by.method.count, sel.by.method.decreasing, ...)
		plot.data.uq <- plot.data.uq[which(names(plot.data.uq) %in% tmp.sel.gpairs)]
		plot.data.uq.df <- inside.trans.coords.uq(plot.data.uq, tmp.sel.gpairs, prioritize.cluster.groups)
	} else {
		plot.data.uq.notm.list <- list()
		tmp.uq.cnt.valid.list <- integer()
		for (iuq in sel.uq.cnt.options) {
			tmp.spgenes <- inside.collect.uq.cnt.each(this.spgenes, iuq, show.topn.inside.onepair)
			tmp.sel.gpairs <- inside.fetch.sel.genepairs(tmp.spgenes, VEinfos, sel.gene.pairs, sel.gene.pairs.method, sel.by.method.count, sel.by.method.decreasing, ...)
			tmp.spgenes <- tmp.spgenes[which(names(tmp.spgenes) %in% tmp.sel.gpairs)]
			if (length(tmp.spgenes) == 0) {
				next  # not added to the result list
			}
			tmp.spgenes.trans <- inside.trans.coords.uq(tmp.spgenes, tmp.sel.gpairs, prioritize.cluster.groups)
			plot.data.uq.notm.list <- c(plot.data.uq.notm.list, list(tmp.spgenes.trans))
			tmp.uq.cnt.valid.list <- c(tmp.uq.cnt.valid.list, iuq)
		}
		names(plot.data.uq.notm.list) <- paste("uq.cnt=", as.character(tmp.uq.cnt.valid.list), sep = ".")
	}

	inside.uq.bar.plot <- function(
		plot.data,
		prioritize.cluster.groups,
		plot.name.X.to.Y,
		bar.colour
	) {
		# re-align cluster orders, conv or rev
		if (plot.name.X.to.Y == FALSE) {
			uq.name.split.df <- Tool.SplitToGenDataFrame(plot.data[, "uq.name"],
				to.split.by = kGenesSplit, res.colnames = c("gene.A", "gene.B"))
			uq.label.split.df <- Tool.SplitToGenDataFrame(plot.data[, "uq.label"],
				to.split.by = kClustersSplit, res.colnames = c("cluster.A", "cluster.B"))
			plot.data$uq.name <- paste(uq.name.split.df$gene.B, uq.name.split.df$gene.A, sep = kGenesSplit)
			plot.data$uq.label <- paste(uq.label.split.df$cluster.B, uq.label.split.df$cluster.A, sep = kClustersSplit)
			#
			# prioritized cluster group changes as well
			tmp.rev.pcg <- Tool.SplitToGenDataFrame(prioritize.cluster.groups,
				to.split.by = kClustersSplit, res.colnames = c("cluster.A", "cluster.B"))
			prioritize.cluster.groups <- paste(tmp.rev.pcg$cluster.B, tmp.rev.pcg$cluster.A, sep = kClustersSplit)
		}  # not change if set as conv

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
		plot.data$uq.name.align <- factor(plot.data$uq.name, levels = unique(plot.data$uq.name))
		# the plot
		this.plot <- ggplot(plot.data, aes(x = uq.label, y = uq.y.axis))
		this.plot <- this.plot + 
			geom_col(aes(fill = uq.label), width = bar.width) + 
			facet_grid(cols = vars(uq.name.align), scales = bar.facet.scales, space = bar.facet.space) + 
			scale_x_discrete(breaks = plot.data$uq.label, 
				limits = tiny.cg.prior,   # function to change the limit
				labels = plot.data$uq.label) + 
			scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # keep the default 5% in y top
			scale_fill_manual(name = "Cluster Group", 
				values = colour.sel.it, 
				breaks = colour.sel.cor.breaks) + 
			labs(x = "Cluster Groups", y = "Power")
		this.plot <- this.plot + 
			theme_half_open(font_size = plot.font.size.base) + 
			theme(axis.text.x = axis.text.x.pattern,
				strip.text.x = bar.facet.text.x,
				strip.background = bar.facet.background) + 
			theme(legend.position = "none")  # remove the legends
		return(this.plot)
	}

	inside.uq.dot.plot <- function(
		plot.data,
		prioritize.cluster.groups,
		plot.name.X.to.Y,
		bound.to.use,
		func.to.use,
		dot.colour.seq,
		dot.colour.value.seq,
		dot.size.range
	) {
		# re-align cluster orders, conv or rev
		if (plot.name.X.to.Y == FALSE) {
			uq.name.split.df <- Tool.SplitToGenDataFrame(plot.data[, "uq.name"],
				to.split.by = kGenesSplit, res.colnames = c("gene.A", "gene.B"))
			uq.label.split.df <- Tool.SplitToGenDataFrame(plot.data[, "uq.label"],
				to.split.by = kClustersSplit, res.colnames = c("cluster.A", "cluster.B"))
			plot.data$uq.name <- paste(uq.name.split.df$gene.B, uq.name.split.df$gene.A, sep = kGenesSplit)
			plot.data$uq.label <- paste(uq.label.split.df$cluster.B, uq.label.split.df$cluster.A, sep = kClustersSplit)
			#
			# prioritized cluster group changes as well
			tmp.rev.pcg <- Tool.SplitToGenDataFrame(prioritize.cluster.groups,
				to.split.by = kClustersSplit, res.colnames = c("cluster.A", "cluster.B"))
			prioritize.cluster.groups <- paste(tmp.rev.pcg$cluster.B, tmp.rev.pcg$cluster.A, sep = kClustersSplit)
		}  # not change if set as conv


		# move the data range
		#plot.data$uq.y.axis  # Power
		#plot.data$uq.z.axis  # confidence
		# get and recalc Power range
		tmp.power.range <- bound.to.use[["Power"]]
		tmp.power <- plot.data$uq.y.axis
		# modify
		tmp.power <- func.to.use[["Power"]](tmp.power)
		tmp.power[which(tmp.power > tmp.power.range[2])] <- tmp.power.range[2]
		tmp.power[which(tmp.power < tmp.power.range[1])] <- tmp.power.range[1]
		plot.data$plot.dot.colour <- tmp.power

		# get and Confidence range
		tmp.confidence.range <- bound.to.use[["Confidence"]]
		tmp.confidence <- plot.data$uq.z.axis
		# modify
		tmp.confidence <- func.to.use[["Confidence"]](tmp.confidence)
		tmp.confidence[which(tmp.confidence > tmp.confidence.range[2])] <- tmp.confidence.range[2]
		tmp.confidence[which(tmp.confidence < tmp.confidence.range[1])] <- tmp.confidence.range[1]
		plot.data$plot.dot.size <- tmp.confidence

		# before plot, get gene pairs ordered correctly
		if (dot.y.order.in.alphabet == TRUE) {
			plot.data$uq.name <- factor(plot.data$uq.name, levels = unique(plot.data$uq.name[order(plot.data$uq.name)]))
		} else {
			plot.data$uq.name <- factor(plot.data$uq.name, levels = unique(plot.data$uq.name))
		}
		
		# add additional re-order step, as dot plot is not the same as bar plot.
		tmp.uq.label <- unique(plot.data$uq.label)
		tmp.inds.prior <- match(prioritize.cluster.groups, tmp.uq.label)
		tmp.uq.label.inferior <- tmp.uq.label[setdiff(seq_along(tmp.uq.label), tmp.inds.prior)]
		tmp.align.uq.label <- c(tmp.uq.label[tmp.inds.prior], tmp.uq.label.inferior[order(tmp.uq.label.inferior)])
		plot.data$uq.label <- factor(plot.data$uq.label, levels = tmp.align.uq.label)
		# the plot
		this.plot <- ggplot(plot.data, aes(x = uq.label, y = uq.name))
		this.plot <- this.plot + 
			geom_point(aes(colour = plot.dot.colour, size = plot.dot.size)) + 
			scale_colour_gradientn(name = "Power", colours = dot.colour.seq, values = dot.colour.value.seq) + 
			scale_size(name = "Confidence", range = dot.size.range) + 
			labs(x = "Cluster Groups", y = "Gene Pairs")
		this.plot <- this.plot + 
			theme_half_open(font_size = plot.font.size.base) + 
			theme_bw() + 
			theme(axis.text.x = axis.text.x.pattern) + 
			theme(legend.position = "right")  # remove the legends
		return(this.plot)		
	}

	inside.gen.ret.table <- function(
		ret.data
	) {
		ret.res <- ret.data[, c("uq.name", "uq.label", "uq.cnt", "uq.y.axis", "uq.z.axis")]
		colnames(ret.res) <- c("gene.pairs", "cluster.groups", "uq.cnt", "gp.power", "gp.confidence")
		ret1.genepairs.df <- Tool.SplitToGenDataFrame(ret.res[, "gene.pairs"], 
			to.split.by = kGenesSplit, 
			res.colnames = c("inter.GeneName.A", "inter.GeneName.B"))
		ret2.clustergroup.df <- Tool.SplitToGenDataFrame(ret.res[, "cluster.groups"],
			to.split.by = kClustersSplit, 
			res.colnames = c("Cluster.X", "Cluster.Y"))
		ret.res <- cbind(ret.res[, "gene.pairs", drop = FALSE], 
			ret1.genepairs.df, ret.res[, "cluster.groups", drop = FALSE], 
			ret2.clustergroup.df, ret.res[, c("uq.cnt", "gp.power", "gp.confidence")],
			stringsAsFactors = FALSE)
		return(ret.res)
	}

	## process: draw plots
	if (barplot.or.dotplot == TRUE) {
		if (plot.uq.cnt.merged == TRUE) {
			this.plot.mg <- inside.uq.bar.plot(plot.data.uq.df, prioritize.cluster.groups, plot.name.X.to.Y, bar.colour)
			return(list(plot = this.plot.mg, table = inside.gen.ret.table(plot.data.uq.df)))
		} else {
			this.notm.plot.list <- list()
			this.notm.ret.table.list <- list()
			for (i.item in names(plot.data.uq.notm.list)) {
				this.notm.plot.list <- c(this.notm.plot.list, list(inside.uq.bar.plot(plot.data.uq.notm.list[[i.item]], prioritize.cluster.groups, plot.name.X.to.Y, bar.colour)))
				this.notm.ret.table.list <- c(this.notm.ret.table.list, list(inside.gen.ret.table(plot.data.uq.notm.list[[i.item]])))
			}
			this.notm.plot <- plot_grid(plotlist = this.notm.plot.list, ncol = grid.plot.ncol)

			return(list(plot = this.notm.plot, table = this.notm.ret.table.list))
		}
	} else {
		if (plot.uq.cnt.merged == TRUE) {
			this.plot.mg <- inside.uq.dot.plot(plot.data.uq.df, prioritize.cluster.groups, plot.name.X.to.Y, bound.to.use, func.to.use, dot.colour.seq, dot.colour.value.seq, dot.size.range)
			return(list(plot = this.plot.mg, table = inside.gen.ret.table(plot.data.uq.df)))
		} else {
			this.notm.plot.list <- list()
			this.notm.ret.table.list <- list()
			for (i.item in names(plot.data.uq.notm.list)) {
				this.notm.plot.list <- c(this.notm.plot.list, list(inside.uq.dot.plot(plot.data.uq.notm.list[[i.item]], prioritize.cluster.groups, plot.name.X.to.Y, bound.to.use, func.to.use, dot.colour.seq, dot.colour.value.seq, dot.size.range)))
				this.notm.ret.table.list <- c(this.notm.ret.table.list, list(inside.gen.ret.table(plot.data.uq.notm.list[[i.item]])))
			}
			this.notm.plot <- plot_grid(plotlist = this.notm.plot.list, ncol = grid.plot.ncol)

			return(list(plot = this.notm.plot, table = this.notm.ret.table.list))
		}
	}
}


