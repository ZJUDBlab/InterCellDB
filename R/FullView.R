

#' Analyze interaction network in full view
#'
#' @description
#' This function analyzes count and power of interaction pairs among all given clusters.
#'
#' @param object [TODO] InterCell Object
#' @param sel.clusters.X [TODO]
#' @param sel.clusters.Y [TODO]
#' @param sel.exprs.change [TODO]
#' @param sel.some.genes.X [TODO]
#' @param sel.some.genes.Y [TODO]
#' @param sel.genes.option [TODO]
#' @param sel.gene.pairs [TODO] data.frame
#' @param force.process [TODO]
#' @param verbose [TODO]
#'
#' @details
#'
#'
#' @return A list.
#'
#' @examples
#'
#'
#' @importFrom dplyr bind_rows left_join
#' @import progress
#'
#' @export
#'
AnalyzeInterInFullView <- function(
	object,
	sel.clusters.X = NULL,
	sel.clusters.Y = NULL,
	sel.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"),
	sel.some.genes.X = NULL,
	sel.some.genes.Y = NULL,
	sel.genes.option = "intersect", 
	sel.gene.pairs = NULL,  # [TODO] modify and set an standard AddUserPairs function
	force.process = FALSE,  # control if no gene selection applied
	verbose = TRUE
) {
	# check sel.clusters.*
	inside.check.sel.clusters <- function(allowed.clusters, sel.clusters, suffix = c("X", "Y")) {
		if (is.null(sel.clusters)) {
			sel.clusters <- allowed.clusters
			sel.clusters <- sel.clusters[order(sel.clusters)]
		} else {
			not.valid.clusters <- setdiff(sel.clusters, allowed.clusters)
			if (length(not.valid.clusters) > 0) {
				warning("Given Undefined clusters in `sel.clusters.", suffix, "`: ", paste0(not.valid.clusters, collapse = ", "), ". ")
			}
			sel.clusters <- intersect(sel.clusters, allowed.clusters)
		}
		return(unique(sel.clusters))
	}
	# check
	all.clusters <- unique(object@fgenes$cluster)
	sel.clusters.X <- inside.check.sel.clusters(all.clusters, sel.clusters.X, suffix = "X")
	sel.clusters.Y <- inside.check.sel.clusters(all.clusters, sel.clusters.Y, suffix = "Y")

	# check selection on gene pairs or some genes
	if (!is.null(sel.gene.pairs)) {
		# check data.frame
		if (!is.data.frame(sel.gene.pairs)) {
			stop("User-selected gene pairs should be given in table (data.frame format)!")
		}
		# check if it is standardized data.frame
		std.colnames.1t4 <- c("inter.GeneID.A", "inter.GeneID.B", "inter.GeneName.A", "inter.GeneName.B")
		if (!identical(colnames(sel.gene.pairs)[1:4], std.colnames.1t4)) {
			stop("Non-standardized table of gene pairs are given, please use `[TODO] <function>` firstly to get standardized one!")
		}
	}
	if (!is.null(sel.some.genes.X)) {
		if (!is.character(sel.some.genes.X)) {
			sel.some.genes.X <- as.character(sel.some.genes.X)
		}
	}
	if (!is.null(sel.some.genes.Y)) {
		if (!is.character(sel.some.genes.Y)) {
			sel.some.genes.Y <- as.character(sel.some.genes.Y)
		}
	}
	# gene pairs if given, genes given will be overwritten
	if (is.null(sel.gene.pairs)) {
		if (is.null(sel.some.genes.X) && is.null(sel.some.genes.Y)) {
			if (force.process == FALSE) {
				stop("No restriction applied on neither genes or gene pairs.",
					" Program will get long time to be finished.",
					" If insisting, set `force.process` to be TRUE.")
			}
		} else {
			if (is.null(sel.some.genes.X) || is.null(sel.some.genes.Y)) {
				print("Processing with partial selection on genes.")
			}
		}
	} else {
		print(paste("Processing with given", nrow(sel.gene.pairs), "gene pairs."))
		if (!is.null(sel.some.genes.X) || !is.null(sel.some.genes.Y)) {
			warning("Given selection on genes is overrided by given selection on gene pairs!")
		}
	}
	# check sel.genes.option
	std.sel.genes.option <- c("intersect", "union")
	sel.genes.func <- base::intersect  # in default using 'intersect'
	sel.genes.option <- sel.genes.option[1]  # only length 1 allowed
	valid.option <- intersect(sel.genes.option, std.sel.genes.option)
	sel.genes.func <- switch(valid.option, 
		"intersect" = base::intersect,
		"union" = base::union,
		warning("Given invalid user option on selecting genes: ", paste0(sel.genes.option, collapse = ", "), ". ", 
			"Using default 'intersect' option!")
	)
	

	# check exprs change
	not.valid.exprs.c <- setdiff(sel.exprs.change, kexprs.change)
	if (length(not.valid.exprs.c) > 0) {
		warning("Given invalid options for expression change: ", paste0(not.valid.exprs.c, collapse = ", ", ". "))
	}
	sel.exprs.change <- intersect(sel.exprs.change, kexprs.change)
	if (length(sel.exprs.change) == 0) {
		stop("No expression change selection is given in parameter `sel.exprs.change`. Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'.")
	}


	## analyze interaction network
	interact.pairs.all <- list(
		list.clusters = list(x.axis = sel.clusters.X, y.axis = sel.clusters.Y),
		data.allpairs = list(),
		anno.allpairs = list(location.A = list(), location.B = list(), type.A = list(), type.B = list()),
		name.allpairs = character(),
		cnt.allpairs = integer(),
		strength.allpairs = single()
	)

	if (verbose == TRUE) {
		prog.bar.fv <- progress::progress_bar$new(total = length(sel.clusters.X) * length(sel.clusters.Y))
		prog.bar.fv$tick(0)
	}
	for (ix in sel.clusters.X) {
		for (jy in sel.clusters.Y) {
			interact.name <- paste0(ix, object@tool.vars$cluster.split, jy)
			fgenes.t.X <- object@fgenes[which(object@fgenes$cluster == ix), ]
			fgenes.t.Y <- object@fgenes[which(object@fgenes$cluster == jy), ]
			genes.X <- fgenes.t.X$gene
			genes.Y <- fgenes.t.Y$gene

			# generate all possible gene pairs
			ref.GeneA <- object@database@pairs.db$inter.GeneName.A
			ref.GeneB <- object@database@pairs.db$inter.GeneName.B
			# conv
			inds.gpairs.conv <- intersect(which(ref.GeneA %in% genes.X), which(ref.GeneB %in% genes.Y))
			gpairs.conv <- object@database@pairs.db[inds.gpairs.conv, ]
			# rev
			inds.gpairs.rev <- intersect(which(ref.GeneB %in% genes.X), which(ref.GeneA %in% genes.Y))
			gpairs.rev <- object@database@pairs.db[inds.gpairs.rev, 
				c(ReverseOddEvenCols(std.index.colname.end.dual), (std.index.colname.end.dual+1):ncol(object@database@pairs.db))]
			colnames(gpairs.rev) <- colnames(object@database@pairs.db)
			# result
			gpairs.result <- rbind(gpairs.conv, gpairs.rev)

			# prioritized usage of selection to process
			if (!is.null(sel.gene.pairs)) {  # use gene pairs
#[TODO]
			} else {  # use selected genes
				inds.sel.X <- inds.sel.Y <- seq_len(nrow(gpairs.result))
				if (!is.null(sel.some.genes.X)) {
					inds.sel.X <- which(gpairs.result$inter.GeneName.A %in% sel.some.genes.X)
				}
				if (!is.null(sel.some.genes.Y)) {
					inds.sel.Y <- which(gpairs.result$inter.GeneName.B %in% sel.some.genes.Y)
				}
				# use select gene merge option and get result
				gpairs.result <- gpairs.result[sel.genes.func(inds.sel.X, inds.sel.Y), ]
			}

			# add logfc, pval data from fgenes
			tmp.sel.cols <- setdiff(object@misc$musthave.colnames, "cluster")
			gpairs.result <- left_join(gpairs.result, fgenes.t.X[, tmp.sel.cols], by = c("inter.GeneName.A" = "gene"))
			tmp.change.cols <- ncol(gpairs.result) - (length(tmp.sel.cols) - 2):0
			colnames(gpairs.result)[tmp.change.cols] <- paste("inter", colnames(gpairs.result)[tmp.change.cols], "A", sep = ".")
			gpairs.result <- left_join(gpairs.result, fgenes.t.Y[, tmp.sel.cols], by = c("inter.GeneName.B" = "gene"))
			tmp.change.cols <- ncol(gpairs.result) - (length(tmp.sel.cols) - 2):0
			colnames(gpairs.result)[tmp.change.cols] <- paste("inter", colnames(gpairs.result)[tmp.change.cols], "B", sep = ".")
			# reorder the result table
			tmp.paired.order.cols <- paste("inter", rep(c("GeneID", "GeneName", "LogFC", "PVal"), each = 2), c("A", "B"), sep = ".")
			gpairs.result <- gpairs.result[, c(tmp.paired.order.cols, setdiff(colnames(gpairs.result), tmp.paired.order.cols))]
			
			# check exprs change
			inds.subg.logfc <- integer()
			inds.A.up <- which(gpairs.result$inter.LogFC.A > 0)
			inds.A.dn <- which(gpairs.result$inter.LogFC.A <= 0)
			inds.B.up <- which(gpairs.result$inter.LogFC.B > 0)
			inds.B.dn <- which(gpairs.result$inter.LogFC.B <= 0)
			if ("Xup.Yup" %in% sel.exprs.change) {
				inds.subg.logfc <- c(inds.subg.logfc, intersect(inds.A.up, inds.B.up))
			}
			if ("Xup.Ydn" %in% sel.exprs.change) {
				inds.subg.logfc <- c(inds.subg.logfc, intersect(inds.A.up, inds.B.dn))
			}
			if ("Xdn.Yup" %in% sel.exprs.change) {
				inds.subg.logfc <- c(inds.subg.logfc, intersect(inds.A.dn, inds.B.up))
			}
			if ("Xdn.Ydn" %in% sel.exprs.change) {
				inds.subg.logfc <- c(inds.subg.logfc, intersect(inds.A.dn, inds.B.dn))
			}
			gpairs.result <- gpairs.result[inds.subg.logfc, ]

			# further unique result
			gpairs.result <- DoPartUnique(gpairs.result, c(1,2))

			# get result for other attributes
			# location
			this.A.locations <- object@database@anno.location.db[which(object@database@anno.location.db$GeneID %in% gpairs.result$inter.GeneID.A), ]
			this.B.locations <- object@database@anno.location.db[which(object@database@anno.location.db$GeneID %in% gpairs.result$inter.GeneID.B), ]
			# type
			this.A.types <- object@database@anno.type.db[which(object@database@anno.type.db$GeneID %in% gpairs.result$inter.GeneID.A), ]
			this.B.types <- object@database@anno.type.db[which(object@database@anno.type.db$GeneID %in% gpairs.result$inter.GeneID.B), ]
			
			# get result saved
			interact.pairs.all$data.allpairs[[interact.name]] <- gpairs.result
			interact.pairs.all$anno.allpairs$location.A[[interact.name]] <- this.A.locations
			interact.pairs.all$anno.allpairs$location.B[[interact.name]] <- this.B.locations
			interact.pairs.all$anno.allpairs$type.A[[interact.name]] <- this.A.types
			interact.pairs.all$anno.allpairs$type.B[[interact.name]] <- this.B.types
			interact.pairs.all$name.allpairs <- c(interact.pairs.all$name.allpairs, interact.name)
			interact.pairs.all$cnt.allpairs <- c(interact.pairs.all$cnt.allpairs, nrow(gpairs.result))
			interact.pairs.all$strength.allpairs <- c(interact.pairs.all$strength.allpairs,
				object@formulae$FULLVIEW(gpairs.result, c("inter.LogFC.A", "inter.LogFC.B")))
			if (verbose == TRUE) {
				prog.bar.fv$tick()
			}
		}
	}

	# result
	object <- setFullViewResult(object, interact.pairs.all)
	return(object)
}





#' Get result of analysis in full view
#' 
#' @description
#' This function focuses on some subset of interaction pairs, and get rather fine grained
#' result from those.
#'
#' @param object [TODO] 
#' @param show.clusters.in.x Vector. Clusters(use cluster names) that are chosen to show in x-axis by users.
#' @param show.clusters.in.y Vector. Clusters(use cluster names) that are chosen to show in y-axis by users.
#' @param plot.cluster.group.method [TODO]
#' @param power.max.limit Numeric. Specify the upper limit of power, whose value is highly user-defined and data-dependent.
#' @param power.min.limit Numeric. Specify the lower limit of power, like \code{power.max.limit}.
#' @param plot.power.range [TODO] This will extend the range of actual value for plotting adjustment.
#' @param hide.power.label Logic. If TRUE, the label appended for power value will be hidden, otherwise, the label will be kept.
#' @param cnt.max.limit Numeric. Specify the upper limit of count, whose value is highly user-defined and data-dependent.
#' @param cnt.min.limit Numeric. Specify the lower limit of count, like \code{cnt.max.limit}.
#' @param plot.cnt.range [TODO] This will extend the range of actual value for plotting adjustment
#' @param hide.cnt.label Logic. If TRUE, the label appended for count value will be hidden, otherwise, the label will be kept.
#' @param nodes.size.range Numeric. It specifies the range of sizes of the nodes in the graph.
#' @param nodes.colour.seq Character. Given colours will be used to generate colour gradient for plotting.
#' @param nodes.colour.value.seq Numeric. If set NULL, evenly place each colour in \code{nodes.colour.seq} vector.
#' Otherwise, numeric values with the same length of \code{nodes.colour.seq} should be given to specify the positions corresponding to each colour.
#' See parameter \code{values} in \code{ggplot2::scale_colour_gradientn} for details.
#' @param label.power.options List. Options are hjust, vjust, nudge.x, nudge.y, size. See \code{ggplot2::geom_text} for details.
#' \itemize{
#'   \item hjust:   horizontal justification. Use either a string in ("left", "center", "right"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item vjust:   vertical justification. Use either a string in ("top", "middle", "bottom"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item nudge.x: horizontal adjustment to nudge labels by.
#'   \item nudge.y: vertical adjustment to nudge labels by.
#'   \item size:    label size.
#' }
#' @param label.cnt.options List. Options are like \code{label.power.options}.
#' @param plot.axis.x.name Character. X-axis name when plot is shown.
#' @param plot.axis.y.name Character. Y-axis name when plot is shown.
#' 
#' @details
#' This function uses some subset of interaction pairs, and does calculations in cluster level.
#' It calculates:
#' \itemize{
#'   \item count: count of interaction pairs.
#'   \item power: strength of interaction pairs, which is calculated by formula: SUM(abs(LogFC[geneA] * LogFC[geneB])), 
#'                as the function defines: \code{FullView.Evaluation.func.default}.
#' }
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @examples
#'
#'
#'
#' @import ggplot2
#'
#' @export
#'
GetResultFullView <- function(
	object,
	show.clusters.in.x = NULL,
	show.clusters.in.y = NULL,
	plot.cluster.group.method = "all",  # also support diagonal or diagonal-2
	power.max.limit = NULL,
	power.min.limit = NULL,
	plot.power.range = NULL, 
	hide.power.label = FALSE,
	cnt.max.limit = NULL,
	cnt.min.limit = NULL,
	plot.cnt.range = NULL, 
	hide.cnt.label = FALSE,
	nodes.size.range = c(1, 6),
	nodes.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
	nodes.colour.value.seq = c(0.0, 0.5, 1.0),
	label.power.options = list(hjust = "middle", vjust = "top", nudge.x = 0, nudge.y = -0.3, size = 2),
	label.cnt.options = list(hjust = "left", vjust = "center", nudge.x = 0.3, nudge.y = 0, size = 2),
	plot.axis.x.name = "clusters-x",
	plot.axis.y.name = "clusters-y"
) {
	interact.pairs.acted <- getFullViewResult(object)
	kClustersSplit <- getClusterSplit(object)
	# param user settings
	default.label.power.options <- list(hjust = "middle", vjust = "top", nudge.x = 0, nudge.y = -0.3, size = 2)
	default.label.cnt.options <- list(hjust = "left", vjust = "center", nudge.x = 0.3, nudge.y = 0, size = 2)
	user.label.power.options <- as.list(label.power.options)
	user.label.cnt.options <- as.list(label.cnt.options)
	# missing items in user settings
	missed.names.power.options <- setdiff(names(default.label.power.options), names(user.label.power.options))
	missed.names.cnt.options <- setdiff(names(default.label.cnt.options), names(user.label.cnt.options))
	# add missing items in user settings
	for (miss.p in missed.names.power.options) {
		user.label.power.options[[miss.p]] <- default.label.power.options[[miss.p]]
	}
	for (miss.cn in missed.names.cnt.options) {
		user.label.cnt.options[[miss.cn]] <- default.label.cnt.options[[miss.cn]]
	}
	# process
	fac.clusters <- interact.pairs.acted$list.clusters  # all clusters related
	fac.x.clusters <- fac.clusters$x.axis
	fac.y.clusters <- fac.clusters$y.axis
	if (!is.null(show.clusters.in.x) && length(show.clusters.in.x) != 0) {  # select part of clusters to be shown in x-axis
		fac.x.clusters <- intersect(show.clusters.in.x, fac.x.clusters)
		print(paste0("Reading clusters shown in X-axis: ", paste0(fac.x.clusters, collapse = ", "), "."))
	}
	if (!is.null(show.clusters.in.y) && length(show.clusters.in.y) != 0) {  # select part of clusters to be shown in y-axis
		fac.y.clusters <- intersect(show.clusters.in.y, fac.y.clusters)
		print(paste0("Reading clusters shown in Y-axis: ", paste0(fac.y.clusters, collapse = ", "), "."))
	}
	if (length(fac.x.clusters) < 1 || length(fac.y.clusters) < 1) {  # check
		stop("Error: X-Y plot needs at least one item in each axis!")
	}
	### data process for selected part of data
	inside.diagonal.selection <- function(vec.a, vec.b, identical.rm = FALSE) {
		if (sum(levels(factor(vec.a)) != levels(factor((vec.b)))) != 0) {
			stop("Diagonal method only allows clusters in x and y to be of same value range!")
		}
		vec.b <- rev(vec.b[match(vec.a, vec.b)])  # re-order to be the same order
		res.va <- res.vb <- NULL
		for (i in seq_along(vec.a)) {
			tmp.len <- ifelse(identical.rm == FALSE, length(vec.a) - i + 1, length(vec.a) - i)
			res.va <- c(res.va, rep(vec.a[i], times = tmp.len))
			res.vb <- c(res.vb, vec.b[seq_len(tmp.len)])
		}
		return(list(xsel = res.va, ysel = res.vb))
	}
	select.cluster.groups <- switch(plot.cluster.group.method,
		"all" = {
				col.x.data <- rep(fac.x.clusters, each = length(fac.y.clusters))
				col.y.data <- rep(fac.y.clusters, times = length(fac.x.clusters))
				list(xsel = col.x.data, ysel = col.y.data)
			},
		"diagonal" = inside.diagonal.selection(fac.x.clusters, fac.y.clusters),
		"diagonal-2" = inside.diagonal.selection(fac.x.clusters, fac.y.clusters, identical.rm = TRUE),
		stop("Undefined method for selecting cluster groups. Please re-check the given parameter: plot.cluster.group.method")
	)
	col.x.data <- select.cluster.groups$xsel
	col.y.data <- select.cluster.groups$ysel
	names.part.data <- paste(col.x.data, col.y.data, sep = kClustersSplit)
	ind.names.m <- match(names.part.data, interact.pairs.acted$name.allpairs)
	# ------
	# for cnt & power, use limit functions
	Limit.Max.inside <- function(x, eval.max) {
		res <- if (x > eval.max) eval.max else x
		res
	}
	Limit.Min.inside <- function(x, eval.min) {
		res <- if (x < eval.min) eval.min else x
		res
	}
	# ------
	# cnt
	cnt.xy.data <- interact.pairs.acted$cnt.allpairs
	cnt.part.data <- cnt.xy.data[ind.names.m]
	# cnt limit
	cnt.limit.part.data <- cnt.part.data
	if (!is.null(cnt.max.limit)) {
		cnt.limit.part.data <- sapply(cnt.limit.part.data, eval.max = cnt.max.limit, Limit.Max.inside)
	}
	if (!is.null(cnt.min.limit)) {
		cnt.limit.part.data <- sapply(cnt.limit.part.data, eval.min = cnt.min.limit, Limit.Min.inside)
	}
	# power
	power.xy.data <- interact.pairs.acted$strength.allpairs
	power.part.data <- power.xy.data[ind.names.m]
	# power limit
	power.limit.part.data <- power.part.data
	if (!is.null(power.max.limit)) {
		power.limit.part.data <- sapply(power.limit.part.data, eval.max = power.max.limit, Limit.Max.inside)
	}
	if (!is.null(power.min.limit)) {
		power.limit.part.data <- sapply(power.limit.part.data, eval.min = power.min.limit, Limit.Min.inside)    
	}
	# plot data preparation
	pairs.plot.db <- data.frame(
		pair.name = names.part.data, 
		x = col.x.data, 
		y = col.y.data,
		cnt.orig = as.integer(cnt.part.data),
		cnt.limit = as.integer(cnt.limit.part.data),
		power.orig = power.part.data,
		power.limit = power.limit.part.data,
		stringsAsFactors = FALSE
	)
	# check if it needs x-axis text rotation
	if.need.x.axis.rotate <- FALSE
	x.axis.text.len <- sapply(pairs.plot.db$x, function(x) {
		nchar(as.character(x))
	})
	if (max(x.axis.text.len) > 4) {
		if.need.x.axis.rotate <- TRUE
	}
	# limits range modification function
	inside.limits.adj.func <- function(user.limits.adj) {
		force(user.limits.adj)
		function(limits) {
			res.limits <- limits
			if (!is.null(user.limits.adj)) {
				res.limits <- user.limits.adj
			}
			return(res.limits)
		}
	}
	# plot process
	plot.res <- ggplot(pairs.plot.db, aes(x, y))
	plot.res <- plot.res + labs(x = plot.axis.x.name, y = plot.axis.y.name)
	plot.res <- plot.res + geom_point(aes(size = cnt.limit, colour = power.limit)) + 
												 scale_x_discrete(limits = unique(col.x.data), breaks = unique(col.x.data)) + 
												 scale_y_discrete(limits = unique(col.y.data), breaks = unique(col.y.data)) + 
												 scale_size(name = "Count", range = nodes.size.range,
													limits = inside.limits.adj.func(plot.cnt.range)) + 
												 scale_colour_gradientn(name = "Power", colours = nodes.colour.seq, values = nodes.colour.value.seq,
													limits = inside.limits.adj.func(plot.power.range), na.value = NA)
	# add labels for those out of range ( > cnt.max.limit or < cnt.min.limit)
	part.on.cnt.db <- data.frame(x = pairs.plot.db$x, y = pairs.plot.db$y, cnt.orig = pairs.plot.db$cnt.orig, stringsAsFactors = FALSE)
	part.on.cnt.ex.max.db <- NULL
	if (!is.null(cnt.max.limit)) {
		part.on.cnt.ex.max.db <- part.on.cnt.db[which(part.on.cnt.db$cnt.orig > cnt.max.limit), ]
	}
	part.on.cnt.ex.min.db <- NULL
	if (!is.null(cnt.min.limit)) {
		part.on.cnt.ex.min.db <- part.on.cnt.db[which(part.on.cnt.db$cnt.orig < cnt.min.limit), ]
	}
	part.on.cnt.limit.db <- rbind(part.on.cnt.ex.max.db, part.on.cnt.ex.min.db)
	if (!is.null(part.on.cnt.limit.db) && !hide.cnt.label) {
		plot.res <- plot.res + 
								geom_text(aes(label = cnt.orig), data = part.on.cnt.limit.db,
													hjust = user.label.cnt.options$hjust, vjust = user.label.cnt.options$vjust,
													nudge_x = user.label.cnt.options$nudge.x, nudge_y = user.label.cnt.options$nudge.y,
													size = user.label.cnt.options$size)
	}
	# add labels for those out of range (> power.max.limit or < power.min.limit)
	part.on.power.db <- data.frame(x = pairs.plot.db$x, y = pairs.plot.db$y, power.orig = pairs.plot.db$power.orig, stringsAsFactors = FALSE)
	part.on.power.ex.max.db <- NULL
	if (!is.null(power.max.limit)) {
		part.on.power.ex.max.db <- part.on.power.db[which(part.on.power.db$power.orig > power.max.limit), ]
	}
	part.on.power.ex.min.db <- NULL
	if (!is.null(power.min.limit)) {
		part.on.power.ex.min.db <- part.on.power.db[which(part.on.power.db$power.orig < power.min.limit), ]
	}
	part.on.power.limit.db <- rbind(part.on.power.ex.max.db, part.on.power.ex.min.db)
	if (!is.null(part.on.power.limit.db) && !hide.power.label) {
		plot.res <- plot.res + 
								geom_text(aes(label = round(power.orig, 2)), 
													data = part.on.power.limit.db,
													hjust = user.label.power.options$hjust, vjust = user.label.power.options$vjust,
													nudge_x = user.label.power.options$nudge.x, nudge_y = user.label.power.options$nudge.y,
													size = user.label.power.options$size)
	}
	plot.res <- plot.res + theme_classic()  # change theme
	# rotate labels in x-axis when needed
	if (if.need.x.axis.rotate) {
		plot.res <- plot.res + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	}
	## construct result tables
	# table raw.res
	res.raw.table <- data.frame(
		pair.name = pairs.plot.db$pair.name,
		x = pairs.plot.db$x,
		y = pairs.plot.db$y,
		cnt = pairs.plot.db$cnt.orig,
		power = pairs.plot.db$power.orig,
		stringsAsFactors = FALSE
	)
	# table cnt & power
	len.row <- length(fac.y.clusters)
	len.col <- length(fac.x.clusters)
	res.cnt.table <- matrix(data = c(0), nrow = len.row, ncol = len.col)
	res.power.table <- matrix(data = c(0.0), nrow = len.row, ncol = len.col)
	for (i in 1:len.col) {
		for (j in 1:len.row) {
			this.index <- len.row * (i - 1) + j
			res.cnt.table[j, i] <- pairs.plot.db$cnt.orig[this.index]
			res.power.table[j, i] <- round(pairs.plot.db$power.orig[this.index], 2)
		}
	}
	# set rownames & colnames
	rownames(res.cnt.table) <- as.character(fac.y.clusters)
	colnames(res.cnt.table) <- as.character(fac.x.clusters)
	rownames(res.power.table) <- as.character(fac.y.clusters)
	colnames(res.power.table) <- as.character(fac.x.clusters)
	#end# return 
	list(plot = plot.res, table = list(raw.res.table = res.raw.table, cnt.table = res.cnt.table, power.table = res.power.table))
}





# %%%%%%%%%%%%%%%%%%
# Other functions
# %%%%%%%%%%%%%%%%%%


FetchClusterGroups <- function(
	object,
	use.former = FALSE,
	cluster.former = character(),  # to X
	use.latter = FALSE,
	cluster.latter = character()  # to Y
) {
	this.fullview <- getFullViewResult(object)
	this.allowed.names <- this.fullview$name.allpairs
	this.former.allowed <- this.fullview$list.clusters$x.axis
	this.latter.allowed <- this.fullview$list.clusters$y.axis

	# pre-check
	if ((use.former == FALSE && length(cluster.former) > 0) ||
		  (use.former == TRUE && length(cluster.former) == 0)) {
		stop("To fetch the cluster relating to former, please set `use.former` = TRUE, and give valid clusters in `cluster.former`.")
	}
	if ((use.latter == FALSE && length(cluster.latter) > 0) ||
		  (use.latter == TRUE && length(cluster.latter) == 0)) {
		stop("To fetch the cluster relating to latter, please set `use.latter` = TRUE, and give valid clusters in `cluster.latter`.")
	}

	res.cg <- character()
	if (use.former == TRUE) {
		not.valid.cluster.X <- setdiff(cluster.former, this.former.allowed)
		if (length(not.valid.cluster.X) > 0) {
			warning("Given undefined clusters in X: ", paste0(not.valid.cluster.X, collapse = ", "), ". ")
		}
		cluster.former <- unique(intersect(cluster.former, this.former.allowed))
		former.res.list <- lapply(cluster.former, ref.names = this.allowed.names, function(x, ref.names) {
				grep(paste0("^", x), ref.names, value = TRUE)
			}
		)
		res.cg <- c(res.cg, as.character(unlist(former.res.list)))
	}

	if (use.latter == TRUE) {
		not.valid.cluster.Y <- setdiff(cluster.latter, this.latter.allowed)
		if (length(not.valid.cluster.Y) > 0) {
			warning("Given undefined clusters in Y: ", paste0(not.valid.cluster.Y, collapse = ", "), ". ")
		}
		cluster.latter <- unique(intersect(cluster.latter, this.latter.allowed))
		latter.res.list <- lapply(cluster.latter, ref.names = this.allowed.names, function(x, ref.names) {
				grep(paste0(x, "$"), ref.names, value = TRUE)
			}
		)
		res.cg <- c(res.cg, as.character(unlist(latter.res.list)))
	}

	res.cg <- unique(res.cg)
	print(paste("Fetch", length(res.cg), "cluster groups."))
	return(res.cg)
}

