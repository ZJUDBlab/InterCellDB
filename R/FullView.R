

#' Perform Network Analysis
#'
#' @description
#' This function works for network analysis, which calculates
#' count and power of interaction pairs among all given clusters.
#'
#' @inheritParams InsideObjectInterCell
#' @param sel.clusters.X Defining one part of interacting clusters. The options can be 
#'  found from \code{\link{ListAllClusters}}. See details for help.
#' @param sel.clusters.Y Defining the other part of interacting clusters. The options can be 
#'  found from \code{\link{ListAllClusters}}. See details for help.
#' @param sel.exprs.change Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'.
#'  It gives the corresponding expression change of every interacting gene pair.
#' @param sel.some.genes.X It gives the genes to be expected to be expressed in clusters given  
#'  by parameter \code{sel.clusters.X}. 
#' @param sel.some.genes.Y It gives the genes to be expected to be expressed in clusters given  
#'  by parameter \code{sel.clusters.Y}. 
#' @param sel.genes.option Options are 'intersect' or 'union'. 'intersect' strictly restricts gene pair to 
#'  have one gene partner in \code{sel.some.genes.X} and the other in \code{sel.some.genes.Y}. 'union' restricts 
#'  gene pair to have at least one gene either in \code{sel.some.genes.X} or \code{sel.some.genes.Y}.
#' @param sel.gene.pairs Directly specify the desired gene pairs. It should be given in standard table that is generated 
#'  by \code{\link{FormatCustomGenePairs}}. To note, it's strictly aligned to clusters, see details.
#' @param force.process It stops the program when no subset of genes are selected either by \code{sel.some.genes.X} and \code{sel.some.genes.Y}
#'  or by \code{sel.gene.pairs}, which may take long time to process. To force process, set this to TRUE.
#' @param verbose If set TRUE, the progress bar will be shown.
#'
#' @details
#' This function performs network analysis and calculates:
#' \itemize{
#'   \item count: count of interaction.
#'   \item power: strength of interaction, which is calculated by formula: SUM(abs(LogFC[gene.X] * LogFC[gene.Y])).
#' }
#'
#' \bold{sel.clusters.X} and \bold{sel.clusters.Y}: 
#' Interactions are defined by pairs of interacting clusters, or say cluster groups. 
#' For one interaction, e.g. cluster-Myeloid ~ cluster-T_cell, the fromer one will be 
#' restricted to clusters given in \code{sel.clusters.X} and the latter one will be 
#' retricted to clusters given in \code{sel.clusters.Y}.
#' The clusters given by \code{ListAllClusters} list all available clusters, and users can 
#' manually pick some or all of them. 
#'
#' \bold{sel.exprs.change}: 
#' The *up means the gene is up-regulated, i.e. its LogFC (log fold change) > 0.
#' The *dn means the gene is down-regulated. Take 'Xup.Ydn' for example,
#' it gets to explore the interaction between cluster X and cluster Y by taking up-regulated genes from X and 
#' down-regulated genes from Y. 
#'
#' \bold{sel.gene.pairs}: 
#' The gene pairs given in this parameter are strictly aligned to clusters given by 
#' parameter \code{sel.clusters.X} and \code{sel.clusters.Y}. After standardizing by \code{FormatCustomGenePairs},
#' the gene pairs are given in 4 columns, and 2 of them are named 'inter.GeneName.A', 'inter.GeneName.B'.
#' The corresponding relation is that genes listed in 'inter.GeneName.A' are used to compare with genes expressed by clusters given in \code{sel.clusters.X} 
#' and genes listed in 'inter.GeneName.B' are used to compare with genes expressed by clusters given in \code{sel.clusters.Y}.
#' For example, if the user gives C3~C3ar1 in \code{sel.gene.pairs}, Myeloid_cell and T_cell in \code{sel.clusters.X}, and 
#' fibroblast and B_cell in \code{sel.clusters.Y}, then C3 will be tested in Myeloid_cell and T_cell but not fibroblast and B_cell,
#' and C3ar1 will be only tested in fibroblast and B_cell. If the user need C3~C3ar1 to be tested in the opposite way, 
#' one more row C3ar1~C3 should be given in \code{sel.gene.pairs}. 
#'
#' @return A \code{InterCell} object.
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
	sel.gene.pairs = NULL,
	force.process = FALSE,
	verbose = TRUE
) {
	kGenesSplit <- getGenePairSplit(object)
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
			stop("Non-standardized table of gene pairs are given, please use function `FormatCustomGenePairs` first to get standardized one!")
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

	##  perform network analysis
	# if prog.bar shown
	if (verbose == TRUE) {
		prog.bar.fv <- progress::progress_bar$new(total = length(sel.clusters.X) * length(sel.clusters.Y))
		prog.bar.fv$tick(0)
	}
	# the used gene pairs 
	used.proc.pairs.db <- object@database@pairs.db
	# pre-slim for used database
	if (!is.null(sel.gene.pairs)) {
		tmp.ref.gpairs <- paste0(used.proc.pairs.db$inter.GeneName.A, used.proc.pairs.db$inter.GeneName.B, sep = kGenesSplit)
		#
		tmp.used.conv <- paste0(sel.gene.pairs$inter.GeneName.A, sel.gene.pairs$inter.GeneName.B, sep = kGenesSplit)
		tmp.used.rev  <- paste0(sel.gene.pairs$inter.GeneName.B, sel.gene.pairs$inter.GeneName.A, sep = kGenesSplit)
		# set the pairs database
		used.proc.pairs.db <- sel.gene.pairs[union(which(tmp.used.conv %in% tmp.ref.gpairs), which(tmp.used.rev %in% tmp.ref.gpairs)), ]
	} else {  # slim for some genes
		tmp.used.genes <- c(sel.some.genes.X, sel.some.genes.Y)
		if (length(tmp.used.genes) > 0) {
			used.proc.pairs.db <- used.proc.pairs.db[union(which(used.proc.pairs.db$inter.GeneName.A %in% tmp.used.genes), 
				which(used.proc.pairs.db$inter.GeneName.B %in% tmp.used.genes)), ]
		}
	}

	# the process
	for (ix in sel.clusters.X) {
		for (jy in sel.clusters.Y) {
			interact.name <- paste0(ix, object@tool.vars$cluster.split, jy)
			fgenes.t.X <- object@fgenes[which(object@fgenes$cluster == ix), ]
			fgenes.t.Y <- object@fgenes[which(object@fgenes$cluster == jy), ]
			genes.X <- fgenes.t.X$gene
			genes.Y <- fgenes.t.Y$gene

			# generate all possible gene pairs
			ref.GeneA <- used.proc.pairs.db$inter.GeneName.A
			ref.GeneB <- used.proc.pairs.db$inter.GeneName.B
			# conv
			inds.gpairs.conv <- intersect(which(ref.GeneA %in% genes.X), which(ref.GeneB %in% genes.Y))
			gpairs.result <- used.proc.pairs.db[inds.gpairs.conv, ]
			
			if (is.null(sel.gene.pairs)) {  # using reference database, then rev. ones should be collected
				# rev
				inds.gpairs.rev <- intersect(which(ref.GeneB %in% genes.X), which(ref.GeneA %in% genes.Y))
				cols.rev <- ReverseOddEvenCols(std.index.colname.end.dual)
				gpairs.rev <- used.proc.pairs.db[inds.gpairs.rev, c(cols.rev, setdiff(seq_len(ncol(used.proc.pairs.db)), cols.rev))]
				colnames(gpairs.rev) <- colnames(used.proc.pairs.db)
				# merge result
				gpairs.result <- rbind(gpairs.result, gpairs.rev)
			}
			
			# prioritized usage of selection to process
			if (is.null(sel.gene.pairs)) {  # use selected genes
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





#' Get Result for Network Analysis
#' 
#' @description
#' This function summarizes the result of network analysis in graph and table, and gives clues to 
#' most active intercellular communications.
#'
#' @inheritParams InsideObjectInterCell
#' @param show.clusters.in.x Clusters chosen to show in x-axis, corresponding to paramter \code{sel.clusters.X} in \code{\link{AnalyzeInterInFullView}}.
#' @param show.clusters.in.y Clusters chosen to show in y-axis, corresponding to paramter \code{sel.clusters.Y} in \code{\link{AnalyzeInterInFullView}}.
#' @param sel.cluster.group.method Options are 'all', 'diagonal' and 'diagonal-2'. It controls 
#'  'all' is recommended for most cases, and 
#'  'diagonal*' is recommended when same restrictions are applied to clusters in X and Y.
#' @param power.max.limit One number that specifies the upper limit of power, and power larger than this value will be coerced to this value.
#' @param power.min.limit One number that specifies the lower limit of power, and power lower than this value will be coerced to this value.
#' @param cnt.max.limit One number that specifies the upper limit of count, and larger ones will be coerced to this.
#' @param cnt.min.limit Numeric. Specify the lower limit of count, and lower ones will be coerced to this.
#' @param plot.power.range This will extend the plotting range beyond actual value range for power.
#' @param hide.power.label If set TRUE, the label appended for power value will be hidden, otherwise, the label will be kept.
#' @param plot.cnt.range This will extend the plotting range beyond actual value range for count.
#' @param hide.cnt.label If set TRUE, the label appended for count value will be hidden, otherwise, the label will be kept.
#' @param nodes.size.range 2 numbers that specifies the range of sizes of the nodes shown in the graph.
#' @param nodes.colour.seq Given colours will be used to generate colour gradient for plotting.
#' @param nodes.colour.value.seq Numeric. If set NULL, colours given in \code{nodes.colour.seq} will be evenly placed.
#' Otherwise, numeric values with the same length of \code{nodes.colour.seq} should be given to specify the positions corresponding to each colour.
#' See parameter \code{values} in \code{ggplot2::scale_colour_gradientn} for details.
#' @param label.power.options List of 'power' label appearance control parameters, including hjust, vjust, nudge.x, nudge.y, size. See details for help.
#' @param label.cnt.options List of label appearance control parameters like those in \code{label.power.options}.
#' @param plot.axis.x.name X-axis name when plotting.
#' @param plot.axis.y.name Y-axis name when plotting.
#' 
#' @details
#' \code{label.power.options} and \code{label.cnt.options}: 
#' The meanings for each item are recommended to refer to \code{ggplot2::geom_text} for details.
#' \itemize{
#'   \item hjust:   horizontal justification. Use either a string in ("left", "center", "right"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item vjust:   vertical justification. Use either a string in ("top", "middle", "bottom"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item nudge.x: horizontal adjustment to nudge labels by.
#'   \item nudge.y: vertical adjustment to nudge labels by.
#'   \item size:    label size.
#' }
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @import ggplot2
#'
#' @export
#'
GetResultFullView <- function(
	object,
	show.clusters.in.x = NULL,
	show.clusters.in.y = NULL,
	sel.cluster.group.method = "all",
	power.max.limit = NULL,
	power.min.limit = NULL,
	cnt.max.limit = NULL,
	cnt.min.limit = NULL,
	plot.power.range = NULL, 
	hide.power.label = FALSE,
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
	select.cluster.groups <- switch(sel.cluster.group.method,
		"all" = {
				col.x.data <- rep(fac.x.clusters, each = length(fac.y.clusters))
				col.y.data <- rep(fac.y.clusters, times = length(fac.x.clusters))
				list(xsel = col.x.data, ysel = col.y.data)
			},
		"diagonal" = inside.diagonal.selection(fac.x.clusters, fac.y.clusters),
		"diagonal-2" = inside.diagonal.selection(fac.x.clusters, fac.y.clusters, identical.rm = TRUE),
		stop("Undefined method for selecting cluster groups. Please re-check the given parameter: sel.cluster.group.method")
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





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Other functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' List Cluster Groups
#'
#' This function is to help on listing all or some of cluster groups (i.e. interacting clusters).
#'
#' @inheritParams InsideObjectInterCell
#' @param use.former It controls whether to list cluster groups relating to one or some clusters(.X) given in 
#'  parameter \code{sel.clusters.X} of \code{link{AnalyzeInterInFullView}}.
#' @param cluster.former The name of clusters, which is used find all X~? interactions (X is cluters given in these parameter). 
#' @param use.latter It controls whether to list cluster groups relating to one or some clusters(.Y) given in 
#'  parameter \code{sel.clusters.Y} of \code{link{AnalyzeInterInFullView}}.
#' @param cluster.latter The name of clusters, which is used find all ?~Y interactions (Y is cluters given in these parameter). 
#'
#' @return Character. Names of cluster groups.
#'
#' @examples
#' \dontrun{
#'  ListClusterGroups(object, use.former = TRUE, cluster.former = "T_cell")
#'  # the output of this will be "T_cell~Macrophage", "T_cell~B_cell" and things like that, depending on available clusters.
#' 
#'  ListClusterGroups(object, use.latter = TRUE, cluster.latter = "T_cell")
#'  # the output of this will be "Macrophage~T_cell", "B_cell~T_cell" and things like that, depending on available clusters.
#' }
#'
#' @export
#' 
ListClusterGroups <- function(
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

	# if nothing is specified return all cluster groups
	if (use.former == FALSE && length(cluster.former) == 0 &&
		use.latter == FALSE && length(cluster.latter) == 0) {
		res.cg <- this.allowed.names
	}

	res.cg <- unique(res.cg)
	return(res.cg)
}

