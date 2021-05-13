
#' Analysis Composition of Actions
#' 
#' @description
#' This function gets to group gene pairs by their expression change, action mode and 
#' action effect.
#'
#' @inheritParams InsideObjectInterCell
#' @param limits.exprs.change Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'.
#'  It collects the gene pairs matching to given expression change.
#' @param limits.action.mode Options are listed in \code{kpred.action.mode}. It specifies 
#'  the range of used action mode.
#' @param limits.action.effect Options are listed in \code{kpred.action.effect}. It specifies 
#'  the range of used action effect.
#' 
#' @return A \code{InterCell} object.
#'
#' @importFrom dplyr left_join
#'
#' @export
#'
AnalyzeInterInAction <- function(
	object, 
	limits.exprs.change = kexprs.change, 
	limits.action.mode = kpred.action.mode,
	limits.action.effect = kpred.action.effect
) {
	VEinfos <- getTgVEInfo(object)

	limits.exprs.change <- CheckParamStd(limits.exprs.change, kexprs.change, "expression change", stop.on.zero = TRUE)
	limits.action.mode <- CheckParamStd(limits.action.mode, kpred.action.mode, "action mode", stop.on.zero = TRUE)
	limits.action.effect <- CheckParamStd(limits.action.effect, kpred.action.effect, "action effect", stop.on.zero = TRUE)

	# transform input to be kaction.id
	limits.action.effect <- TransActionEffectToActId(limits.action.effect)

	# predefining grouping rules - last level - mode
	group.proto.action.mode <- lapply(limits.action.mode, function(x) { data.frame() })
	names(group.proto.action.mode) <- limits.action.mode
	# predefining grouping rules - second level - action
	group.proto.action.effects <- lapply(limits.action.effect, use.sub = group.proto.action.mode,
		function(x, use.sub) {
			group.proto.action.mode
		})
	names(group.proto.action.effects) <- limits.action.effect
	# grouping - first level - expression change
	group.proto.exprs.variations <- lapply(limits.exprs.change, use.sub = group.proto.action.effects,
		function(x, use.sub) {
			group.proto.action.effects
		})
	names(group.proto.exprs.variations) <- limits.exprs.change

	# create result box
	orig.clusternames <- getOrigClusterNameTgVEInfo(object)
	group.act.res <- list(clusters.name = list(cluster.X = orig.clusternames$cluster.name.A, cluster.Y = orig.clusternames$cluster.name.B),
		result = group.proto.exprs.variations)

	## usage for process
	# cluster names
	op.clustername <- VEinfos$cluster.name.A
	rv.clustername <- VEinfos$cluster.name.B
	# all gene pairs
	allpairs <- VEinfos$edges.infos
	tmp.m.cols <- c("UID", "GeneName", "ClusterName", "LogFC")
	tmp.former.cols <- paste("former", tmp.m.cols[2:length(tmp.m.cols)], sep = ".")
	tmp.latter.cols <- paste("latter", tmp.m.cols[2:length(tmp.m.cols)], sep = ".")
	allpairs <- left_join(allpairs, VEinfos$vertices.infos[, tmp.m.cols], by = c("from" = "UID"))
	colnames(allpairs)[(ncol(allpairs) - length(tmp.m.cols) + 2) : ncol(allpairs)] <- tmp.former.cols
	allpairs <- left_join(allpairs, VEinfos$vertices.infos[, tmp.m.cols], by = c("to" = "UID"))
	colnames(allpairs)[(ncol(allpairs) - length(tmp.m.cols) + 2) : ncol(allpairs)] <- tmp.latter.cols
	allpairs <- allpairs[, setdiff(colnames(allpairs), c("from", "to"))]
	# reorder the columns
	tmp.paired.cols <- as.character(unlist(lapply(seq_along(tmp.former.cols), former.cols = tmp.former.cols, latter.cols = tmp.latter.cols, 
		function(x, former.cols, latter.cols) {
			c(former.cols[x], latter.cols[x])
		})))
	allpairs <- allpairs[, c(tmp.paired.cols, setdiff(colnames(allpairs), tmp.paired.cols))]

	# process
	prog.bar.dgsa <- progress::progress_bar$new(total = length(limits.exprs.change) * length(limits.action.effect))
	prog.bar.dgsa$tick(0)
	for (i.exc in limits.exprs.change) {
		for (j.act in limits.action.effect) {
			tmp.action.effect <- switch(j.act, 
				"A---B" = list(direction = 1, action = "undirected"),  # use one way direction, 1 or -1
				"A-->B" = list(direction = 1, action = "positive"),
				"A<--B" = list(direction = -1, action = "positive"),
				"A--|B" = list(direction = 1, action = "negative"),
				"A|--B" = list(direction = -1, action = "negative"),
				"A--oB" = list(direction = 1, action = "unspecified"),
				"Ao--B" = list(direction = -1, action = "unspecified"),
				stop("Undefined action effects in given parameter!")
				)

			## select subset 
			tmp.sub.1 <- cbind(allpairs, UPDN.group = i.exc, stringsAsFactors = FALSE)
			# select subset by matching action.effect and also re-align data, make former.* <-> cluster.A, latter.* <-> cluster.B
			if (tmp.action.effect$direction == 1) {
				tmp.sub.1 <- tmp.sub.1[intersect(which(tmp.sub.1$former.ClusterName == op.clustername), which(tmp.sub.1$latter.ClusterName == rv.clustername)), ]
			} else {
				if (tmp.action.effect$direction == -1) {
					tmp.sub.1 <- tmp.sub.1[intersect(which(tmp.sub.1$former.ClusterName == rv.clustername), which(tmp.sub.1$latter.ClusterName == op.clustername)), ]
					tmp.colnames <- colnames(tmp.sub.1)
					rev.some.cols <- colnames(tmp.sub.1)[ReverseOddEvenCols(length(tmp.paired.cols))]
					tmp.sub.1 <- tmp.sub.1[, c(rev.some.cols, setdiff(colnames(tmp.sub.1), tmp.paired.cols))]
					colnames(tmp.sub.1) <- tmp.colnames
				} else {
					stop("Undefined direction ID!")
				}
			}
			tmp.sub.1 <- tmp.sub.1[which(tmp.sub.1$action.effect == tmp.action.effect$action), ]
			
			# select subset by up.dn changes
			tmp.sub.1 <- switch(i.exc,
					"Xup.Yup" = tmp.sub.1[intersect(which(tmp.sub.1[, "former.LogFC"] > 0), which(tmp.sub.1[, "latter.LogFC"] > 0)), ],
					"Xup.Ydn" = tmp.sub.1[intersect(which(tmp.sub.1[, "former.LogFC"] > 0), which(tmp.sub.1[, "latter.LogFC"] < 0)), ],
					"Xdn.Yup" = tmp.sub.1[intersect(which(tmp.sub.1[, "former.LogFC"] < 0), which(tmp.sub.1[, "latter.LogFC"] > 0)), ],
					"Xdn.Ydn" = tmp.sub.1[intersect(which(tmp.sub.1[, "former.LogFC"] < 0), which(tmp.sub.1[, "latter.LogFC"] < 0)), ],
					stop("Undefined expression change group in given parameter!")
			)

			# to count on mode
			for (k.mode in limits.action.mode) {
				group.act.res$result[[i.exc]][[j.act]][[k.mode]] <- unique(tmp.sub.1[which(tmp.sub.1$mode == k.mode), ])
			}
			prog.bar.dgsa$tick()
		}
	}
	#end# return
	object <- setTgActionComp(object, group.act.res)
	return(object)
}





#' Get Result for Composition of Action Mode
#'
#' @description
#' This function focuses on one interaction pair and its detailed information(expression changes, 
#' gene-gene action mode), and get result from it.
#'
#' @inheritParams InsideObjectInterCell
#' @param limits.exprs.change Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'. It selects 
#'  the part of result to be shown.
#' @param limits.action.mode Options are listed in \code{kpred.action.mode}. It specifies 
#'  the range of shown action mode.
#' @param color.action.mode The colors to be used. It is one-to-one matching with \code{limits.action.mode}.
#' @param legend.title The content of legend title.
#' @param legend.title.size The size of legend title.
#' @param legend.key.size The size of keys in legend. 
#' @param legend.text.size The size of legend text.
#' @param legend.box.margin The margin of legend box, and should be given in \code{margin()}. See \code{?margin} for further help.
#' @param show.note If set TRUE, the note will be shown in the graph.
#' @param note.label.size Numeric. The size of note label.
#' @param note.lineheight Numeric. The lineheight of note.
#' @param note.postion.xy The position of note, should be 2 numbers that give the coordinate. Allowed range is (0~1, 0~1).
#' @param note.hjust The horizontal alignment of note. Allowed range is 0~1.
#' @param note.vjust The vertical alignment of note. Allowed range is 0~1.
#'
#' @return A list. Use \code{Tool.ShowGraph()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: NULL.
#'   \item grid.plot: the graph.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @import dplyr ggplot2
#' @importFrom cowplot draw_label get_legend ggdraw plot_grid
#'
#' @export
#'
GetResultPieActionMode <- function(
	object,
	limits.exprs.change = kexprs.change,
	limits.action.mode = kpred.action.mode,
	color.action.mode = kpred.color.mode,
	legend.title = "Action Mode",
	legend.title.size = element_text(size = 14),
	legend.key.size = unit(7, "mm"),
	legend.text.size = element_text(size = 12),
	legend.box.margin = margin(0, 0, 0, 6),
	show.note = TRUE,
	note.label.size = 3,
	note.lineheight = 1,
	note.postion.xy = c(0.5, 0.5),
	note.hjust = 0.5,
	note.vjust = 0.5
) {
	VEinfos <- getTgVEInfo(object)
	ActionComp <- getTgActionComp(object)
	kGenesSplit <- getGenePairSplit(object)

	# check expression change option
	limits.exprs.change <- CheckParamStd(limits.exprs.change, kexprs.change, opt.name = "expression change", stop.on.zero = TRUE)
	# check given action mode
	limits.action.mode <- CheckParamStd(limits.action.mode, kpred.action.mode, opt.name = "action mode", stop.on.zero = TRUE)

	# check color used in action mode
	if (is.character(color.action.mode) == TRUE) {
		if (length(color.action.mode) < length(limits.action.mode)) {
			stop("Given insufficient number of colors for action mode. Please check parameter `limits.action.mode` and `color.action.mode`. ")
		}
		names(color.action.mode) <- limits.action.mode
		color.action.mode <- scale_fill_manual(values = color.action.mode, breaks = limits.action.mode)
	} else {
		stop("Given invalid color in parameter `color.action.mode`!")		
	}
	
	## collect all action mode
	# data structure
	collect.mode.result <-  lapply(limits.exprs.change, function(x) { data.frame() })
	names(collect.mode.result) <- limits.exprs.change
	# collect process
	for (i.exc in limits.exprs.change) {
		this.exc.list <- ActionComp$result[[i.exc]]
		for (k.mode in limits.action.mode) {
			tmp.ret <- list()
			for (j.item in seq_along(this.exc.list)) {
				tmp.ret <- c(tmp.ret, list(this.exc.list[[j.item]][[k.mode]]))
			}
			collect.mode.result[[i.exc]] <- rbind(collect.mode.result[[i.exc]], bind_rows(tmp.ret))
		}
		collect.mode.result[[i.exc]] <- DoPartUnique(collect.mode.result[[i.exc]], 
			match(c("former.GeneName", "latter.GeneName", "mode"), colnames(collect.mode.result[[i.exc]])))
	}

	# plot data preparation
	plot.data.mode.result <- data.frame(exprs.change = character(),
		action.mode = character(),
		interact.cnt = integer(),
		stringsAsFactors = FALSE
		)
	for (i.exc in limits.exprs.change) {
		tmp.df <- collect.mode.result[[i.exc]]
		for (k.mode in limits.action.mode) {
			tmp.cnt <- nrow(tmp.df[which(tmp.df$mode == k.mode), ])
			if (tmp.cnt > 0) {
				plot.data.mode.result <- rbind(plot.data.mode.result,
					data.frame(exprs.change = i.exc, 
						action.mode = k.mode,
						interact.cnt = tmp.cnt,
						stringsAsFactors = FALSE
						)
				)
			}
		}
	}

	inside.PiePlot.mode <- function(data, color.palette) {
		if (nrow(data) == 0) {
			return(NULL)
		}
		data <- data[order(data$action.mode), ]
		data$cnt.percent <- data$interact.cnt / sum(data$interact.cnt) * 100
		data$label.y <- cumsum(data$cnt.percent) - 0.5 * data$cnt.percent
		data$label.y <- 100 - data$label.y
		#
		gplot.sgp <- ggplot(data, aes(x = exprs.change))
		gplot.sgp <- gplot.sgp + geom_bar(aes(y = cnt.percent, fill = action.mode), stat = "identity", position = "stack") +
			geom_text(aes(y = label.y, label = interact.cnt)) + 
			color.palette + 
			coord_polar(theta = "y") + 
			scale_x_discrete(breaks = NULL) + 
			scale_y_continuous(breaks = NULL) + 
			theme(panel.background = element_blank())
		# return
		gplot.sgp
	}

	plot.mode.list <- list()
	for (i.exc in limits.exprs.change) {
		tmp.plot <- inside.PiePlot.mode(plot.data.mode.result[which(plot.data.mode.result$exprs.change == i.exc), ],
				color.palette = color.action.mode)
		if (!is.null(tmp.plot)) {
			tmp.plot <- tmp.plot + xlab("") + ylab(i.exc) + theme(legend.position = "none")
			plot.mode.list <- c(plot.mode.list, list(tmp.plot))
		}
	}

	# add shared x lab
	if (length(plot.mode.list) > 0) {
		plot.mode.list[[1]] <- plot.mode.list[[1]] + xlab("pairs count")
	}

	# merged legend
	dummy.plot.for.legend <- ggplot(data = plot.data.mode.result) + 
		geom_bar(aes(x = action.mode, y = interact.cnt, fill = action.mode), stat = "identity", position = "stack") + 
		color.action.mode + 
		guides(fill = guide_legend(title = legend.title)) +
		theme(legend.position = "right", 
			legend.key.size = legend.key.size,
			legend.title = legend.title.size,
			legend.text = legend.text.size,
			legend.box.margin = legend.box.margin
		)	
	dummy.plot.obj <- ggplotGrob(dummy.plot.for.legend)
	dummy.the.legend <- dummy.plot.obj$grobs[[which(sapply(dummy.plot.obj$grobs, function(x) x$name) == "guide-box")]]

	res.plot <- arrangeGrob(grobs = c(plot.mode.list, list(dummy.the.legend)), nrow = 1, 
		widths = c(rep(4, times = length(plot.mode.list)), 1 * length(plot.mode.list))
	)

	if (show.note == TRUE) {
		note.text.it <- paste0("X,Y: clusters that gene pairs from.\nX: ", 
			ActionComp$clusters.name$cluster.X, ", Y: ", ActionComp$clusters.name$cluster.Y,
			"\nA: genes from cluster X, B: genes from cluster Y.")
		note.text.df <- data.frame(x = c(0, 0, 1, 1, note.postion.xy[1]),
			y = c(0, 1, 1, 0, note.postion.xy[2]),
			label = c("", "", "", "", note.text.it),
			stringsAsFactors = FALSE)

		text.ggplot <- ggplot(note.text.df) + 
			geom_text(aes(x = x, y = y, label = label),
				size = note.label.size,
				lineheight = note.lineheight,
				hjust = note.hjust,
				vjust = note.vjust) + 
			theme_void()

		res.plot <- arrangeGrob(res.plot, ggplotGrob(text.ggplot), ncol = 1,
			heights = c(4, 1))
	}

	#print(paste0("Generating plot for composition of action mode for gene pairs between ", ActionComp$clusters.name$cluster.X,
	#	"(X) and ", ActionComp$clusters.name$cluster.Y, "(Y). Genes from X are denoted as A, while those from Y as B."))
	
	#end# return
	return(list(plot = NULL, grid.plot = res.plot, table = c(list(plot.data.mode.result), collect.mode.result)))
}





#' Get Result for Composition of Action Effect
#'
#' @description
#' This function focuses on one interaction pair and its detailed information(expression changes, 
#' gene-gene action effect), and get result from it.
#'
#' @inheritParams InsideObjectInterCell
#' @param limits.exprs.change Options are 'Xup.Yup', 'Xup.Ydn', 'Xdn.Yup', 'Xdn.Ydn'. It selects 
#'  the part of result to be shown.
#' @param limits.ext.action.effect Options are listed in \code{kpred.ext.action.effect}. It specifies 
#'  the range of shown extended action effect.
#' @param color.ext.action.effect The colors to be used. It is one-to-one matching with \code{limits.ext.action.effect}.
#' @param legend.title The content of legend title.
#' @param legend.title.size The size of legend title.
#' @param legend.key.size The size of keys in legend. 
#' @param legend.text.size The size of legend text.
#' @param legend.box.margin The margin of legend box, and should be given in \code{margin()}. See \code{?margin} for further help.
#' @param show.note If set TRUE, the note will be shown in the graph.
#' @param note.label.size Numeric. The size of note label.
#' @param note.lineheight Numeric. The lineheight of note.
#' @param note.postion.xy The position of note, should be 2 numbers that give the coordinate. Allowed range is (0~1, 0~1).
#' @param note.hjust The horizontal alignment of note. Allowed range is 0~1.
#' @param note.vjust The vertical alignment of note. Allowed range is 0~1.
#'
#' @return A list. Use \code{Tool.ShowGraph()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: NULL.
#'   \item grid.plot: the graph.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @importFrom scales colour_ramp
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @importFrom cowplot draw_label get_legend ggdraw plot_grid
#'
#' @export
#'
GetResultPieActionEffect <- function(
	object,
	limits.exprs.change = kexprs.change,
	limits.ext.action.effect = kpred.ext.action.effect,
	color.ext.action.effect = kpred.color.ext.effect,
	legend.title = "Action Effect",
	legend.title.size = element_text(size = 14),
	legend.key.size = unit(7, "mm"),
	legend.text.size = element_text(size = 12),
	legend.box.margin = margin(0, 0, 0, 6),
	show.note = TRUE,
	note.label.size = 3,
	note.lineheight = 1,
	note.postion.xy = c(0.5, 0.5),
	note.hjust = 0.5,
	note.vjust = 0.5
) {
	VEinfos <- getTgVEInfo(object)
	ActionComp <- getTgActionComp(object)
	kGenesSplit <- getGenePairSplit(object)

	# check expression change option
	limits.exprs.change <- CheckParamStd(limits.exprs.change, kexprs.change, opt.name = "expression change", stop.on.zero = TRUE)
	# check given action effect
	limits.ext.action.effect <- CheckParamStd(limits.ext.action.effect, kpred.ext.action.effect, opt.name = "action effect", stop.on.zero = TRUE)

	# check color used in action effect
	if (is.character(color.ext.action.effect) == TRUE) {
		if (length(color.ext.action.effect) < length(limits.ext.action.effect)) {
			stop("Given insufficient number of colors for action effect. Please check parameter `limits.ext.action.effect` and `color.ext.action.effect`. ")
		}
		names(color.ext.action.effect) <- limits.ext.action.effect
		color.ext.action.effect <- scale_fill_manual(values = color.ext.action.effect, breaks = limits.ext.action.effect)
	} else {
		stop("Given invalid color in parameter `color.ext.action.effect`!")		
	}
	
	## collect all action mode
	# data structure
	collect.effect.result <-  lapply(limits.exprs.change, function(x) { data.frame() })
	names(collect.effect.result) <- limits.exprs.change
	# collect process
	for (i.exc in limits.exprs.change) {
		this.exc.list <- ActionComp$result[[i.exc]]
		for (j.act in limits.ext.action.effect) {
			tmp.ret <- list()
			for (k.item in seq_along(this.exc.list[[j.act]])) {
				tmp.df <- this.exc.list[[j.act]][[k.item]]
				tmp.df <- cbind(tmp.df, act.id = rep(j.act, times = nrow(tmp.df)), stringsAsFactors = FALSE)
				tmp.ret <- c(tmp.ret, list(tmp.df))
			}
			collect.effect.result[[i.exc]] <- rbind(collect.effect.result[[i.exc]], bind_rows(tmp.ret))
		}
		collect.effect.result[[i.exc]] <- DoPartUnique(collect.effect.result[[i.exc]], 
			match(c("former.GeneName", "latter.GeneName", "action.effect", "act.id"), colnames(collect.effect.result[[i.exc]])))
	}

	# check plot parameters
	if ((is.numeric(note.postion.xy) || is.integer(note.postion.xy)) 
		&& length(note.postion.xy) == 2) {
		if (note.postion.xy[[1]] > 1 || note.postion.xy[[2]] > 1
		 || note.postion.xy[[1]] < 0 || note.postion.xy[[2]] < 0) {
			stop("Please give 2 numeric number within 0~1 for parameter `note.postion.xy`.")
		}
	} else {
		stop("Please give 2 numeric number within 0~1 for parameter `note.postion.xy`.")
	}

	# plot data preparation
	plot.data.effect.result <- data.frame(exprs.change = character(),
		action.effect = character(),
		interact.cnt = integer(),
		stringsAsFactors = FALSE
		)
	for (i.exc in limits.exprs.change) {
		tmp.df <- collect.effect.result[[i.exc]]
		for (j.act in limits.ext.action.effect) {
			tmp.cnt <- nrow(tmp.df[which(tmp.df$act.id == j.act), ])
			if (tmp.cnt > 0) {
				plot.data.effect.result <- rbind(plot.data.effect.result,
					data.frame(exprs.change = i.exc, 
						action.effect = j.act,
						interact.cnt = tmp.cnt,
						stringsAsFactors = FALSE
						)
				)
			}
		}
	}

	inside.PiePlot.effect <- function(data, color.palette) {
		if (nrow(data) == 0) {
			return(NULL)
		}
		data <- data[order(data$action.effect), ]
		data$cnt.percent <- data$interact.cnt / sum(data$interact.cnt) * 100
		data$label.y <- cumsum(data$cnt.percent) - 0.5 * data$cnt.percent
		data$label.y <- 100 - data$label.y
		#
		gplot.sgp <- ggplot(data, aes(x = exprs.change))
		gplot.sgp <- gplot.sgp + geom_bar(aes(y = cnt.percent, fill = action.effect), stat = "identity", position = "stack") +
			geom_text(aes(y = label.y, label = interact.cnt)) + 
			color.palette + 
			coord_polar(theta = "y") + 
			scale_x_discrete(breaks = NULL) + 
			scale_y_continuous(breaks = NULL) + 
			theme(panel.background = element_blank())
		# return
		gplot.sgp
	}

	plot.effect.list <- list()
	for (i.exc in limits.exprs.change) {
		tmp.plot <- inside.PiePlot.effect(plot.data.effect.result[which(plot.data.effect.result$exprs.change == i.exc), ],
				color.palette = color.ext.action.effect)
		if (!is.null(tmp.plot)) {
			tmp.plot <- tmp.plot + xlab("") + ylab(i.exc) + theme(legend.position = "none")
			plot.effect.list <- c(plot.effect.list, list(tmp.plot))
		}
	}

	# add shared x lab
	if (length(plot.effect.list) > 0) {
		plot.effect.list[[1]] <- plot.effect.list[[1]] + xlab("pairs count")
	}

	# merged legend
	dummy.plot.for.legend <- ggplot(data = plot.data.effect.result) + 
		geom_bar(aes(x = action.effect, y = interact.cnt, fill = action.effect), stat = "identity", position = "stack") + 
		color.ext.action.effect + 
		guides(fill = guide_legend(title = legend.title)) +
		theme(legend.position = "right", 
			legend.key.size = legend.key.size,
			legend.title = legend.title.size,
			legend.text = legend.text.size,
			legend.box.margin = legend.box.margin
		)	
	dummy.plot.obj <- ggplotGrob(dummy.plot.for.legend)
	dummy.the.legend <- dummy.plot.obj$grobs[[which(sapply(dummy.plot.obj$grobs, function(x) x$name) == "guide-box")]]

	res.plot <- arrangeGrob(grobs = c(plot.effect.list, list(dummy.the.legend)), nrow = 1, 
		widths = c(rep(4, times = length(plot.effect.list)), 1 * length(plot.effect.list))
	)

	if (show.note == TRUE) {
		note.text.it <- paste0("X,Y: clusters that gene pairs from.\nX: ", 
			ActionComp$clusters.name$cluster.X, ", Y: ", ActionComp$clusters.name$cluster.Y,
			"\nA: genes from cluster X, B: genes from cluster Y.")
		note.text.df <- data.frame(x = c(0, 0, 1, 1, note.postion.xy[1]),
			y = c(0, 1, 1, 0, note.postion.xy[2]),
			label = c("", "", "", "", note.text.it),
			stringsAsFactors = FALSE)

		text.ggplot <- ggplot(note.text.df) + 
			geom_text(aes(x = x, y = y, label = label),
				size = note.label.size,
				lineheight = note.lineheight,
				hjust = note.hjust,
				vjust = note.vjust) + 
			theme_void()

		res.plot <- arrangeGrob(res.plot, ggplotGrob(text.ggplot), ncol = 1,
			heights = c(4, 1))
	}

	# return
	return(list(plot = NULL, grid.plot = res.plot, table = c(list(plot.data.effect.result), collect.effect.result)))
}


