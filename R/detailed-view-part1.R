


#' Collect gene pairs by expression change and action effects
#' 
#' @description
#' This function gets to collect gene pairs by their expression change (total 4 category) and 
#' action effects (total 4 category).
#'
#' @inheritParams Inside.DummyVEinfos
#' @param limits.exprs.change It specifies the range of expression change status that needs to be collected. 
#' @param limits.action.effects It specifies the range of action effects that need to be collected. 
#' 
#' @details
#' This function is to collect hierachical information inside the raw data, which is collecting gene pairs 
#' and grouping them by expression change status and action effects.
#'
#' @return
#' List of 2 layers. The first is the expression change status, and the list is length 4 and names of this list are 
#' "Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn". The second level of this list is the action effects, and every sublist gets
#' length 7 and its names are "A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B". Every single item is a table defined in 
#' class \code{data.frame}.
#'
#' @importFrom dplyr left_join
#'
#' @export
#'
CollectHierachyOnepairClusters <- function(
  VEinfos, 
  limits.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"), 
  limits.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B")
) {
  if (nrow(VEinfos$edges.infos) <= 0) {
    stop("Nither specific nor enough interaction pairs are given! Please recheck input!")
  }
  # predefining grouping rules
  group.proto.action.effects <- list(
    "A---B" = data.frame(),  # "A---B"  #1
    "A-->B" = data.frame(),  # "A-->B"  #2
    "A<--B" = data.frame(),  # "A<--B"  #3
    "A--|B" = data.frame(),  # "A--|B"  #4
    "A|--B" = data.frame(),  # "A|--B"  #5
    "A--oB" = data.frame(),  # "A--oB"  #6
    "Ao--B" = data.frame()   # "Ao--B"  #7
  )
  op.clustername <- VEinfos$cluster.name.A
  rv.clustername <- VEinfos$cluster.name.B
  group.proto.exprs.variations <- list(
    clusters.name = list(cluster.X = op.clustername, cluster.Y = rv.clustername),
    Xup.Yup = group.proto.action.effects,  # group.proto.action.effects, and the same for all below g.C*.D*
    Xup.Ydn = group.proto.action.effects, 
    Xdn.Yup = group.proto.action.effects,
    Xdn.Ydn = group.proto.action.effects
  )
  # process
  group.act.res <- group.proto.exprs.variations
  prog.bar.dgsa <- progress::progress_bar$new(total = length(limits.exprs.change) * length(limits.action.effects))
  prog.bar.dgsa$tick(0)
  for (i.exc in limits.exprs.change) {
    for (j.act in limits.action.effects) {
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
      # select subset by matching action.effect
      tmp.sub.1 <- VEinfos$edges.infos[which(VEinfos$edges.infos$action.effect == tmp.action.effect$action), ]
      tmp.m.cols <- c("UID", "ClusterName", "GeneName", "LogFC")
      tmp.sub.1 <- left_join(tmp.sub.1, VEinfos$vertices.infos[, tmp.m.cols], by = c("from" = "UID"))
      colnames(tmp.sub.1)[(ncol(tmp.sub.1) - length(tmp.m.cols) + 2) : ncol(tmp.sub.1)] <- paste("from", tmp.m.cols[2:length(tmp.m.cols)], sep = ".")
      tmp.sub.1 <- left_join(tmp.sub.1, VEinfos$vertices.infos[, tmp.m.cols], by = c("to" = "UID"))
      colnames(tmp.sub.1)[(ncol(tmp.sub.1) - length(tmp.m.cols) + 2) : ncol(tmp.sub.1)] <- paste("to", tmp.m.cols[2:length(tmp.m.cols)], sep = ".")
      tmp.sub.1 <- tmp.sub.1[, setdiff(colnames(tmp.sub.1), c("from", "to", "mode"))]
      # do select on action effect
      if (tmp.action.effect$direction == 1) {
        tmp.sub.1 <- tmp.sub.1[intersect(which(tmp.sub.1$from.ClusterName == op.clustername), which(tmp.sub.1$to.ClusterName == rv.clustername)), ]
      } else {
        if (tmp.action.effect$direction == -1) {
          tmp.sub.1 <- tmp.sub.1[intersect(which(tmp.sub.1$from.ClusterName == rv.clustername), which(tmp.sub.1$to.ClusterName == op.clustername)), ]
        } else {
          stop("Undefined direction ID!")
        }
      }
      
      # select subset by up.dn changes
      inside.updn.select <- function(df, UPDN.group, direction) {
        if (direction == 1) {  # direction change the actual mapping from X**.Y** and the action
          col.1 <- "from.LogFC"
          col.2 <- "to.LogFC"
        } else {
          col.1 <- "to.LogFC"
          col.2 <- "from.LogFC"
        }   
        res.df <- switch(UPDN.group,
          "Xup.Yup" = df[intersect(which(df[, col.1] > 0), which(df[, col.2] > 0)), ],
          "Xup.Ydn" = df[intersect(which(df[, col.1] > 0), which(df[, col.2] < 0)), ],
          "Xdn.Yup" = df[intersect(which(df[, col.1] < 0), which(df[, col.2] > 0)), ],
          "Xdn.Ydn" = df[intersect(which(df[, col.1] < 0), which(df[, col.2] < 0)), ],
          stop("Undefined expression change group in given parameter!")
        )
        return(res.df)
      }
      tmp.sub.1 <- inside.updn.select(tmp.sub.1, i.exc, tmp.action.effect$direction)
      # re-align the result
      tmp.sub.1 <- tmp.sub.1[, c("from.GeneName", "to.GeneName", "from.LogFC", "to.LogFC", "from.ClusterName", "to.ClusterName")]
      tmp.sub.1.conv <- tmp.sub.1[intersect(which(tmp.sub.1$from.ClusterName == op.clustername), which(tmp.sub.1$to.ClusterName == rv.clustername)), ]
      colnames(tmp.sub.1.conv) <- c("act.C.genename", "act.D.genename", "act.C.logfc", "act.D.logfc", "act.C.ClusterName", "act.D.ClusterName")
      tmp.sub.1.rev  <- tmp.sub.1[intersect(which(tmp.sub.1$from.ClusterName == rv.clustername), which(tmp.sub.1$to.ClusterName == op.clustername)), ]
      tmp.sub.1.rev  <- tmp.sub.1.rev[, ReverseOddEvenCols(6)]
      colnames(tmp.sub.1.rev) <- c("act.C.genename", "act.D.genename", "act.C.logfc", "act.D.logfc", "act.C.ClusterName", "act.D.ClusterName")
      tmp.sub.1 <- rbind(tmp.sub.1.conv, tmp.sub.1.rev)
      # add it to result
      if (nrow(tmp.sub.1) > 0) {
        tmp.sub.1 <- cbind(tmp.sub.1, UPDN.group = i.exc, ACT.TYPE = j.act, stringsAsFactors = FALSE)
        # do unique
        tmp.sub.1 <- DoPartUnique(tmp.sub.1, c(1:2))
        # add to result
        group.act.res[[i.exc]][[j.act]] <- tmp.sub.1
      }
      prog.bar.dgsa$tick()
    }
  }
  #end# return
  group.act.res
}





# [inside usage]
# This function is to draw plot for showing the result that contains different
# expression level changes and action effects.
# 
# @param onepair.dgsa A list. Return value of Inside.CollectHierachyOnepairClusters().
# @param show.exprs.change SEE GetResult.SummaryOnepairClusters for details.
# @param show.action.effects SEE GetResult.SummaryOnepairClusters for details.
# other @param The same as in GetResult.SummaryOnepairClusters. 
#
#
# @import ggplot2
# @import cowplot
#
Inside.PlotShowGrouping <- function(
  onepair.dgsa,
  show.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"),
  show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"),
  scale.fill.discrete,
  legend.key.size,
  legend.title,
  legend.text,
  legend.box.margin,
  caption.label.size,
  caption.label.x,
  caption.label.y,
  caption.label.hjust,
  caption.label.vjust,
  ...
) {
  # result allocation
  gsub.result <- vector("list", length = length(show.exprs.change))
  names(gsub.result) <- as.character(show.exprs.change)
  # process
  for (ix in show.exprs.change) {
    tmp.CD.dgsa <- onepair.dgsa[[ix]]
    tmp.CD.cnt <- 0  # total interact count
    for (i.type in show.action.effects) {
      tmp.CD.cnt <- tmp.CD.cnt + nrow(tmp.CD.dgsa[[i.type]])
    }
    gsub.result[[ix]] <- data.frame()
    for (i.type in show.action.effects) {  # split into their collections
      gsub.result[[ix]] <- rbind(gsub.result[[ix]],
        data.frame(dummy.x = ix,
          interact.cnt = nrow(tmp.CD.dgsa[[i.type]]),
          percent.cnt = nrow(tmp.CD.dgsa[[i.type]]) / tmp.CD.cnt,
          genes.action.effects = i.type))
    }
    gsub.result[[ix]]$genes.action.effects <- factor(gsub.result[[ix]]$genes.action.effects)
  }
  inside.PiePlot.shared <- function(data) {
    gplot.sgp <- ggplot(data, aes(x = dummy.x, y = interact.cnt, fill = genes.action.effects))
    gplot.sgp <- gplot.sgp + geom_col(position = "stack") +
      scale.fill.discrete + 
      coord_polar(theta = "y") + 
      theme(panel.background = element_blank(), plot.margin = margin(6, 0, 6, 0)) + 
      scale_x_discrete(breaks = NULL) + 
      scale_y_continuous(breaks = NULL)
    # return
    gplot.sgp
  }
  gplot.res.list <- list()
  for (ix in show.exprs.change) {
    tmp.gplot <- inside.PiePlot.shared(gsub.result[[ix]])
    gplot.res.list <- c(gplot.res.list, list(tmp.gplot))
  }
  # add unified legend
  res.legend <- get_legend(gplot.res.list[[1]] + 
    guides(fill = guide_legend(title = "Action Effect")) + 
    theme(legend.position = "right", 
      legend.key.size = legend.key.size,
      legend.title = legend.title,
      legend.text = legend.text,
      legend.box.margin = legend.box.margin)
    )
  # add additional plot modification
  for (ix in show.exprs.change) {
    ind.sel <- which(show.exprs.change == ix)
    gplot.res.list[[ind.sel]] <- gplot.res.list[[ind.sel]] + 
        xlab(NULL) + ylab(ix) + 
        theme(legend.position = "none")
  }
  gplot.res.list[[1]] <- gplot.res.list[[1]] + 
          xlab("interact pairs count") +  # add only 1 xlab
          theme(plot.margin = margin(6, 0, 6, 10))
  # grid the plots in one row
  res.plot.no.legend <- plot_grid(plotlist = gplot.res.list, align = "vh", nrow = 1)
  # result plot
  relative.width <- 0.3
  if (length(show.exprs.change) > 2) {
    relative.width <- 0.2
  }
  res.plot <- plot_grid(res.plot.no.legend, res.legend, nrow = 1, rel_widths = c(1, relative.width))
  # add joint caption
  caption.label <- paste0("---symbols---\nA,B: genes, X,Y: clusters.\nX: ", 
    onepair.dgsa$clusters.name[["cluster.X"]], ", Y: ", onepair.dgsa$clusters.name[["cluster.Y"]])
  caption.plot <- ggdraw() + 
    draw_label(caption.label, size = caption.label.size, 
      x = caption.label.x, hjust = caption.label.hjust,
      y = caption.label.y, vjust = caption.label.vjust) + 
    theme(plot.margin = margin(0, 0, 0, 0))
  final.plot <- plot_grid(res.plot, caption.plot, ncol = 1, rel_heights = c(1, 0.2))      
  #end# return
  final.plot
}






#' Get result of analysis on one pair
#'
#' @description
#' This function focuses on one interaction pair and its detailed information(expression changes, 
#' gene-gene action effects), and get result from it.
#'
#' @param onepair.dgsa List. Return value of \code{\link{CollectHierachyOnepairClusters}}.
#' @param show.exprs.change Character. Use exprssion level change of clusters to select part of data to be shown.
#' @param show.action.effects Character. Select some action effects to be put in the result.
#' @param scale.fill.discrete A ggplot object. It gives the colour strategy for plotting.
#' @param legend.key.size The size of keys in legend. It should be in unit format, see \code{?unit} for further help.
#' @param legend.title It sets the attributes of legend title, and should be \code{element_text()}.
#' @param legend.text It sets the attributes of legend annotation texts, and should be \code{element_text()}.
#' @param legend.box.margin It sets the margin of legend box, and should be \code{margin()}, see \code{?margin} for further help.
#' @param caption.label.size It sets the size of caption label.
#' @param caption.label.x It sets the position of caption label in x-axis.
#' @param caption.label.y It sets the position of caption label in y-axis.
#' @param caption.label.hjust It sets the alignment type of caption label in horizontal level, see \code{vignette("ggplot2-specs")}
#' for further help.
#' @param caption.label.vjust It sets the alignment type of caption label in vertical level, see \code{vignette("ggplot2-specs")}
#' for further help.
#'
#' @return A list. Use \code{Tool.ShowGraph()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @examples
#' # select only Xup.Yup
#' \dontrun{
#' GetResult.SummaryOnepairClusters(
#'   onepair.gmoc,
#'   show.exprs.change = c("Xup.Yup")
#' )
#' }
#'
#'
#'
#' @import ggplot2
#' @importFrom cowplot draw_label get_legend ggdraw plot_grid
#'
#' @export
#'
GetResult.SummaryOnepairClusters <- function(
  onepair.dgsa,
  show.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"),
  show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"),
  scale.fill.discrete = scale_fill_brewer(palette = "Set3"),
  legend.key.size = unit(6, "mm"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 6),
  caption.label.size = 12,
  caption.label.x = 0.5,
  caption.label.y = 1,
  caption.label.hjust = 0.5,
  caption.label.vjust = 1
) {
  # pre-result check
  notin.exprs.c <- setdiff(show.exprs.change, c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"))
  if (length(notin.exprs.c) > 0) {
    warning("Given @param show.exprs.change has undefined values: ", 
      paste0(notin.exprs.c, collapse = ", "), "\n and will be ignored!")
  }
  show.exprs.change <- setdiff(show.exprs.change, notin.exprs.c)
  if (length(show.exprs.change) == 0) {
    stop("Select available data less than 1!")
  }
  # pre-check 2
  notin.act.effect <- setdiff(show.action.effects, c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"))
  if (length(notin.act.effect) > 0) {
    warning("Given @param show.action.effects has undefined values: ",
      paste0(notin.act.effect, collapse = ", "), "\n and will be ignored!")
  }
  show.action.effects <- setdiff(show.action.effects, notin.act.effect)
  if (length(show.action.effects) == 0) {
    stop("Select action effects less than 1. Unable to fetch the result!")
  }
  ## result plot
  plot.res <- Inside.PlotShowGrouping(onepair.dgsa, 
    show.exprs.change = show.exprs.change, 
    show.action.effects = show.action.effects, 
    scale.fill.discrete = scale.fill.discrete,
    legend.key.size = legend.key.size,
    legend.title = legend.title,
    legend.text = legend.text,
    legend.box.margin = legend.box.margin,
    caption.label.size = caption.label.size,
    caption.label.x = caption.label.x,
    caption.label.y = caption.label.y,
    caption.label.hjust = caption.label.hjust,
    caption.label.vjust = caption.label.vjust
    )
  ## result table
  res.table.list <- list(
    Xup.Yup = data.frame(),
    Xup.Ydn = data.frame(),
    Xdn.Yup = data.frame(),
    Xdn.Ydn = data.frame()
  )
  this.cluster.X <- onepair.dgsa$clusters.name$cluster.X
  this.cluster.Y <- onepair.dgsa$clusters.name$cluster.Y
  for (i.exc in show.exprs.change) {
    this.exc <- onepair.dgsa[[i.exc]]
    for (i.type in show.action.effects) {
      if (nrow(this.exc[[i.type]]) == 0) {
        next
      }
      res.table.list[[i.exc]] <- rbind(res.table.list[[i.exc]],
        data.frame(cluster.X = this.cluster.X,
               gene.A = this.exc[[i.type]]$act.C.genename,
               logfc.A = this.exc[[i.type]]$act.C.logfc,
               cluster.Y = this.cluster.Y,
               gene.B = this.exc[[i.type]]$act.D.genename,
               logfc.B = this.exc[[i.type]]$act.D.logfc,
               action.effects = i.type,
               stringsAsFactors = FALSE
               )
        )
    }
  }
  res.final.table <- list()
  for (i.res in names(res.table.list)) {
    if (nrow(res.table.list[[i.res]]) > 0) {
      res.final.table[[i.res]] <- res.table.list[[i.res]]
    }
  }
  #end# return
  list(plot = plot.res, table = res.final.table)
}






# new added at 2021.03.15, consider deprecated at 2021.03.16
# for collect mode
# other mode means the given gene pairs has no mode
#
#
# @param use.prioritized.mode if TRUE, \code{limits.action.mode} will be prioritized one-after-one to get exclusive collection
#
# @import dplyr
# @import ggplot2
# @import cowplot
#
PiePlot.action.mode <- function(
  VEinfos,
  limits.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"),
  limits.action.mode = c(kpred.mode, "other"),
  use.prioritized.mode = FALSE
) {
  # pre-check
    # [TODO] check exprs.change
    # [TODO] check action mode

  # 
  this.cluster.A <- VEinfos$cluster.name.A
  this.cluster.B <- VEinfos$cluster.name.B
  this.vertices <- VEinfos$vertices.infos
  this.edges <- VEinfos$edges.infos
  # collect mode infos
  tmp.to.change.names <- c("ClusterName", "GeneName", "LogFC")
  this.pie.df <- left_join(this.edges, this.vertices[, c("UID", "ClusterName", "GeneName", "LogFC")], by = c("from" = "UID"))
  tmp.inds.from <- match(tmp.to.change.names, colnames(this.pie.df))
  colnames(this.pie.df)[tmp.inds.from] <- paste("from", tmp.to.change.names, sep = ".")
  this.pie.df <- left_join(this.pie.df, this.vertices[, c("UID", "ClusterName", "GeneName", "LogFC")], by = c("to" = "UID"))
  tmp.inds.to <- match(tmp.to.change.names, colnames(this.pie.df))
  colnames(this.pie.df)[tmp.inds.to] <- paste("to", tmp.to.change.names, sep = ".")
  # adjust the columns
  this.pie.df <- this.pie.df[, c("from.GeneName", "to.GeneName", "from.ClusterName", "to.ClusterName", "from.LogFC", "to.LogFC", "mode")]
  # use conv and rev to align the clusters
  tmp.conv.df <- this.pie.df[intersect(which(this.pie.df$from.ClusterName == this.cluster.A),
    which(this.pie.df$to.ClusterName == this.cluster.B)), ]
  tmp.rev.df <- this.pie.df[intersect(which(this.pie.df$to.ClusterName == this.cluster.A),
    which(this.pie.df$from.ClusterName == this.cluster.B)), ]
  tmp.rev.df <- tmp.rev.df[, c(ReverseOddEvenCols(6), 7:ncol(tmp.rev.df))]
  colnames(tmp.rev.df) <- colnames(tmp.conv.df)
  # merge the result
  this.pie.df <- rbind(tmp.conv.df, tmp.rev.df)
  this.pie.df <- DoPartUnique(this.pie.df, match(c("from.GeneName", "to.GeneName", "mode"), colnames(this.pie.df)))

  # inside plot function
  this.color.pal <- c("#FB8072", "#B3DE69", "#80B1D3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDB462", "#FCCDE5")
  names(this.color.pal) <- c(kpred.mode, "other")
  inside.PiePlot.for.mode <- function(data, color.palette) {
    #data <- data[order(data$exprs.change), ]
    data <- data[order(data$mode), ]
    data$cnt.percent <- data$cnt / sum(data$cnt) * 100
    data$label.y <- cumsum(data$cnt.percent) - 0.5 * data$cnt.percent
    data$label.y <- 100 - data$label.y
    #
    gplot.sgp <- ggplot(data, aes(x = ""))
    gplot.sgp <- gplot.sgp + geom_bar(aes(y = cnt.percent, fill = mode), stat = "identity", position = "stack") + 
      geom_text(aes(y = label.y, label = cnt)) +
      scale_fill_manual(values = color.palette) + 
      coord_polar(theta = "y") + 
      theme(panel.background = element_blank(), plot.margin = margin(6, 0, 6, 0)) + 
      scale_x_discrete(breaks = NULL) + 
      scale_y_continuous(breaks = NULL)
    # return
    gplot.sgp
  }

  # result prototype
  pie.plot.set <- vector(mode = "list", length = length(limits.exprs.change))
  names(pie.plot.set) <- limits.exprs.change
  #
  pie.exc.collect.Xup <- which(this.pie.df$from.LogFC > 0)
  pie.exc.collect.Xdn <- which(this.pie.df$from.LogFC <= 0)
  pie.exc.collect.Yup <- which(this.pie.df$to.LogFC > 0)
  pie.exc.collect.Ydn <- which(this.pie.df$to.LogFC <= 0)
  for (i.exc in limits.exprs.change) {
    this.exc.inds <- switch(i.exc,
      "Xup.Yup" = intersect(pie.exc.collect.Xup, pie.exc.collect.Yup),
      "Xup.Ydn" = intersect(pie.exc.collect.Xup, pie.exc.collect.Ydn),
      "Xdn.Yup" = intersect(pie.exc.collect.Xdn, pie.exc.collect.Yup),
      "Xdn.Ydn" = intersect(pie.exc.collect.Xdn, pie.exc.collect.Ydn),
      stop("Undefined expression changing status! Please check given parameter 'limits.exprs.change'!"))
    tmp.raw.df <- this.pie.df[this.exc.inds, ]
    # use prioritized mode to cover others
    if (use.prioritized.mode == TRUE) {
      tmp.gpairs <- paste(tmp.raw.df$from.GeneName, tmp.raw.df$to.GeneName, sep = kGenesSplit)
      tmp.inds.gpairs <- tapply(seq_along(tmp.gpairs), tmp.gpairs, function(x) { x })
      tmp.mode.cnt.raw <- lapply(tmp.inds.gpairs, 
        ref.df = tmp.raw.df, ref.mode.sets = limits.action.mode,
        function(x, ref.df, ref.mode.sets) {
          tmp.inds <- match(ref.mode.sets, ref.df[x, "mode"])
          tmp.inds <- tmp.inds[which(!is.na(tmp.inds))]
          c(ref.df[x[tmp.inds][1], "mode"], 
            paste(ref.df$from.GeneName[x[1]], ref.df$to.GeneName[x[1]], sep = "->"))
        })
      tmp.mode.cnt.raw <- as.character(unlist(tmp.mode.cnt.raw))
  print(tmp.mode.cnt.raw)
      tmp.mode.cnt.collect <- tapply(seq_along(tmp.mode.cnt.raw), tmp.mode.cnt.raw, length)
      if (length(tmp.mode.cnt.collect) > 0) {
        pie.prior.set <- data.frame(mode = names(tmp.mode.cnt.collect), cnt = tmp.mode.cnt.collect, exprs.change = i.exc)
        pie.prior.set <- pie.prior.set[which(pie.prior.set$mode %in% limits.action.mode), ]
  print(pie.prior.set)
        pie.plot.set[[i.exc]] <- inside.PiePlot.for.mode(pie.prior.set, this.color.pal)
      }
    } else {
      tmp.mode.res <- tapply(seq_len(nrow(tmp.raw.df)), tmp.raw.df$mode, length)
      if (length(tmp.mode.res) > 0) {
        pie.res.set <- data.frame(mode = names(tmp.mode.res), cnt = tmp.mode.res, exprs.change = i.exc)
        pie.res.set <- pie.res.set[which(pie.res.set$mode %in% limits.action.mode), ]
        pie.plot.set[[i.exc]] <- inside.PiePlot.for.mode(pie.res.set, this.color.pal)
      }
    }
  }

  ## plot pie
  # get legend
  pie.legend <- get_legend(pie.plot.set[[1]] + 
    guides(fill = guide_legend(title = "Action Mode")) + 
    theme(legend.position = "right",
      legend.key.size = unit(7, "mm"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.box.margin = margin(0, 0, 0, 6)))
  # add additional plot modification
  for (i.exc in limits.exprs.change) {
    if (!is.null(pie.plot.set[[i.exc]])) {
      pie.plot.set[[i.exc]] <- pie.plot.set[[i.exc]] + 
        xlab(NULL) + ylab(i.exc) + 
        theme(legend.position = "none")
    }
  }
  pie.plot.set[[1]] <- pie.plot.set[[1]] + 
          xlab("interact pairs count") +  # add only 1 xlab
          theme(plot.margin = margin(6, 0, 6, 10))
  # grid the plots in one row
  res.plot.no.legend <- plot_grid(plotlist = pie.plot.set, align = "vh", nrow = 1)
  # result plot
  relative.width <- 0.3
  if (length(limits.exprs.change) > 2) {
    relative.width <- 0.2
  }
  res.plot <- plot_grid(res.plot.no.legend, pie.legend, nrow = 1, rel_widths = c(1, relative.width))
  # add joint caption
  caption.label <- paste0("---symbols---\nX,Y: clusters that gene pairs from.\nX: ", 
    this.cluster.A, ", Y: ", this.cluster.B)
  caption.plot <- ggdraw() + 
    draw_label(caption.label, size = 12, 
      x = 0.5, hjust = 0.5,
      y = 0.8, vjust = 1) + 
   theme(plot.margin = margin(0, 0, 0, 0))
  final.plot <- plot_grid(res.plot, caption.plot, ncol = 1, rel_heights = c(1, 0.2))      
  
  #end# return
  return(list(plot = final.plot, table = this.pie.df))
}


