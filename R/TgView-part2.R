

#' Filter Important Gene Pairs in One Interaction
#'
#' @description
#' This function help filter and rank the important gene pairs by evaluating both power and confidence, 
#' which is calculated from 'LogFC' and 'PVal' respectively.
#'
#' @inheritParams InsideObjectInterCell
#' @param direction.X.to.Y Options are 'NULL', 'TRUE', 'FALSE'. It selects subset of data based on direction.
#'  The 'NULL' will keep 2-way interacting pairs, 'TRUE' keeps the X-to-Y pairs and 'FALSE' keeps the Y-to-X pairs. 
#'  See details for help.
#' @param bound.to.use It specifies the user specified bound for evaluation params. The values out of bound will be coerced to 
#'  either lower bound or upper bound. Default no bound is set, i.e. [-Inf, +Inf]
#' @param func.to.use The function used to further tranform the values, e.g. \code{log1p}. 
#' @param plot.X.to.Y The clusters drawn in x-axis and y-axis are in default aligned with the network analysis.
#'  If set FALSE, switch the clusters drawn in x-axis and y-axis.
#' @param axis.order.xy It determines how the gene names will be ordered in the axis when plotting. 
#'  Default is \code{AlphaBet}, options are \code{Power} and \code{Specificity}.
#' @param axis.order.xy.decreasing It determines whether the orders are of decreasing pattern or increasing pattern.
#' @param plot.font.size.base It gives the font size of texts such as labels and titles. 
#' @param nodes.colour.seq It specifies the colour sequence of the nodes.
#' @param nodes.colour.value.seq It is along with the param \code{nodes.colour.seq}, and changes the colour expansion.
#' @param nodes.size.range It specifies the size range of the nodes.
#' @param axis.text.x.pattern It defines the axis text style in x-axis. 
#'
#' @details
#' The meaning for 2 used values is:
#' \itemize{
#'   \item LogFC: the log fold change, which indicates the relative gene expression.
#'   \item PVal: the confidence of discovering the gene as differently expressed genes. 
#'               If it is generated from Seurat, it is orginally calculated by bonferroni correction.
#' }
#'
#' Illustration for \code{direction.X.to.Y}:
#' When running this function, gene pairs have been combined with their actions. As a result, gene pairs 
#' get to have direction for their action, e.g. IL6->IL6R, which means IL6 gets to activate IL6R, and the direction
#' should be IL6 to IL6R. Suppose IL6 is expressed by cell cluster X, IL6R is expressed by cell cluster Y, then 
#' IL6->IL6R will be reserved if \code{direction.X.to.Y} is set either 'NULL' or 'TRUE', but not 'FALSE'. 
#' The cluster X and Y is aligned with what users specified in \code{link{AnalyzeInterInFullView}}, and X corresponds to 
#' those clusters shown in x-axis and Y corresponds to those in y-axis. More closely, the X and Y are corresponding to 
#' \code{cluster.x} and \code{cluster.y} given by \code{\link{FetchInterOI}}.
#'
#' @return
#' List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#'
#' @export
#'
GetResultTgCrosstalk <- function(
  object,
  direction.X.to.Y = NULL,
  bound.to.use = list("Power" = c(-Inf, Inf), "Specificity" = c(-Inf, Inf)),
  func.to.use = list("Power" = function(x) {x}, "Specificity" = function(x) {1/(1+x)}),
  plot.X.to.Y = TRUE,
  axis.order.xy = c("AlphaBet", "AlphaBet"),  # how to order axis in final plot. Can also be one of colnames.to.cmp
  axis.order.xy.decreasing = c(TRUE, TRUE),  # order direction
  plot.font.size.base = 12, 
  nodes.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
  nodes.colour.value.seq = c(0.0, 0.5, 1.0),
  nodes.size.range = c(2, 8),
  axis.text.x.pattern = element_text(angle = 30, hjust = 1)
) {
  ITinfos <- getTgItInfo(object)
  VEinfos <- getTgVEInfo(object)
  this.musthave.colnames <- object@misc$musthave.colnames
  fgenes.remapped.all <- object@fgenes
  columns.to.use <- c("Power", "Specificity")

  # pre-check
  tmp.valid.boundnames <- CheckParamStd(names(bound.to.use), columns.to.use, "Bound names", stop.on.zero = FALSE)
  bound.to.use <- bound.to.use[tmp.valid.boundnames]
  tmp.valid.funcnames <- CheckParamStd(names(func.to.use), columns.to.use, "Tranform functions", stop.on.zero = FALSE)
  func.to.use <- func.to.use[tmp.valid.funcnames]

  # add default bounds
  tmp.not.in.bound <- setdiff(columns.to.use, tmp.valid.boundnames)
  if (length(tmp.not.in.bound) != 0) {
    for (tmp.i in tmp.not.in.bound) {
      bound.to.use <- c(bound.to.use, list(c(-Inf, Inf)))
      names(bound.to.use)[length(bound.to.use)] <- tmp.i
    }
  }
  # add default functions
  tmp.not.in.func <- setdiff(columns.to.use, tmp.valid.funcnames)
  if (length(tmp.not.in.func) != 0) {
    for (tmp.j in tmp.not.in.func) {
      func.to.use <- c(func.to.use, list(function(x) {x}))
      names(func.to.use)[length(func.to.use)] <- tmp.j
    }
  }

  # settle with columns names
  colnames.to.cmp <- paste0("inter.", columns.to.use)
  names(bound.to.use) <- paste0("inter.", names(bound.to.use))
  names(func.to.use) <- paste0("inter.", names(func.to.use))

  # check node color
  if (length(nodes.colour.seq) != length(nodes.colour.value.seq)) {
    stop("Colors and their gradient values should be of same length! Check parameter `nodes.colour.seq` and `nodes.colour.value.seq`.")
  }

  # check order settings
  allowed.order.xy <- c("AlphaBet", "Power", "Specificity")
  axis.order.xy[1] <- CheckParamStd(axis.order.xy[1], allowed.order.xy, "Order on X axis", stop.on.zero = TRUE)
  axis.order.xy[2] <- CheckParamStd(axis.order.xy[2], allowed.order.xy, "Order on Y axis", stop.on.zero = TRUE)


  # go with VEinfos
  act.A.clustername <- VEinfos$cluster.name.A
  act.B.clustername <- VEinfos$cluster.name.B
  edges.infos <- VEinfos$edges.infos
  vertices.infos <- VEinfos$vertices.infos
  # make one portable table and re-sync edges with vertices
  tmp.e2.col <- c("ClusterName", "GeneName", "LogFC")
  edges.fullinfos <- left_join(edges.infos, vertices.infos[, c("UID", tmp.e2.col)], by = c("from" = "UID"))
  tmp.ind.new3 <- match(c("ClusterName", "GeneName", "LogFC"), colnames(edges.fullinfos))
  colnames(edges.fullinfos)[tmp.ind.new3] <- c("from.cluster", "from.gene", "from.LogFC")
  edges.fullinfos <- left_join(edges.fullinfos, vertices.infos[, c("UID", tmp.e2.col)], by = c("to" = "UID"))
  tmp.ind.new3 <- match(c("ClusterName", "GeneName", "LogFC"), colnames(edges.fullinfos))
  colnames(edges.fullinfos)[tmp.ind.new3] <- c("to.cluster", "to.gene", "to.LogFC")
  edges.fullinfos <- edges.fullinfos[, paste0(rep(c("from.", "to."), times = length(tmp.e2.col)), rep(c("cluster", "gene", "LogFC"), each = 2))]
  
  # restrict to only one direction
  if (!is.null(direction.X.to.Y)) {
    if (direction.X.to.Y == TRUE) {
      tmp.inds <- intersect(which(edges.fullinfos[, "from.cluster"] == act.A.clustername), which(edges.fullinfos[, "to.cluster"] == act.B.clustername))
      edges.fullinfos <- edges.fullinfos[tmp.inds, ]
    } else {
      tmp.inds <- intersect(which(edges.fullinfos[, "to.cluster"] == act.A.clustername), which(edges.fullinfos[, "from.cluster"] == act.B.clustername))
      edges.fullinfos <- edges.fullinfos[tmp.inds, ]
    }    
  } else {
    # merge the bidirectional pairs, and from.cluster is set as act.A.clustername, to.cluster == act.B.clustername
    tmp.inds.conv <- intersect(which(edges.fullinfos[, "from.cluster"] == act.A.clustername), which(edges.fullinfos[, "to.cluster"] == act.B.clustername))
    tmp.edges.conv <- edges.fullinfos[tmp.inds.conv, ]
    tmp.inds.rev <- intersect(which(edges.fullinfos[, "to.cluster"] == act.A.clustername), which(edges.fullinfos[, "from.cluster"] == act.B.clustername))
    tmp.edges.rev <- edges.fullinfos[tmp.inds.rev, ReverseOddEvenCols(6)]
    colnames(tmp.edges.rev) <- colnames(tmp.edges.conv)  # treat reverse one as the conv one, that is to remove the direction meanings of gene pairs
    edges.fullinfos <- rbind(tmp.edges.conv, tmp.edges.rev)
  }
  edges.fullinfos <- unique(edges.fullinfos)

  # change X and Y when plotting (from.* will be shown in x-axis, to.* will be shown in y-axis)
  if (!is.null(plot.X.to.Y) && plot.X.to.Y == FALSE) {  # switch the x-axis and y-axis
    tmp.colname <- colnames(edges.fullinfos)
    edges.fullinfos <- edges.fullinfos[, ReverseOddEvenCols(6)]
    colnames(edges.fullinfos) <- tmp.colname
  }

  # use itinfo and add cluster info
  tmp.cluster.A <- ITinfos$clusters.name[1]
  tmp.cluster.B <- ITinfos$clusters.name[2]
  itused.infos <- ITinfos$bt.pairs
  itused.infos$inter.Cluster.A <- tmp.cluster.A
  itused.infos$inter.Cluster.B <- tmp.cluster.B
  # pack infos
  tmp.used.m.cols <- c(paste("inter", c("GeneName", "Cluster"), rep(c("A", "B"), each = 2), sep = "."), "inter.Power", "inter.Specificity")
  pack1.infos <- left_join(subset(edges.fullinfos, from.cluster == tmp.cluster.A & to.cluster == tmp.cluster.B), 
    itused.infos[, tmp.used.m.cols],
    by = c("from.cluster" = "inter.Cluster.A", "to.cluster" = "inter.Cluster.B",
      "from.gene" = "inter.GeneName.A", "to.gene" = "inter.GeneName.B"))
  pack2.infos <- left_join(subset(edges.fullinfos, from.cluster == tmp.cluster.B & to.cluster == tmp.cluster.A), 
    itused.infos[, tmp.used.m.cols],
    by = c("from.cluster" = "inter.Cluster.B", "to.cluster" = "inter.Cluster.A",
      "from.gene" = "inter.GeneName.B", "to.gene" = "inter.GeneName.A"))
  packed.infos <- rbind(pack1.infos, pack2.infos)
  
  ## apply user-defined formula
  for (i in seq_along(colnames.to.cmp)) {
    tmp.colname <- colnames.to.cmp[i]
    tmp.res <- as.numeric(func.to.use[[tmp.colname]](packed.infos[, tmp.colname]))
    # use max and min limit
    tmp.min <- bound.to.use[[tmp.colname]][1]  # min
    tmp.max <- bound.to.use[[tmp.colname]][2]  # max
    # avoid to be too far away from the real value
    tmp.min <- ifelse(tmp.min < min(tmp.res), min(tmp.res), tmp.min)
    tmp.max <- ifelse(tmp.max > max(tmp.res), max(tmp.res), tmp.max)
    # apply the bound
    if (is.numeric(tmp.min) && !is.na(tmp.min) && 
      is.numeric(tmp.max) && !is.na(tmp.max)) {
      tmp.res <- unlist(lapply(tmp.res, min = tmp.min, max = tmp.max, function(x, max, min) {
        ifelse(x < min, min, ifelse(x > max, max, x))
      }))
    }
    packed.infos[, tmp.colname] <- tmp.res
  }

  ## draw the graph
  if (nrow(packed.infos) == 0) {
    stop("No available data!")
  }
  x.axis.data <- packed.infos[, setdiff(colnames(packed.infos), grep("^to\\.", colnames(packed.infos)))]
  x.axis.names <- x.axis.data[, "from.gene"]
  y.axis.data <- packed.infos[, setdiff(colnames(packed.infos), grep("^from\\.", colnames(packed.infos)))] 
  y.axis.names <- y.axis.data[, "to.gene"]

  # set order for names of x-y-axises
  x.axis.names <- switch(axis.order.xy[1],
    "AlphaBet"    = x.axis.data[order(x.axis.data[, "from.gene"], decreasing = axis.order.xy.decreasing[1]), "from.gene"],
    "Power"       = x.axis.data[order(x.axis.data[, "inter.Power"], decreasing = axis.order.xy.decreasing[1]), "from.gene"],
    "Specificity" = x.axis.data[order(x.axis.data[, "inter.Specificity"], decreasing = axis.order.xy.decreasing[1]), "from.gene"],
    stop("Invalid order name is given upon x-axis!")
    )
  x.axis.names <- unique(x.axis.names)
  #
  y.axis.names <- switch(axis.order.xy[2],
    "AlphaBet"    = y.axis.data[order(y.axis.data[, "to.gene"], decreasing = axis.order.xy.decreasing[2]), "to.gene"],
    "Power"       = y.axis.data[order(y.axis.data[, "inter.Power"], decreasing = axis.order.xy.decreasing[2]), "to.gene"],
    "Specificity" = y.axis.data[order(y.axis.data[, "inter.Specificity"], decreasing = axis.order.xy.decreasing[2]), "to.gene"],
    stop("Invalid order name is given upon x-axis!")
    )
  y.axis.names <- unique(y.axis.names)


  ##
  gp.res <- ggplot(packed.infos, aes(x = from.gene, y = to.gene))
  tmp.sym.size <- sym("inter.Power")
  tmp.sym.colour <- sym("inter.Specificity")
  gp.res <- gp.res + 
      geom_point(aes(size = !!tmp.sym.size, colour = !!tmp.sym.colour)) + 
      scale_x_discrete(limits = x.axis.names, breaks = x.axis.names) + 
      scale_y_discrete(limits = y.axis.names, breaks = y.axis.names) + 
      scale_size(name = "Power", range = nodes.size.range) + 
      scale_colour_gradientn(name = "Specificity",
        colours = nodes.colour.seq, values = nodes.colour.value.seq)
  gp.res <- gp.res + 
      labs(x = packed.infos[1, "from.cluster"], y = packed.infos[1, "to.cluster"]) + 
      theme_half_open(font_size = plot.font.size.base) + 
      background_grid() + 
      theme(axis.text.x = axis.text.x.pattern)
  # draw the final graph
  return(list(plot = gp.res, table = packed.infos))
}
