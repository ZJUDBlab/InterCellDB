

#' Evaluation gene pairs for one pair of clusters
#'
#' @description
#' This function evaluates the importance of some specific gene pairs of one pair of clusters, 
#' and the default evaluation params are LogFC and PValAdj.
#'
#' @param object [TODO]
#' @param plot.x.to.y [TODO]
#' @param colnames.to.cmp Character. The colnames to be used as evaluation params, currently only 2 params are supported.
#' The 1st one will be plotted differently by different size of nodes, and 2nd one will be different by colour of nodes.
#' @param range.to.use List. It specifies the user specified ranges for evaluation params.
#' @param axis.order.xy Character. It determines how the gene names will be ordered in the axis when plotting.
#' @param axis.order.xy.decreasing Logic. It determines whether the orders are of decreasing pattern or increasing pattern.
#' @param plot.font.size.base Numeric. It defines the font size of texts such as labels and titles. 
#' @param nodes.colour.seq Character. It specifies the colour sequence of the nodes.
#' @param nodes.colour.value.seq Numeric. It is along with the param \code{nodes.colour.seq}, and changes the colour expansion.
#' @param nodes.size.range Numeric. It specifies the size range of the nodes.
#' @param axis.text.x.pattern It defines the axis text style in x-axis. 
#'
#'
#' @details
#' This function uses user-selected gene pairs, and uses evalution params to calculate the relative importance of each gene pair.
#' It calculates in default settings:
#' \itemize{
#'   \item LogFC: the log of fold changes, which indicates the relative gene expression.
#'   \item PValAdj: the confidence of discovering the gene as differently expressed genes. In Seurat, it uses bonferroni correction.
#' }
#'
#'
#'
#' @return
#' List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
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
  plot.x.to.y = TRUE,
  colnames.to.cmp = c("LogFC", "PVal"),
  range.to.use = list("LogFC" = c(-Inf, Inf), "PVal" = c(-Inf, Inf)),
  axis.order.xy = c("AlphaBet", "AlphaBet"),  # how to order axis in final plot. Can also be one of colnames.to.cmp
  axis.order.xy.decreasing = c(TRUE, TRUE),  # order direction
  plot.font.size.base = 12, 
  nodes.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
  nodes.colour.value.seq = c(0.0, 0.5, 1.0),
  nodes.size.range = c(2, 8),
  axis.text.x.pattern = element_text(angle = 30, hjust = 1)
) {
  VEinfos <- getTgVEInfo(object)
  fgenes.remapped.all <- object@fgenes
  formula.to.use <- object@formulae[c("TG.LOGFC", "TG.PVAL")]
  # pre-check
  if (!all(colnames.to.cmp %in% c("LogFC", "PVal"))) {
    stop("Invalid colnames detected, only 'LogFC' and 'PVal' are supported yet.")
  }
  tmp.range.not.valid <- setdiff(names(range.to.use), colnames.to.cmp)
  if (length(tmp.range.not.valid) != 0) {
    warning("Find UNVALID names in given range: ", paste0(tmp.range.not.valid, collapse = ", "))
  }
  tmp.not.in.range <- setdiff(colnames.to.cmp, colnames.to.cmp[which(colnames.to.cmp %in% names(range.to.use))])
  if (length(tmp.not.in.range) != 0) {
    for (tmp.i in tmp.not.in.range) {
      range.to.use <- c(range.to.use, list(c(-Inf, Inf)))
      names(range.to.use)[length(range.to.use)] <- tmp.i
    }
    warning("Auto add ranges on: ", paste0(tmp.not.in.range, collapse = ", "))
  }
  

  # VEinfos
  act.A.clustername <- VEinfos$cluster.name.A
  act.B.clustername <- VEinfos$cluster.name.B
  edges.infos <- VEinfos$edges.infos
  vertices.infos <- VEinfos$vertices.infos

  ## make one portable table
  # re-sync edges with vertices
  tmp.e2.col <- c("ClusterName", "GeneName")
  edges.fullinfos <- left_join(edges.infos, vertices.infos[, c("UID", tmp.e2.col)], by = c("from" = "UID"))
  tmp.ind.new3 <- match(c("ClusterName", "GeneName"), colnames(edges.fullinfos))
  colnames(edges.fullinfos)[tmp.ind.new3] <- c("from.cluster", "from.gene")
  edges.fullinfos <- left_join(edges.fullinfos, vertices.infos[, c("UID", tmp.e2.col)], by = c("to" = "UID"))
  tmp.ind.new3 <- match(c("ClusterName", "GeneName"), colnames(edges.fullinfos))
  colnames(edges.fullinfos)[tmp.ind.new3] <- c("to.cluster", "to.gene")
  edges.fullinfos <- edges.fullinfos[, c("from.cluster", "to.cluster", "from.gene", "to.gene")]
  # restrict to only one direction
  if (!is.null(plot.x.to.y)) {
    if (plot.x.to.y == TRUE) {
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
    tmp.edges.rev <- edges.fullinfos[tmp.inds.rev, ReverseOddEvenCols(4)]
    colnames(tmp.edges.rev) <- colnames(tmp.edges.conv)  # treat reverse one as the conv one, that is to remove the direction meanings of gene pairs
    edges.fullinfos <- rbind(tmp.edges.conv, tmp.edges.rev)
  }
  edges.fullinfos <- unique(edges.fullinfos)
  # compact on vertices.infos
  vertices.selinfos <- vertices.infos[, c(tmp.e2.col, colnames.to.cmp)]
  vertices.selinfos <- DoPartUnique(vertices.selinfos, match(tmp.e2.col, colnames(vertices.selinfos)))
  # merge all
  packed.infos <- left_join(edges.fullinfos, vertices.selinfos, by = c("from.cluster" = "ClusterName", "from.gene" = "GeneName"))
  tmp.add2.col <- match(colnames.to.cmp, colnames(packed.infos))
  colnames(packed.infos)[tmp.add2.col] <- paste("from", colnames.to.cmp, sep = ".")
  packed.infos <- left_join(packed.infos, vertices.selinfos, by = c("to.cluster" = "ClusterName", "to.gene" = "GeneName"))
  tmp.add2.col <- match(colnames.to.cmp, colnames(packed.infos))
  colnames(packed.infos)[tmp.add2.col] <- paste("to", colnames.to.cmp, sep = ".")
  ## apply user-defined formula
  for (i in 1:length(colnames.to.cmp)) {
    tmp.colname <- colnames.to.cmp[i]
    tmp.sel.cols <- paste(c("from", "to"), tmp.colname, sep = ".")
    tmp.res <- formula.to.use[[tmp.colname]](packed.infos[, tmp.sel.cols[1]], packed.infos[, tmp.sel.cols[2]])
    tmp.res <- as.numeric(tmp.res)
    tmp.res.coln <- paste("res", tmp.colname, sep = ".")
    # use max and min limit
    tmp.min <- range.to.use[[tmp.colname]][1]  # min
    tmp.max <- range.to.use[[tmp.colname]][2]  # max
    # avoid to be too far away from the real value
    tmp.min <- ifelse(tmp.min < min(tmp.res), min(tmp.res), tmp.min)
    tmp.max <- ifelse(tmp.max > max(tmp.res), max(tmp.res), tmp.max)
    # apply the range
    if (is.numeric(tmp.min) && !is.na(tmp.min) && 
      is.numeric(tmp.max) && !is.na(tmp.max)) {
      tmp.res <- unlist(lapply(tmp.res, min = tmp.min, max = tmp.max, function(x, max, min) {
        ifelse(x < min, min, ifelse(x > max, max, x))
      }))
    }
    packed.infos[, tmp.res.coln] <- tmp.res
  }

  ## draw the graph
  if (nrow(packed.infos) == 0) {
    stop("No available data!")
  }
  x.axis.data <- packed.infos[, grep("^from\\.", colnames(packed.infos))]
  x.axis.names <- x.axis.data[, "from.gene"]
  y.axis.data <- packed.infos[, grep("^to\\.", colnames(packed.infos))]
  y.axis.names <- y.axis.data[, "to.gene"]

  # set order for names of x-y-axises 
  if (axis.order.xy[1] == "AlphaBet") {
    x.axis.names <- x.axis.data[order(x.axis.data[, "from.gene"], decreasing = axis.order.xy.decreasing[1]), "from.gene"]
  } else {
    if (axis.order.xy[1] %in% colnames.to.cmp) {
      x.axis.names <- x.axis.data[order(x.axis.data[, paste("from", axis.order.xy[1], sep = ".")], decreasing = axis.order.xy.decreasing[1]), "from.gene"]
    } else {
      stop("Invalid order name is given upon x-axis!")
    }
  }
  x.axis.names <- unique(x.axis.names)
  if (axis.order.xy[2] == "AlphaBet") {
    y.axis.names <- y.axis.data[order(y.axis.data[, "to.gene"], decreasing = axis.order.xy.decreasing[2]), "to.gene"]
  } else {
    if (axis.order.xy[2] %in% colnames.to.cmp) {
      y.axis.names <- y.axis.data[order(y.axis.data[, paste("to", axis.order.xy[2], sep = ".")], decreasing = axis.order.xy.decreasing[2]), "to.gene"]
    } else {
      stop("Invalid order name is given upon y-axis!")
    }
  }
  y.axis.names <- unique(y.axis.names)
  ##
  gp.res <- ggplot(packed.infos, aes(x = from.gene, y = to.gene))
  tmp.sym.size <- sym(paste("res", colnames.to.cmp[1], sep = "."))
  tmp.sym.colour <- sym(paste("res", colnames.to.cmp[2], sep = "."))
  gp.res <- gp.res + 
      geom_point(aes(size = !!tmp.sym.size, colour = !!tmp.sym.colour)) + 
      scale_x_discrete(limits = x.axis.names, breaks = x.axis.names) + 
      scale_y_discrete(limits = y.axis.names, breaks = y.axis.names) + 
      scale_size(name = colnames.to.cmp[1], range = nodes.size.range) + 
      scale_colour_gradientn(name = colnames.to.cmp[2], colours = nodes.colour.seq, values = nodes.colour.value.seq)
  gp.res <- gp.res + 
      labs(x = packed.infos[1, "from.cluster"], y = packed.infos[1, "to.cluster"]) + 
      theme_half_open(font_size = plot.font.size.base) + 
      background_grid() + 
      theme(axis.text.x = axis.text.x.pattern)
  # draw the final graph
  return(list(plot = gp.res, tables = packed.infos))
}
