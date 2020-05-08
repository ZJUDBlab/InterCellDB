

#' Extract one pair of interaction
#'
#' @description
#' `ExtractTargetOnepairClusters` returns one pair of interaction from user-given cluster names.
#'
#' @param interact.pairs.acted A list. The return value of \code{\link{AnalyzeClustersInteracts}}.
#' @param cluster.name.C Character. Name of one cluster.
#' @param cluster.name.D Character. Name of one cluster.
#'
#' @details
#' The direction of this function is C-D, corresponding to coordinates X-Y in plot, and 
#' the subset options applied before will be inherited, which means that options
#' applied on X-axis clusters will be also on C, so as options on Y-axis will be on D.
#'
#' @return A list.
#' \itemize{
#'   \item clusters.name: names of clusters.
#'   \item bt.pairs: a list of all interaction pairs between specified 2 clusters (\code{given in clusters.name}).
#'   \item anno.infos: a list of lists. The sublists are 
#'         \itemize{
#'           \item location.A: it records locations of A in gene pairs formatted as A-B
#'           \item location.B: it records locations of B in gene pairs formatted as A-B
#'           \item type.A: it records molecular functions of A in gene pairs formatted as A-B
#'           \item type.B: it records molecular functions of B in gene pairs formatted as A-B
#'         }
#' }
#'
#'
#'
#' @export
#'
ExtractTargetOnepairClusters <- function(
  interact.pairs.acted,
  cluster.name.C,
  cluster.name.D
) {
  if (!(cluster.name.C %in% interact.pairs.acted$list.clusters &&
      cluster.name.D %in% interact.pairs.acted$list.clusters)) {  # check if given cluster.name.* are in clusters.name list
    stop("Error: Target one-pair clusters not found, as giving undefined name of clusters.")
  }
  this.clusters.name <- c(cluster.name.C, cluster.name.D)
  this.sel.name <- paste0(cluster.name.C, kClustersSplit, cluster.name.D) 
  ## slim - locations
  this.loc.slim.tg.cols <- c("Gene.name", "GO.Term.target", "score")
  this.loc.save.cols <- c("GeneID", "Gene.name", "GO.Term.target")
  # loc.A
  tmp.A.raw.loc <- interact.pairs.acted$anno.allpairs$location.A[[this.sel.name]]
  this.A.loc <- DoPartUnique(tmp.A.raw.loc, cols.select = match(this.loc.slim.tg.cols, colnames(tmp.A.raw.loc)))
  this.A.loc <- this.A.loc[, this.loc.save.cols]
  # loc.B
  tmp.B.raw.loc <- interact.pairs.acted$anno.allpairs$location.B[[this.sel.name]]
  this.B.loc <- DoPartUnique(tmp.B.raw.loc, cols.select = match(this.loc.slim.tg.cols, colnames(tmp.B.raw.loc)))
  this.B.loc <- this.B.loc[, this.loc.save.cols]
  ## slim - types (already done in the database level)
  #
  #end# res
  list(clusters.name = this.clusters.name, 
       bt.pairs = interact.pairs.acted$data.allpairs[[this.sel.name]],
       anno.infos = list(location.A = this.A.loc, 
                         location.B = this.B.loc, 
                         type.A = interact.pairs.acted$anno.allpairs$type.A[[this.sel.name]],
                         type.B = interact.pairs.acted$anno.allpairs$type.B[[this.sel.name]])
      )
}





# [inside usage]
# This function is to interpret meanings hidden in action.<taxonomy>.ref.db,
# and uses columns(mode, is_directional, a_is_acting) in the database.
# This function will only be used inside GenerateMapDetailOnepairClusters()
#
Inside.CollectActionMapping <- function(
  onerow.info,
  act.A.genename,
  act.B.genename
) {
  res.mode <- onerow.info["mode"]
  res.actionid <- 1
  ifisdir <- if (onerow.info["is_directional"] == 't') TRUE else FALSE
  if (!ifisdir) {
    if (onerow.info["a_is_acting"] == 't') {
      stop("Exception: database is broken by mannually modification, as is_directional is f but a_is_acting is t!")
    } # else res.actionid <- 1 (stay default)
  } else {
    if (onerow.info["a_is_acting"] == 'f') {
      # no change, res.actionid <- 1 (stay default)
    } else {
      # is_directional TRUE, a_is_acting TRUE, confirm that it is A-act-upon-B
      ifconv <- if (onerow.info["inter.GeneName.A"] == act.A.genename) 0 else 1
      if (onerow.info["action"] == "") {
        res.actionid <- match("A--oB", kaction.id.mapped) + ifconv
      } else {
        if (onerow.info["action"] == "activation") {  # the database record positive actions in "activation"
          res.actionid <- match("A-->B", kaction.id.mapped) + ifconv
        } else {
          if (onerow.info["action"] == "inhibition") {  # the database record negative actions in "inhibition"
            res.actionid <- match("A--|B", kaction.id.mapped) + ifconv
          } else {
            stop("Exception: database is broken by mannually modification, as undefined values appear in column('action')!")
          }
        }
      }
    }
  }
  #end# return
  c(as.character(res.mode), as.character(res.actionid))
}





#' Generate interaction pairs with their actions
#' 
#' @description
#' This function uses actions.ref.db to distinguish as well as collect pairs whose direction of interact actions
#' are known in some degree and are recorded, and pairs with no detailed informations.
#'
#' @param clusters.onepair.select  A list. Return value of \code{\link{ExtractTargetOnepairClusters}}.
#' @inheritParams Inside.DummyActionsRefDB
#'
#' @return A list.
#' \itemize{
#'   \item clusters.name: names of clusters involved. Its length is of 2.
#'   \item anno.infos: a list of lists. The sublists are 
#'         \itemize{
#'           \item location.A: it records locations of A in gene pairs formatted as A-B
#'           \item location.B: it records locations of B in gene pairs formatted as A-B
#'           \item type.A: it records molecular functions of A in gene pairs formatted as A-B
#'           \item type.B: it records molecular functions of B in gene pairs formatted as A-B
#'         }
#'   \item actions.detailed: A list of detailed information about interaction pairs whose actions are recorded in actions.ref.db.
#' }
#'
#'
#'
#' @import progress
#'
#' @export
#'
GenerateMapDetailOnepairClusters <- function(
  clusters.onepair.select,
  actions.ref.db
) {
  bt.pairs <- clusters.onepair.select$bt.pairs
  if (nrow(bt.pairs) == 0) {
    stop("Error: No pairs available in given parameter.",
      "Generating from ", paste0(clusters.onepair.select$clusters.name, collapse = " & "), " is failed."
    )
  }
  bt.pairs.result <- list(
    clusters.name = clusters.onepair.select$clusters.name,
    anno.infos = clusters.onepair.select$anno.infos,
    actions.detailed = list()  # pairs that have known actions in actions.ref.db
  )
  print(paste0("Generating from ", paste0(clusters.onepair.select$clusters.name, collapse = " and "), "."))
  this.act.detailed.id <- 0
  prog.bar.gmoc <- progress::progress_bar$new(total = nrow(bt.pairs))
  prog.bar.gmoc$tick(0)
  for (i in 1:nrow(bt.pairs)) {
    onerow.tmp <- bt.pairs[i, ]
    inds.conv.A <- which(actions.ref.db$inter.GeneName.A == onerow.tmp$inter.GeneName.A)
    inds.conv.B <- which(actions.ref.db$inter.GeneName.B == onerow.tmp$inter.GeneName.B)
    inds.conv.res <- intersect(inds.conv.A, inds.conv.B)
    inds.rev.A <- which(actions.ref.db$inter.GeneName.B == onerow.tmp$inter.GeneName.A)
    inds.rev.B <- which(actions.ref.db$inter.GeneName.A == onerow.tmp$inter.GeneName.B)
    inds.rev.res <- intersect(inds.rev.A, inds.rev.B)
    onerow.pairs <- actions.ref.db[union(inds.conv.res, inds.rev.res), ]  # all actions records for this pair
    # extracting
    action.A.genename <- onerow.tmp$inter.GeneName.A
    action.B.genename <- onerow.tmp$inter.GeneName.B
    action.A.exprs <- onerow.tmp$inter.Exprs.A
    action.B.exprs <- onerow.tmp$inter.Exprs.B
    action.A.logfc <- onerow.tmp$inter.LogFC.A
    action.B.logfc <- onerow.tmp$inter.LogFC.B
    if (nrow(onerow.pairs) > 0) {  # pairs that have known actions in actions.ref.db
      action.infos <- apply(onerow.pairs, MARGIN = 1, 
        act.A.genename = action.A.genename, act.B.genename = action.B.genename,
        function(x, act.A.genename, act.B.genename) {
          Inside.CollectActionMapping(x, act.A.genename, act.B.genename)
        }
      )
      action.infos <- t(action.infos)
      ## data trimming
      # deliminate those 1s when higher IDs exist, which means if directional is determined, leave out old ambiguous ones.
      action.infos.df <- data.frame(
        mode = as.character(action.infos[, 1]), 
        actionid = as.integer(action.infos[, 2]),
        stringsAsFactors = FALSE
      )
      inds.del.1 <- which(action.infos.df$actionid == 1)  # to examine if more specific mode has been recorded
      inds.rest <- which(action.infos.df$actionid != 1)
      logic.del.1 <- action.infos.df[inds.del.1, "mode"] %in% action.infos.df[inds.rest, "mode"]
      action.infos.df.exm <- rbind(action.infos.df[inds.del.1[which(logic.del.1 == FALSE)], ],
                     action.infos.df[inds.rest, ]
                   )
      action.infos.df.exm <- unique(action.infos.df.exm)
      onerow.res <- list(
        act.A.genename = action.A.genename,
        act.B.genename = action.B.genename,
        act.A.exprs = action.A.exprs,
        act.B.exprs = action.B.exprs,
        act.A.logfc = action.A.logfc,
        act.B.logfc = action.B.logfc,
        action.infos = action.infos.df.exm
      )   
    } else {  # pairs that dont's have recorded known actions
      onerow.res <- list(
        act.A.genename = action.A.genename,
        act.B.genename = action.B.genename,
        act.A.exprs = action.A.exprs,
        act.B.exprs = action.B.exprs,
        act.A.logfc = action.A.logfc,
        act.B.logfc = action.B.logfc,
        action.infos = data.frame(mode = "other", actionid = 1, stringsAsFactors = FALSE)
      )
    }
    this.act.detailed.id <- this.act.detailed.id + 1
    bt.pairs.result$actions.detailed[[this.act.detailed.id]] <- onerow.res
    prog.bar.gmoc$tick()
  }
  #end# return
  bt.pairs.result
}





# [inside usage]
# This function is to collect hierachical information inside the raw data.
#
# @param onepair.gmoc List. Return value of GenerateMapDetailOnepairClusters().
#
Inside.CollectHierachyOnepairClusters <- function(
  onepair.gmoc
) {
  if (length(onepair.gmoc$actions.detailed) <= 0) {
    stop("No specific interaction pairs with detailed action effects recorded. Please recheck input!")
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
  op.clustername <- onepair.gmoc$clusters.name[1]
  rv.clustername <- onepair.gmoc$clusters.name[2]
  group.proto.exprs.variations <- list(
    clusters.name = list(cluster.C = op.clustername, cluster.D = rv.clustername),
    Cup.Dup = group.proto.action.effects,  # group.proto.action.effects, and the same for all below g.C*.D*
    Cup.Ddn = group.proto.action.effects, 
    Cdn.Dup = group.proto.action.effects,
    Cdn.Ddn = group.proto.action.effects
  )
  # process
  group.act.res <- group.proto.exprs.variations
  for (i.list in 1:length(onepair.gmoc$actions.detailed)) {
    this.pair.info <- onepair.gmoc$actions.detailed[[i.list]]
    this.pair.actions <- this.pair.info$action.infos
    # other identity infos
    this.A.genename <- this.pair.info$act.A.genename
    this.B.genename <- this.pair.info$act.B.genename
    this.A.logfc  <- this.pair.info$act.A.logfc
    this.B.logfc  <- this.pair.info$act.B.logfc
    # do grouping upon gene.A & gene.B expression changes(logfc)
    target.AB.in.group <- NULL
    if (this.A.logfc > 0 && this.B.logfc > 0) {
      target.AB.in.group <- "Cup.Dup"
    }
    if (this.A.logfc > 0 && this.B.logfc < 0) {
      target.AB.in.group <- "Cup.Ddn"
    }
    if (this.A.logfc < 0 && this.B.logfc > 0) {
      target.AB.in.group <- "Cdn.Dup"
    }
    if (this.A.logfc < 0 && this.B.logfc < 0) {
      target.AB.in.group <- "Cdn.Ddn"
    }
    # do grouping upon action effects
    this.pair.act.ids <- levels(factor(this.pair.actions[, "actionid"]))
    target.AB.actions <- character()  # set length = 0
    # A---B
    if ((tmp.actionid = 1) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A---B")
    }
    # A-->B
    if ((tmp.actionid = 2) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A-->B")
    }
    # A<--B
    if ((tmp.actionid = 3) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A<--B")
    }
    # A--|B
    if ((tmp.actionid = 4) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A--|B")
    }
    # A|--B
    if ((tmp.actionid = 5) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A|--B")
    }
    # A--oB
    if ((tmp.actionid = 6) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "A--oB")
    }
    # Ao--B
    if ((tmp.actionid = 7) %in% this.pair.act.ids) {
      target.AB.actions <- append(target.AB.actions, "Ao--B")
    }
    for (i.type in target.AB.actions) {
      group.act.res[[target.AB.in.group]][[i.type]] <- rbind(group.act.res[[target.AB.in.group]][[i.type]], 
        data.frame(act.C.genename = this.A.genename,
               act.D.genename = this.B.genename,
               act.C.logfc = this.A.logfc,
               act.D.logfc = this.B.logfc,
               stringsAsFactors = FALSE
      ))
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
#
Inside.PlotShowGrouping <- function(
  onepair.dgsa,
  show.exprs.change = c("Cup.Dup", "Cup.Ddn", "Cdn.Dup", "Cdn.Ddn"),
  show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"),
  scale.fill.discrete,
  legend.key.size,
  legend.title,
  legend.text,
  legend.box.margin,
  caption.label.size,
  caption.label.y,
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
    guides(fill = guide_legend(title = "Action.Effects")) + 
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
  caption.label <- paste0("---symbols---\nA,B: genes, C,D: clusters.\nC: ", 
    onepair.dgsa$clusters.name[["cluster.C"]], ", D: ", onepair.dgsa$clusters.name[["cluster.D"]])
  caption.plot <- ggdraw() + 
    draw_label(caption.label, size = caption.label.size, y = caption.label.y, vjust = caption.label.vjust) + 
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
#' @param onepair.gmoc List. Return value of \code{\link{GenerateMapDetailOnepairClusters}}.
#' @param show.exprs.change Character. Use exprssion level change of clusters to select part of data to be shown.
#' @param show.action.effects Character. Select some action effects to be put in the result.
#' @param plot.axis.x.name Character. X-axis name when plot is shown.
#' @param plot.axis.y.name Character. Y-axis name when plot is shown.
#'
#'
#' @return A list. Use \code{Tool.ShowGraph()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @examples
#' # select only Cup.Dup
#' \dontrun{
#' GetResult.SummaryOnepairClusters(
#'   onepair.gmoc,
#'   show.exprs.change = c("Cup.Dup")
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
  onepair.gmoc,
  show.exprs.change = c("Cup.Dup", "Cup.Ddn", "Cdn.Dup", "Cdn.Ddn"),
  show.action.effects = c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"),
  scale.fill.discrete = scale_fill_brewer(palette = "Set3"),
  legend.key.size = unit(6, "mm"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 6),
  caption.label.size = 12,
  caption.label.y = 1,
  caption.label.vjust = 1
) {
  # process
  onepair.dgsa <- Inside.CollectHierachyOnepairClusters(onepair.gmoc)
  # pre-result check
  notin.exprs.c <- show.exprs.change[-which(show.exprs.change %in% c("Cup.Dup", "Cup.Ddn", "Cdn.Dup", "Cdn.Ddn"))]
  if (length(notin.exprs.c) > 0) {
    warning("Given @param show.exprs.change has undefined values: ", 
      paste0(notin.exprs.c, collapse = ", "), "\n and will be ignored!")
  }
  show.exprs.change <- setdiff(show.exprs.change, notin.exprs.c)
  if (length(show.exprs.change) == 0) {
    stop("Select available data less than 1!")
  }
  # pre-check 2
  notin.act.effect <- show.action.effects[-which(show.action.effects %in% c("A-->B", "A<--B", "A--|B", "A|--B", "A--oB", "Ao--B", "A---B"))]
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
    caption.label.y = caption.label.y,
    caption.label.vjust = caption.label.vjust
    )
  ## result table
  res.table.list <- list(
    Cup.Dup = data.frame(),
    Cup.Ddn = data.frame(),
    Cdn.Dup = data.frame(),
    Cdn.Ddn = data.frame()
  )
  this.cluster.C <- onepair.dgsa$clusters.name$cluster.C
  this.cluster.D <- onepair.dgsa$clusters.name$cluster.D
  for (i.exc in show.exprs.change) {
    this.exc <- onepair.dgsa[[i.exc]]
    for (i.type in show.action.effects) {
      if (nrow(this.exc[[i.type]]) == 0) {
        next
      }
      res.table.list[[i.exc]] <- rbind(res.table.list[[i.exc]],
        data.frame(cluster.C = this.cluster.C,
               gene.A = this.exc[[i.type]]$act.C.genename,
               logfc.A = this.exc[[i.type]]$act.C.logfc,
               cluster.D = this.cluster.D,
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
