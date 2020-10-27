

#' Generate data about vertices and edges
#'
#' @description
#' This function uses detailed informations about one interaction pair(return value of 
#' \code{GenerateMapDetailOnepairClusters()}), to generate data for drawing relation plot.
#'
#' @param onepair.gmoc List. Return value of \code{\link{GenerateMapDetailOnepairClusters}}.
#' @inheritParams Inside.DummyFgenes 
#' @param is.directional Logic. If TRUE, it only generates VEinfos from \code{onepair.gmoc$clusters.name[1]} to \code{onepair.gmoc$clusters.name[2]}, 
#' otherwise bi-directional VEinfos will be generated.
#'
#' @details
#' This function only uses actions that are already known, i.e. recorded thoroughly in some databases,
#' to generate formatted data structure(with vertices and edges).
#'
#' In vertices, all gene informations are well recorded, and every gene is given one unique ID.
#'
#' In edges, it uses unique vertices IDs to contruct the linkes, and records mode and action.effect for every link.
#' 
#'
#' @return A list.
#' \itemize{
#'   \item {\code{cluster.name.A}&\code{cluster.name.B}:} {cluster names involved.}
#'   \item edges.infos: data.frame that records the edges(the interaction pairs).
#'   \item vertices.infos: data.frame that records the vertices(the genes).
#'   \item vertices.apx.type.A: data.frame that records the types(molecular functions) of A in gene pairs formatted as A-B.
#'   \item vertices.apx.type.B: data.frame that records the types(molecular functions) of B in gene pairs formatted as A-B.
#' }
#'
#'
#'
#' @importFrom dplyr left_join bind_rows
#'
#' @export
#'
GenerateVEinfos <- function(
  onepair.gmoc,
  fgenes.remapped.all,
  is.directional = TRUE
) {
  ### generate vertices list and edges list
  list.interact.pairs <- onepair.gmoc$actions.detailed
  anno.infos <- onepair.gmoc$anno.infos
  act.A.clustername <- onepair.gmoc$clusters.name[1]  #
  act.B.clustername <- onepair.gmoc$clusters.name[2]  # 
  if (length(list.interact.pairs) == 0) {  # if no actions.detailed exists, RETURN here
    stop(paste0("Given pair: ", act.A.clustername, "---", act.B.clustername, ", has no explicit actions defined in current settings!"))
  }

  ## --- vertices ---
  vertices.names <- character()
  vertices.A.names <- character()
  vertices.B.names <- character()
  vertices.names <- sapply(list.interact.pairs, function(x) {c(x$act.A.genename, x$act.B.genename)})
  vertices.names <- t(vertices.names)
  vertices.A.names <- unique(as.character(vertices.names[, 1]))
  vertices.B.names <- unique(as.character(vertices.names[, 2]))
  # pack vA vB to be df
  vertices.A.pack.df <- data.frame(GeneName = vertices.A.names, ClusterName = c(act.A.clustername), stringsAsFactors = FALSE)
  vertices.B.pack.df <- data.frame(GeneName = vertices.B.names, ClusterName = c(act.B.clustername), stringsAsFactors = FALSE)
  ## get other attributes about the vertices
  # A# type -single
  vertices.A.apx.types <- anno.infos$type.A[, c("Gene.name", "Keyword.Name")]
  # A# loc
  vertices.A.pack.df <- left_join(vertices.A.pack.df, anno.infos$location.A[, c("Gene.name", "GO.Term.target")], by = c("GeneName" = "Gene.name"))
  #vertices.A.pack.df <- left_join(vertices.A.pack.df, anno.infos$type.A[, c("Gene.name", "Keyword.Name")], by = c("GeneName" = "Gene.name"))
  # !special rescue rule
  vertices.A.pack.df[which(is.na(vertices.A.pack.df[, "GO.Term.target"])), "GO.Term.target"] <- "Other"  # [rescue]
  # A# logfc
  fgenes.part.A <- fgenes.remapped.all[which(fgenes.remapped.all$cluster == act.A.clustername), ]
  vertices.A.pack.df <- left_join(vertices.A.pack.df, fgenes.part.A, by = c("GeneName" = "gene"))
  # B# type -single
  vertices.B.apx.types <- anno.infos$type.B[, c("Gene.name", "Keyword.Name")]
  # B# loc
  vertices.B.pack.df <- left_join(vertices.B.pack.df, anno.infos$location.B[, c("Gene.name", "GO.Term.target")], by = c("GeneName" = "Gene.name"))
  #vertices.B.pack.df <- left_join(vertices.B.pack.df, anno.infos$type.B[, c("Gene.name", "Keyword.Name")], by = c("GeneName" = "Gene.name"))
  # !special rescue rule
  vertices.B.pack.df[which(is.na(vertices.B.pack.df[, "GO.Term.target"])), "GO.Term.target"] <- "Other"  # [rescue] 
  # B# logfc
  fgenes.part.B <- fgenes.remapped.all[which(fgenes.remapped.all$cluster == act.B.clustername), ]
  vertices.B.pack.df <- left_join(vertices.B.pack.df, fgenes.part.B, by = c("GeneName" = "gene"))
  # !! here, special rules will be applied upon if act.A.clustername == act.B.clustername
  afterV.A.clustername <- act.A.clustername
  afterV.B.clustername <- act.B.clustername
  if (act.A.clustername == act.B.clustername) {
    afterV.B.clustername <- paste0(act.B.clustername, ".mirror")  # [attention here!]
    vertices.B.pack.df$ClusterName <- afterV.B.clustername
  }
  vertices.all.infos <- rbind(vertices.A.pack.df, vertices.B.pack.df)
  vertices.all.infos$UID <- 1:nrow(vertices.all.infos)
  rownames(vertices.all.infos) <- NULL
  # change colnames in vertices.all
  tmp.cols.change <- match(c("GO.Term.target", "avg_logFC"), colnames(vertices.all.infos))
  colnames(vertices.all.infos)[tmp.cols.change] <- c("Location", "LogFC")
  tmp.cols.first5 <- c("UID", "ClusterName", "GeneName", "Location", "LogFC")
  vertices.all.infos <- vertices.all.infos[, c(tmp.cols.first5, setdiff(colnames(vertices.all.infos), tmp.cols.first5))]  # rearrange the columns
  # change colnames in apx.*
  colnames(vertices.A.apx.types) <- colnames(vertices.B.apx.types) <- c("GeneName", "Type")

  ## --- edges ---
  # predefined function
  gen.edges.vei.inside <- function(act.part1.UID, act.part2.UID, action.mode, action.effect) {
    # this function is to generate all permutation of act.part1.UID ~ act.part2.UID, e.g. A*2 B*3 will get 2*3 results
    tmp.all.pert <- lapply(act.part1.UID,
      act.part2.UID = act.part2.UID, action.mode = action.mode, action.effect = action.effect,
      FUN = function(x, act.part2.UID, action.mode, action.effect) {
        data.frame(from = x, to = act.part2.UID, action.mode = action.mode, action.effect = action.effect, stringsAsFactors = FALSE)
      }
    )
    tmp.all.pert  # return
  }
  # the process
  tmp.vertices.all.gene.inds <- tapply(1:nrow(vertices.all.infos), vertices.all.infos[, "GeneName"], function(x) {x})
  tmp.vertices.all.cluster.inds <- tapply(1:nrow(vertices.all.infos), vertices.all.infos[, "ClusterName"], function(x) {x})
  prog.bar.edge.collect <- progress::progress_bar$new(total = length(list.interact.pairs))
  prog.bar.edge.collect$tick(0)
  this.act.result <- lapply(list.interact.pairs, vertices.all.infos =  vertices.all.infos, 
    afterV.A.clustername = afterV.A.clustername, afterV.B.clustername = afterV.B.clustername, 
    tmp.gene.inds = tmp.vertices.all.gene.inds, tmp.cluster.inds = tmp.vertices.all.cluster.inds, 
    prog.bar.edge.collect = prog.bar.edge.collect, 
    gen.edges.vei.inside = gen.edges.vei.inside,  # function
    function(x, vertices.all.infos, afterV.A.clustername, afterV.B.clustername, 
      tmp.gene.inds, tmp.cluster.inds, prog.bar.edge.collect, gen.edges.vei.inside) {
      this.list <- x
      act.A.genename <- this.list$act.A.genename
      act.A.UID <- intersect(tmp.gene.inds[[act.A.genename]], tmp.cluster.inds[[afterV.A.clustername]])
      act.B.genename <- this.list$act.B.genename
      act.B.UID <- intersect(tmp.gene.inds[[act.B.genename]], tmp.cluster.inds[[afterV.B.clustername]])
      act.infos <- this.list$action.infos
      tmp.act.res <- list()  # for return
      if (nrow(act.infos) > 0) {
        for (j in 1:nrow(act.infos)) {
          this.row <- act.infos[j, ]
          rownames(this.row) <- NULL
          if (this.row["actionid"] == 1) {  # for undirected one, give two directed edge and special symbol representing those
            tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "undirected"))
            if (!is.directional) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "undirected"))
            }
          } else {
            if (this.row["actionid"] < 2 || this.row["actionid"] > 7) {
              stop(paste0("Undefined actionid from @param onepair.gmoc$actions.detailed[[", i, "]]!"))
            }
            if (this.row["actionid"] == 2) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "positive"))
            } 
            if (this.row["actionid"] == 3 && !is.directional) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "positive"))
            }
            if (this.row["actionid"] == 4) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "negative"))
            }
            if (this.row["actionid"] == 5 && !is.directional) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "negative"))
            }
            if (this.row["actionid"] == 6) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "unspecified"))
            }
            if (this.row["actionid"] == 7 && !is.directional) {
              tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "unspecified"))
            }
          }
        }
      }
      prog.bar.edge.collect$tick()
      bind_rows(tmp.act.res)
  })
  edges.all.infos <- bind_rows(this.act.result)
  #end# return
  VEinfos <- list(cluster.name.A = afterV.A.clustername, cluster.name.B = afterV.B.clustername,
    edges.infos = edges.all.infos, 
    vertices.infos = vertices.all.infos,
    vertices.apx.type.A = vertices.A.apx.types,
    vertices.apx.type.B = vertices.B.apx.types
    )
}



#' [TODO]
#'
#' @description
#' [TODO]
#'
#' @param VEinfos List. It contains informations about vertices and edges, and is exactly return value of
#' \code{GenerateVEinfos()}.
#' @param sel.mode.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.mode}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.mode}.
#' @param sel.action.effect.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.action.effect}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.action.effect}.
#'
#' @details
#' [TODO]
#'
#'
#' @export
#'
TrimVEinfos <- function(
  VEinfos, 
  sel.mode.val = NULL, 
  sel.action.effect.val = NULL
) {
  afterV.A.clustername <- VEinfos$cluster.name.A
  afterV.B.clustername <- VEinfos$cluster.name.B
  vertices.all.infos <- VEinfos$vertices.infos
  edges.all.infos <- VEinfos$edges.infos
  vertices.A.apx.types <- VEinfos$vertices.A.apx.types
  vertices.B.apx.types <- VEinfos$vertices.B.apx.types
  # as it plot either directed or undirected graphs, new definition of action effects are given as below
  # for "A---B",              given type: "undirected"  --- kaction.id.mapped[1]
  # for "A-->B" or "A<--B",   given type: "positive"  --- kaction.id.mapped[c(2,3)]
  # for "A--|B" or "A|--B",   given type: "negative"  --- kaction.id.mapped[c(4,5)]
  # for "A--oB" or "Ao--B",   given type: "unspecified" --- kaction.id.mapped[c(6,7)]
  #
  ### select target edges.part.infos and vertices.part.infos
  ## check if valid, sel.mode.val, sel.action.effect.val
  predefined.mode.list <- kpred.mode
  predefined.action.effect.list <- kpred.action.effect
  if ((sum(sel.mode.val %in% predefined.mode.list) == length(sel.mode.val) ||
     is.null(sel.mode.val)) &&
    (sum(sel.action.effect.val %in% predefined.action.effect.list) == length(sel.action.effect.val) ||
     is.null(sel.action.effect.val))) {
    # --- mode ---
    if (is.null(sel.mode.val)) {
      edges.part.infos <- edges.all.infos
    } else {
      inds.part.select <- match(edges.all.infos[, "mode"], sel.mode.val)
      edges.part.infos <- edges.all.infos[which(!is.na(inds.part.select)), ]
    }
    # --- action.effect ---
    if (!is.null(sel.action.effect.val)) {
      inds.part.select.ex <- match(edges.part.infos[, "action.effect"], sel.action.effect.val)
      edges.part.infos <- edges.part.infos[which(!is.na(inds.part.select.ex)), ]
    }
    # recheck if nrow() > 0
    if (nrow(edges.part.infos) == 0) {
      stop("No given subset of interactions between cluster: ", act.A.clustername, " and cluster: ", act.B.clustername, "!")  # afterV.B.clustername is not used for displaying warnings
    }
    part.select.vertices <- unique(c(levels(factor(edges.part.infos[, "from"])), levels(factor(edges.part.infos[, "to"]))))
    vertices.part.infos <- vertices.all.infos[match(part.select.vertices, vertices.all.infos[, "UID"]), ]
    # remapping UIDs
    vertices.part.infos$UID <- 1:nrow(vertices.part.infos)  # as target vertices may be less than total vertices, so remapping the UID
    inds.part.new.id.from <- match(edges.part.infos[, "from"], rownames(vertices.part.infos))
    inds.part.new.id.to   <- match(edges.part.infos[, "to"], rownames(vertices.part.infos))
    edges.part.infos[, "from"] <- vertices.part.infos$UID[inds.part.new.id.from]
    edges.part.infos[, "to"] <- vertices.part.infos$UID[inds.part.new.id.to]
    rownames(vertices.part.infos) <- NULL  # make rownames be equal to UID
    # set the apx* vars
    inds.part.A.vx <- which(vertices.A.apx.types[, "GeneName"] %in% vertices.part.infos[which(vertices.part.infos$ClusterName == afterV.A.clustername), "GeneName"])
    vertices.A.apx.types <- vertices.A.apx.types[inds.part.A.vx, ]
    inds.part.B.vx <- which(vertices.B.apx.types[, "GeneName"] %in% vertices.part.infos[which(vertices.part.infos$ClusterName == afterV.B.clustername), "GeneName"])
    vertices.B.apx.types <- vertices.B.apx.types[inds.part.B.vx, ]
  } else {
    not.inlist.mode <- sel.mode.val[which(!(sel.mode.val %in% predefined.mode.list))]
    not.inlist.action.effect <- sel.action.effect.val[which(!(sel.action.effect.val %in% predefined.action.effect.list))]
    stop(paste0("Error in given @param, with mode not in list: ", paste0(not.inlist.mode, collapse = ", "), 
      ", with action.effect not in list: ", paste0(not.inlist.action.effect, collapse = ", "),
      ", please recheck these given above!"))
  }
  #end# return
  VEinfos <- list(cluster.name.A = afterV.A.clustername, cluster.name.B = afterV.B.clustername,
    edges.infos = edges.part.infos, 
    vertices.infos = vertices.part.infos,
    vertices.apx.type.A = vertices.A.apx.types,
    vertices.apx.type.B = vertices.B.apx.types
    )
}





# [inside usage]
# This function is a random points scatter algorithm.
# The return value of this function is the @param data.veinfo plus 2 new columns
# that given by @param coords.xy.colnames
#
ScatterSimple.Plot <- function(
  data.veinfo,
  center.xy,  # (x, y)
  ext.len,  # circle, if non-circle is wanted, use coords.trans functions
  sample.gap.degree = 30,  # default will be 30
  density.half.near = 1 / 3.0,  # default will be 1/3
  coords.xy.colnames = c("gx", "gy")
) {
  # pre-check
  if (nrow(data.veinfo) < 1) {
    if (!is.null(colnames(data.veinfo))) {
      tmp.ret <- data.frame(ttx = numeric(0), tty = numeric(0))
      colnames(tmp.ret) <- coords.xy.colnames
      return(cbind(data.veinfo, tmp.ret))
    } else {
      return(NULL)
    }
  }
  #
  if (length(coords.xy.colnames) != 2) {
    stop("Unexpected colnames whose length <!=> 2!")
  }
  #
  this.puts.cnt <- ext.len - 2  # exclude the edges and the center
  if (this.puts.cnt < 1) {
    stop("Too small ext.len")
  }
  # sample
  if (sample.gap.degree < 1 || sample.gap.degree > 360) {
    sample.gap.degree <- 30  # the default
  }
  this.deg.splits <- floor(360 / sample.gap.degree)
  # check capacity
  if (this.puts.cnt * this.deg.splits <= nrow(data.veinfo)) {
    stop("Cannot be located inside!")
  }

  # scatter preparation
  radius.near.center <- ceiling(density.half.near * this.puts.cnt)
  rad.start.away.center <- radius.near.center + 1
  if (rad.start.away.center >= ext.len) {
    stop(paste0("Too small ext.len: ", ext.len))
  }
  tmp.sel.near.c <- floor(nrow(data.veinfo) * density.half.near)
  data.sel.near.c <- 1:ifelse(tmp.sel.near.c >= 1, tmp.sel.near.c, 1)
  data.in.near.c <- data.veinfo[data.sel.near.c, ]
  data.in.away.c <- data.veinfo[setdiff(1:nrow(data.veinfo), data.sel.near.c), ]
  ## scatter capacity
  # in near center
  sp.rad.near.c <- sample(1:radius.near.center, length(1:radius.near.center))
  sp.deg.near.c <- sample(1:this.deg.splits, length(1:this.deg.splits))
  sp.total.near.c <- length(sp.rad.near.c) * length(sp.deg.near.c)
  if (sp.total.near.c < nrow(data.in.near.c)) {
    stop("Capacity error in near center!")
  }
  # in away from center
  sp.rad.away.c <- sample(rad.start.away.center:this.puts.cnt, length(rad.start.away.center:this.puts.cnt))
  sp.deg.away.c <- sample(1:this.deg.splits, length(1:this.deg.splits))
  sp.total.away.c <- length(sp.rad.away.c) * length(sp.deg.away.c)
  if (sp.total.away.c < nrow(data.in.away.c)) {
    stop("Capacity error away from center!")
  }
  ### scatter process (use deg.* as ref and extend it)
  ## get coords
  gen.coords.ssp <- function(data.input, sp.deg.seq, sp.rad.seq, sample.gap.degree, center.xy) {
    tmp.suit.deg.times <- ceiling(nrow(data.input) / length(sp.deg.seq))
    tmp.suit.rad.times <- ceiling(tmp.suit.deg.times * length(sp.deg.seq) / length(sp.rad.seq))
    tmp.sref.deg <- rep(sp.deg.seq, times = tmp.suit.deg.times) * sample.gap.degree
    tmp.sref.rad <- rep(sp.rad.seq, times = tmp.suit.rad.times)
    length(tmp.sref.deg) <- length(tmp.sref.rad) <- nrow(data.input)
    coords.res <- data.frame(x.s = cos(tmp.sref.deg) * tmp.sref.rad + center.xy[1],
                 y.s = sin(tmp.sref.deg) * tmp.sref.rad + center.xy[2])
    coords.res
  }
  # in near center
  coords.near.c <- gen.coords.ssp(data.in.near.c, sp.deg.near.c, sp.rad.near.c, sample.gap.degree, center.xy)
  # in away from center
  coords.away.c <- gen.coords.ssp(data.in.away.c, sp.deg.away.c, sp.rad.away.c, sample.gap.degree, center.xy)
  ## merge back to df
  data.in.near.c[, coords.xy.colnames] <- coords.near.c
  data.in.away.c[, coords.xy.colnames] <- coords.away.c
  data.ret <- rbind(data.in.near.c, data.in.away.c)
  # return
  data.ret
}





# [inside usage]
# enlarge and rotate
# enlarge will be directly put in x-axis and y-axis, so limit its usage.
# 
Inside.TransCoords.Enlarge.Rotate <- function(
  orig.coords,
  enlarge.xy.times,  # c(x.times, y.times)
  rotate.degree,
  rotate.xy.ref = c(0, 0)
) {
  this.xvals <- orig.coords[, 1]
  this.yvals <- orig.coords[, 2]
  # enlarge
  this.xvals <- this.xvals * enlarge.xy.times[1]
  this.yvals <- this.yvals * enlarge.xy.times[2]
  # rotate
  rotate.x.ref <- rotate.xy.ref[1]
  rotate.y.ref <- rotate.xy.ref[2]
  deg.trans <- (pi * rotate.degree / 180)
  val.cos <- cos(deg.trans)
  val.sin <- sin(deg.trans)
  target.xvals <- val.cos * this.xvals - val.sin * this.yvals + (1 - val.cos) * rotate.x.ref + val.sin * rotate.y.ref
  target.yvals <- val.sin * this.xvals + val.cos * this.yvals - val.sin * rotate.x.ref + (1 - val.cos) * rotate.y.ref
  # return prep
  this.ret <- data.frame(x.r = target.xvals, y.r = target.yvals, stringsAsFactors = FALSE)
  colnames(this.ret) <- colnames(orig.coords)
  # return
  this.ret
}





#' Draw CellPlot for one pair of clusters
#' 
#' @description
#' This function analyzes the interaction pairs between two target clusters, and gives one vivid two-cell
#' graph and several result tables in the return values.
#'
#' @param VEinfos List. It contains informations about vertices and edges, and is exactly return value of
#' \code{GenerateVEinfos()} or \code{TrimVEinfos()}.
#' @param area.extend.times Numeric. If a warning given like "Cannot be located inside!" or something else, 
#' one should change this paramemter to be larger to get all vertices allocated.
#' @param hide.locations.A Character. It applies extra limitation on the locations of A in gene pairs formatted as A-B.
#' @param hide.types.A Character. It applies extra limitation on the types(molecular functions) of A in gene pairs formatted as A-B.
#' @param hide.locations.B Character. It applies extra limitation on the locations of B in gene pairs formatted as A-B.
#' @param hide.types.B Character. It applies extra limitation on the types(molecular functions) of B in gene pairs formatted as A-B.
#' @param hide.sole.vertices Character. It hides sole vertices which have no available edges.
#' @param nodes.size Numeric. Size of nodes.
#' @param nodes.colour Character. Colour of nodes.
#' @param nodes.alpha Numeric. Alpha of nodes.
#' @param nodes.shape Vector. Use shape reprentable IDs(in integer) or directly shape description(in character). 
#' See \pkg{ggplot2} for details
#' @param nodes.stroke Numeric. Stroke size of nodes. See details in \pkg{ggplot2}.
#' @param nodes.label.size Numeric. Size of label on nodes.
#' @param nodes.label.colour Character. Colour of label on nodes.
#' @param link.size Numeric. Size of link width.
#' @param link.colour Character. Colour of links, the length should be same as \code{CellTalkDB::kpred.action.effect}.
#' @param link.alpha Numeric. Alpha of link.
#' @param link.arrow.angle Numeric. Angle of link arrow.
#' @param link.arrow.length Numeric. Length of link arrow.
#' @param link.arrow.type Character. Type of link arrow, either \code{open} or \code{closed}.
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#'
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#'
#' @export
#'
GetResult.PlotOnepairClusters.CellPlot <- function(
  VEinfos,
  area.extend.times = 1,  # [TODO]
  hide.locations.A = NULL,
  hide.types.A = NULL,
  hide.locations.B = NULL,
  hide.types.B = NULL,
  hide.sole.vertices = TRUE,  # if TRUE, remove those edges cannot formed vertices
  nodes.size = 3,
  nodes.colour = "grey",
  nodes.alpha = 0.8,
  nodes.shape = 21,
  nodes.stroke = 0.7,
  nodes.label.size = 2,
  nodes.label.colour = "white",
  link.size = 0.1,
  link.colour = c("red", "green", "blue", "lightgrey"),
  link.alpha = 0.7,
  link.arrow.angle = 20,
  link.arrow.length = 10,
  link.arrow.type = "open"
) {
  # precheck
  if (length(link.colour) != length(kpred.action.effect)) {
    length(link.colour) <- length(kpred.action.effect)
    link.colour[which(is.na(link.colour))] <- link.colour[which(!is.na(link.colour))][1]  # use 1st one to all other
    warning("Given link colour are shorter than expected, and are automatically extended.")
  }
  #
  act.A.clustername <- VEinfos$cluster.name.A
  act.B.clustername <- VEinfos$cluster.name.B
  edges.infos <- VEinfos$edges.infos
  vertices.infos <- VEinfos$vertices.infos  # infos included: location
  vertices.apx.type.A <- VEinfos$vertices.apx.type.A
  vertices.apx.type.B <- VEinfos$vertices.apx.type.A

  # hide some vertices
  # A
  tmp.is.A <- which(vertices.infos$ClusterName == act.A.clustername)
  # A loc
  if (is.null(hide.locations.A)) {
    tmp.inds.A.loc <- integer(0)
  } else {
    tmp.inds.A.loc <- intersect(tmp.is.A, which(vertices.infos$Location %in% hide.locations.A))
  }
  # A type
  if (is.null(hide.types.A)) {
    tmp.inds.A.type <- integer(0)
  } else {
    tmp.inds.A.type <- intersect(tmp.is.A, which(vertices.infos$GeneName %in% (vertices.apx.type.A[which(vertices.apx.type.A$Type %in% hide.types.A), "GeneName"]) ))
  }
  # B
  tmp.is.B <- which(vertices.infos$ClusterName == act.B.clustername)
  # B loc
  if (is.null(hide.locations.B)) {
    tmp.inds.B.loc <- integer(0)
  } else {
    tmp.inds.B.loc <- intersect(tmp.is.B, which(vertices.infos$Location %in% hide.locations.B))
  }
  # B type
  if (is.null(hide.types.B)) {
    tmp.inds.B.type <- integer(0)
  } else {
    tmp.inds.B.type <- intersect(tmp.is.B, which(vertices.infos$GeneName %in% (vertices.apx.type.B[which(vertices.apx.type.B$Type %in% hide.types.B), "GeneName"]) ))
  }
  # merge hide inds
  tmp.inds.hide <- c(tmp.inds.A.loc, tmp.inds.A.type, tmp.inds.B.loc, tmp.inds.B.type)
  if (length(tmp.inds.hide) != 0) {  # [TODO] remove from edges
    # so then the vertices
    vertices.infos <- vertices.infos[-tmp.inds.hide, ]
    # hide sole vertices
    if (hide.sole.vertices == TRUE) {
      tmp.m.from <- match(edges.infos[, "from"], rownames(vertices.infos))  # match to old UID
      tmp.m.to <- match(edges.infos[, "to"], rownames(vertices.infos))
      inds.m.edges <- intersect(which(!is.na(tmp.m.from)), which(!is.na(tmp.m.to)))
      tmp.UID.keep <- unique(c(edges.infos[inds.m.edges, "from"], edges.infos[inds.m.edges, "to"]))
      vertices.infos <- vertices.infos[which(vertices.infos[, "UID"] %in% tmp.UID.keep), ]
    }
    # check if all edges were removed
    tmp.is.A <- which(vertices.infos$ClusterName == act.A.clustername)
    tmp.is.B <- which(vertices.infos$ClusterName == act.B.clustername)
    if (length(tmp.is.A) == 0 || length(tmp.is.B) == 0) {
      stop("Hide too much vertices to be NO available edges at all.")
    }
    vertices.infos$UID <- 1:nrow(vertices.infos)
    inds.part.new.id.from <- match(edges.infos[, "from"], rownames(vertices.infos))
    inds.part.new.id.to   <- match(edges.infos[, "to"], rownames(vertices.infos))
    edges.infos[, "from"] <- vertices.infos$UID[inds.part.new.id.from]
    edges.infos[, "to"] <- vertices.infos$UID[inds.part.new.id.to]
    # remove NAs
    edges.infos <- edges.infos[intersect(which(!is.na(edges.infos[, "from"])), which(!is.na(edges.infos[, "to"]))), ]
    rownames(vertices.infos) <- NULL  # make rownames be equal to UID
    # set the apx* vars
    tmp.inds.A.apx.types <- which(vertices.apx.type.A[, "GeneName"] %in% vertices.infos[which(vertices.infos$ClusterName == act.A.clustername), "GeneName"])
    vertices.apx.type.A <- vertices.apx.type.A[tmp.inds.A.apx.types, ]
    tmp.inds.B.apx.types <- which(vertices.apx.type.B[, "GeneName"] %in% vertices.infos[which(vertices.infos$ClusterName == act.B.clustername), "GeneName"])
    vertices.apx.type.B <- vertices.apx.type.B[tmp.inds.B.apx.types, ]
  }



  # special locations when ploting
  pred.loc.special <- c("Plasma Membrane", "Extracellular Region", "Other")
  # others
  pred.loc.common.A <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.A.clustername)])), pred.loc.special)
  pred.loc.common.B <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.B.clustername)])), pred.loc.special)

  ## A
  tmp.inds.is.A <- which(vertices.infos$ClusterName == act.A.clustername)
  # locations
  this.inds.A.pmem <- intersect(which(vertices.infos$Location == "Plasma Membrane"), tmp.inds.is.A)
  this.inds.A.exm <- intersect(which(vertices.infos$Location == "Extracellular Region"), tmp.inds.is.A)
  this.inds.A.other <- intersect(which(vertices.infos$Location == "Other"), tmp.inds.is.A)
  this.inds.A.common <- lapply(pred.loc.common.A, 
    inds.is.A = tmp.inds.is.A, vertices.infos = vertices.infos,
    function(x, inds.is.A, vertices.infos) {
      intersect(which(vertices.infos$Location == x), inds.is.A)
    })
  names(this.inds.A.common) <- pred.loc.common.A
  ## B
  tmp.inds.is.B <- which(vertices.infos$ClusterName == act.B.clustername)
  # locations
  this.inds.B.pmem <- intersect(which(vertices.infos$Location == "Plasma Membrane"), tmp.inds.is.B)
  this.inds.B.exm <- intersect(which(vertices.infos$Location == "Extracellular Region"), tmp.inds.is.B)
  this.inds.B.other <- intersect(which(vertices.infos$Location == "Other"), tmp.inds.is.B)
  this.inds.B.common <- lapply(pred.loc.common.B, 
    inds.is.B = tmp.inds.is.B, vertices.infos = vertices.infos,
    function(x, inds.is.B, vertices.infos) {
      intersect(which(vertices.infos$Location == x), inds.is.B)
    })
  names(this.inds.B.common) <- pred.loc.common.B


  # Simple base plot
  this.base.graph <- ggplot()
  # shared compts
  # comp A in -x axis
  this.plot.A.cell <- data.frame(x = c(-175, -175, -55, -55) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.A.pmem <- data.frame(x = c(-55, -55, -45, -45) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.A.exm <- data.frame(x = c(-45, -45, -15, -15) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.A.other <- data.frame(x = c(-195, -195, -175, -175) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)

  # comp B in +x axis
  this.plot.B.cell <- data.frame(x = c(55, 55, 175, 175) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.pmem <- data.frame(x = c(45, 45, 55, 55) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.exm <- data.frame(x = c(15, 15, 45, 45) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.other <- data.frame(x = c(175, 175, 195, 195) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)

  ## plot the base
  this.graph.raw.base <- this.base.graph +
      geom_polygon(data = this.plot.A.other, aes(x, y), colour = "grey") +
      geom_polygon(data = this.plot.A.cell, aes(x, y)) +
      geom_polygon(data = this.plot.A.pmem, aes(x, y), colour = "orange") +
      geom_polygon(data = this.plot.A.exm, aes(x, y)) +
      geom_polygon(data = this.plot.B.other, aes(x, y), colour = "grey") +
      geom_polygon(data = this.plot.B.cell, aes(x, y)) +
      geom_polygon(data = this.plot.B.pmem, aes(x, y), colour = "orange") +
      geom_polygon(data = this.plot.B.exm, aes(x, y))
  # label the cluster
  this.label.clusters <- data.frame(
    x.lc = c(-105, 105) * area.extend.times,
    y.lc = c(52, 52) * area.extend.times,
    r.label.c = c(act.A.clustername, act.B.clustername),
    r.label.nudge.y.c = c(0, 0) * area.extend.times)
  this.graph.raw.base <- this.graph.raw.base + 
      geom_text(data = this.label.clusters, 
        mapping = aes(x = x.lc, y = y.lc, label = r.label.c),
        colour = "black", size = 4, 
        vjust = -0.9, nudge_y = this.label.clusters$r.label.nudge.y.c)
  # label every big region
  this.label.base.items <- data.frame(
    x.l = c(-115, 115, -50, 50, -30, 30, -185, 185) * area.extend.times,
    y.l = rep(c(45), times = 8) * area.extend.times,
    r.label = rep(c("Cytoplasm", "Plasma Membrane", "Extracellular", "Other"), each = 2),
    r.label.nudge.y = rep(c(0, 2, 0, 0), each = 2) * area.extend.times)
  this.graph.raw.base <- this.graph.raw.base + 
      geom_text(data = this.label.base.items, 
        mapping = aes(x = x.l, y = y.l, label = r.label),
        colour = "black", size = 3, 
        vjust = -1, nudge_y = this.label.base.items$r.label.nudge.y)

  ### plot points
  this.graph.add.ps <- this.graph.raw.base
  this.vx.ext.infos <- data.frame()
  ## plot tools
  # NULL
  ## plot Plasma Membrane [TODO] other way
  # A
  this.pmem.A.ctp.xy <- c(-50, 0) * area.extend.times
  this.pmem.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.pmem, ], this.pmem.A.ctp.xy, 5 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.pmem.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.pmem.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 9), rotate.degree = 0, rotate.xy.ref = this.pmem.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.A.vx.ext)
  # B
  this.pmem.B.ctp.xy <- c(50, 0) * area.extend.times
  this.pmem.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.pmem, ], this.pmem.B.ctp.xy, 5 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.pmem.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.pmem.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 9), rotate.degree = 0, rotate.xy.ref = this.pmem.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.B.vx.ext)

  ## plot Extracellular Region
  # A
  this.exm.A.ctp.xy <- c(-30, 0) * area.extend.times
  this.exm.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.exm, ], this.exm.A.ctp.xy, 10 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.exm.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, rotate.xy.ref = this.exm.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.A.vx.ext)
  # B
  this.exm.B.ctp.xy <- c(30, 0) * area.extend.times
  this.exm.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.exm, ], this.exm.B.ctp.xy, 10 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.exm.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, rotate.xy.ref = this.exm.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.B.vx.ext)

  ## plot Other
  # A
  this.other.A.ctp.xy <- c(-185, 0) * area.extend.times
  this.other.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.other, ], this.other.A.ctp.xy, 10 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.other.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, rotate.xy.ref = this.other.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.A.vx.ext)
  # B
  this.other.B.ctp.xy <- c(185, 0) * area.extend.times
  this.other.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.other, ], this.other.B.ctp.xy, 10 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
  this.other.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, rotate.xy.ref = this.other.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.B.vx.ext)

  ## plot common
  gen.coords.put.ps <- function(x, inds.ref, ctp.ref, vertices.infos) {
    this.res <- NULL
    this.inds <- inds.ref[[x]]
    if (length(this.inds) > 0) {
      tmp.ind.m <- match(x, ctp.ref$Map.Items)
      this.ctp <- as.numeric(ctp.ref[tmp.ind.m, c("tpx", "tpy")])
      this.points <- vertices.infos[this.inds, ]
      this.res <- ScatterSimple.Plot(this.points, this.ctp, 12 * area.extend.times, coords.xy.colnames = c("gx", "gy"))
    }
    this.res
  }

  # A common
  this.common.A.ctp.xy <- data.frame(tpx = rep(c(-70, -100, -130,  -160), each = 3) * area.extend.times,
                     tpy = rep(c(30, 0, -30), times = 4) * area.extend.times)
  # tmp-ly shortcuts  [TODO]
  tmp.plc <- pred.loc.common.A
  length(tmp.plc) <- 12
  this.common.A.ctp.xy$Map.Items <- tmp.plc
  #
  this.common.A.vx.inlist <- lapply(names(this.inds.A.common),
    inds.ref = this.inds.A.common,
    ctp.ref = this.common.A.ctp.xy,
    vertices.infos = vertices.infos,
    FUN = gen.coords.put.ps
    )
  this.vx.ext.infos <- rbind(this.vx.ext.infos, bind_rows(this.common.A.vx.inlist))

  # B common
  this.common.B.ctp.xy <- data.frame(tpx = rep(c(70, 100, 130, 160), each = 3) * area.extend.times,
                     tpy = rep(c(30, 0, -30), times = 4) * area.extend.times)
  # tmp-ly shortcuts  [TODO]
  tmp.plc <- pred.loc.common.B
  length(tmp.plc) <- 12
  this.common.B.ctp.xy$Map.Items <- tmp.plc
  #
  this.common.B.vx.inlist <- lapply(names(this.inds.B.common),
    inds.ref = this.inds.B.common,
    ctp.ref = this.common.B.ctp.xy,
    vertices.infos = vertices.infos,
    FUN = gen.coords.put.ps
    )
  this.vx.ext.infos <- rbind(this.vx.ext.infos, bind_rows(this.common.B.vx.inlist))


  ## plot dash-lined common region
  # use this.common.A.ctp.xy, this.common.B.ctp.xy
  tmp.common.ctp.xy <- rbind(this.common.A.ctp.xy[which(!is.na(this.common.A.ctp.xy$Map.Items)), ],
                 this.common.B.ctp.xy[which(!is.na(this.common.B.ctp.xy$Map.Items)), ])
  if (nrow(tmp.common.ctp.xy) > 0) {
    tmp.common.polygons <- list()
    const.extlen <- 15 * area.extend.times
    for (i in 1:nrow(tmp.common.ctp.xy)) {
      this.cx <- tmp.common.ctp.xy[i, "tpx"]
      this.cy <- tmp.common.ctp.xy[i, "tpy"]
      this.ploygon.xy <- data.frame(x = rep(c(this.cx - 14.5 * area.extend.times, this.cx + 14.5 * area.extend.times), each = 2), 
        y = c(this.cy - 14.5 * area.extend.times, this.cy + 14.5 * area.extend.times, this.cy + 14.5 * area.extend.times, this.cy - 14.5 * area.extend.times))
      tmp.common.polygons <- append(tmp.common.polygons, 
        geom_polygon(data = this.ploygon.xy, aes(x, y),
          colour = "red", linetype = "dashed"
        ))
    }
    this.graph.add.ps <- this.graph.add.ps + tmp.common.polygons

    tmp.common.ctp.xy[, "label.pos.tpy"] <- tmp.common.ctp.xy[, "tpy"] + 13.5 * area.extend.times
    # use short names

    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Endoplasmic Reticulum"), "Map.Items"] <- "ER"
    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Golgi Apparatus"), "Map.Items"] <- "Golgi"
    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Cytoplasm"), "Map.Items"] <- "OtherCytoplasm"
    #
    this.graph.add.ps <- this.graph.add.ps + 
      geom_text(data = tmp.common.ctp.xy, 
        mapping = aes(x = tpx, y = label.pos.tpy, label = Map.Items),
        colour = "red", size = 2)
  }
  

  # -------------
  this.graph.gplot <- this.graph.add.ps
  ## plot vertices
  this.graph.gplot <- this.graph.gplot +
    geom_point(data = this.vx.ext.infos,
      mapping = aes(x = gx, y = gy),
      size = nodes.size,
      colour = nodes.colour,
      alpha = nodes.alpha,
      shape = nodes.shape,
      stroke = nodes.stroke)
  this.graph.gplot <- this.graph.gplot + 
    geom_text(data = this.vx.ext.infos,
      mapping = aes(x = gx, y = gy, label = GeneName),
      size = nodes.label.size,
      colour = nodes.label.colour)

  # -------------

  ## plot edges
  inds.e.from.match <- match(edges.infos$from, this.vx.ext.infos$UID)
  inds.e.to.match <- match(edges.infos$to, this.vx.ext.infos$UID)
  edges.infos[, c("from.gx", "from.gy")] <- this.vx.ext.infos[inds.e.from.match, c("gx", "gy")]
  edges.infos[, c("to.gx", "to.gy")] <- this.vx.ext.infos[inds.e.to.match, c("gx", "gy")]
  #
  pred.mode <- kpred.mode
  pred.action.effect <- kpred.action.effect
  show.mode.val <- levels(factor(edges.infos$mode))
  show.action.effect.val <- levels(factor(edges.infos$action.effect))

  # draw edges (remove colour edges)
  names(link.colour) <- pred.action.effect  # for aesthetics

  draw.add.comps.list <- list()
  # use reverse order to plot positive & negative above the unspecified or undirected
  for (j in rev(1:length(pred.action.effect))) {
    this.edges.infos <- edges.infos[which(edges.infos[, "action.effect"] == pred.action.effect[j]), ]
    for (i in 1:length(pred.mode)) {
      if (!(pred.action.effect[j] %in% show.action.effect.val)) {
        break
      }
      inds.this.mode <- which(this.edges.infos[, "mode"] == pred.mode[i])
      draw.add.comps.list <- append(draw.add.comps.list, 
        geom_segment(mapping = aes(x = from.gx, y = from.gy, xend = to.gx, yend = to.gy, colour = action.effect),
          data = this.edges.infos[inds.this.mode, ],
          alpha = link.alpha, size = link.size,
          arrow = arrow(angle = link.arrow.angle, length = unit(link.arrow.length, "pt"), type = link.arrow.type)
          )
      )
    }
  }
  draw.add.comps.list <- append(draw.add.comps.list, 
    scale_colour_manual(name = "Action.effect", values = link.colour, breaks = names(link.colour), aesthetics = "colour")
  )
  this.graph.gplot <- this.graph.gplot + draw.add.comps.list

  # other settings
  this.graph.gplot <- this.graph.gplot + 
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme(panel.background = element_blank())

  ## res
  list(plot = this.graph.gplot, table = list(vertices.infos = this.vx.ext.infos, edges.infos = edges.infos))
}





#' [TODO]
#'
#' @description
#' [TODO]
#'
#' @details
#' [TODO]
#'
#'
#'
#' @export
#'
Tool.formula.onLogFC.default <- function(
  data.f,  # vector
  data.b   # vector
) {
  if (length(data.f) != length(data.b)) {
    stop("Unexpected non-identical length data input!")
  }
  return(data.f + data.b)
}



#' [TODO]
#'
#' @description
#' [TODO]
#'
#' @details
#' [TODO]
#'
#'
#'
#' @export
#'
Tool.formula.onPValAdj.default <- function(
  data.f,  # vector
  data.b   # vector
) {
  if (length(data.f) != length(data.b)) {
    stop("Unexpected non-identical length data input!")
  }
  tmp.f <- abs(log10(data.f))
  tmp.b <- abs(log10(data.b))
  max.f <- max(tmp.f[which(is.finite(tmp.f))])
  max.b <- max(tmp.b[which(is.finite(tmp.b))])
  tmp.f[which(is.infinite(tmp.f))] <- 10 * max.f
  tmp.b[which(is.infinite(tmp.b))] <- 10 * max.b
  return(log10(tmp.f * tmp.b + 1))
}





#' [TODO]
#'
#' @description
#' [TODO]
#'
#' @details
#' [TODO]
#'
#' @return
#' [TODO]
#'
#'
#'
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
GetResult.PlotOnepairClusters.GeneCrosstalk <- function(
  VEinfos,
  fgenes.remapped.all,
  direction.A.to.B = TRUE,
  colnames.to.cmp = c("LogFC", "p_val_adj"),
  range.to.use = list("LogFC" = c(-Inf, Inf), "p_val_adj" = c(-Inf, Inf)),
  formula.to.use = list(Tool.formula.onLogFC.default, Tool.formula.onPValAdj.default),
  axis.order.xy = c("AlphaBet", "AlphaBet"),  # how to order axis in final plot. Can also be one of colnames.to.cmp
  axis.order.xy.decreasing = c(TRUE, TRUE),  # order direction
  nodes.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
  nodes.colour.value.seq = c(0.0, 0.5, 1.0),
  nodes.size.range = c(2, 8)
) {
  # pre-check
  if (length(colnames.to.cmp) < length(formula.to.use)) {
    warning("Formulae given too many, the last ones will not be used!")
  }
  if (length(colnames.to.cmp) > length(formula.to.use)) {
    stop("Items to be compared are not given sufficient formulae!")
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
  edges.fullinfos <- unique(edges.fullinfos)
  # restrict to only one direction
  if (!is.null(direction.A.to.B)) {
    if (direction.A.to.B == TRUE) {
      tmp.inds <- intersect(which(edges.fullinfos[, "from.cluster"] == act.A.clustername), which(edges.fullinfos[, "to.cluster"] == act.B.clustername))
      edges.fullinfos <- edges.fullinfos[tmp.inds, ]
    } else {
      tmp.inds <- intersect(which(edges.fullinfos[, "to.cluster"] == act.A.clustername), which(edges.fullinfos[, "from.cluster"] == act.B.clustername))
      edges.fullinfos <- edges.fullinfos[tmp.inds, ]
    }
  }  # else use all rows in edges.infos
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
    tmp.res <- formula.to.use[[i]](packed.infos[, tmp.sel.cols[1]], packed.infos[, tmp.sel.cols[2]])
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
      theme_half_open(font_size = 12) + 
      background_grid() + 
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  # draw the final graph
  return(list(plot = gp.res, tables = packed.infos))
}
