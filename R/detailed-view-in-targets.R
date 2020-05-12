

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
#' @param sel.mode.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.mode}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.mode}.
#' @param sel.action.effect.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.action.effect}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.action.effect}.
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
  is.directional = TRUE,
  sel.mode.val = NULL,
  sel.action.effect.val = NULL
) {
  # predefined values extracted from database
  pred.mode <- kpred.mode
  pred.action.effect <- kpred.action.effect
  # as it plot either directed or undirected graphs, new definition of action effects are given as below
  # for "A---B",              given type: "undirected"  --- kaction.id.mapped[1]
  # for "A-->B" or "A<--B",   given type: "positive"  --- kaction.id.mapped[c(2,3)]
  # for "A--|B" or "A|--B",   given type: "negative"  --- kaction.id.mapped[c(4,5)]
  # for "A--oB" or "Ao--B",   given type: "unspecified" --- kaction.id.mapped[c(6,7)]
  #
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
  for (i in 1:length(list.interact.pairs)) {
    this.list <- list.interact.pairs[[i]]
    vertices.A.names <- append(vertices.A.names, this.list$act.A.genename)
    vertices.B.names <- append(vertices.B.names, this.list$act.B.genename)
  }
  vertices.A.names <- unique(vertices.A.names)
  vertices.B.names <- unique(vertices.B.names)
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
  vertices.A.pack.df <- left_join(vertices.A.pack.df, fgenes.part.A[, c("gene", "avg_logFC")], by = c("GeneName" = "gene"))
  # B# type -single
  vertices.B.apx.types <- anno.infos$type.B[, c("Gene.name", "Keyword.Name")]
  # B# loc
  vertices.B.pack.df <- left_join(vertices.B.pack.df, anno.infos$location.B[, c("Gene.name", "GO.Term.target")], by = c("GeneName" = "Gene.name"))
  #vertices.B.pack.df <- left_join(vertices.B.pack.df, anno.infos$type.B[, c("Gene.name", "Keyword.Name")], by = c("GeneName" = "Gene.name"))
  # !special rescue rule
  vertices.B.pack.df[which(is.na(vertices.B.pack.df[, "GO.Term.target"])), "GO.Term.target"] <- "Other"  # [rescue] 
  # B# logfc
  fgenes.part.B <- fgenes.remapped.all[which(fgenes.remapped.all$cluster == act.B.clustername), ]
  vertices.B.pack.df <- left_join(vertices.B.pack.df, fgenes.part.B[, c("gene", "avg_logFC")], by = c("GeneName" = "gene"))
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
  vertices.all.infos <- vertices.all.infos[, c("UID", "ClusterName", "GeneName", "Location", "LogFC")]  # rearrange the columns
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
  this.act.result <- list()  # list of data.frame
  for (i in 1:length(list.interact.pairs)) {
    this.list <- list.interact.pairs[[i]]
    act.A.genename <- this.list$act.A.genename
    act.A.UID <- intersect(which(vertices.all.infos[, "GeneName"] == act.A.genename), which(vertices.all.infos[, "ClusterName"] == afterV.A.clustername))
    act.B.genename <- this.list$act.B.genename
    act.B.UID <- intersect(which(vertices.all.infos[, "GeneName"] == act.B.genename), which(vertices.all.infos[, "ClusterName"] == afterV.B.clustername))
    act.infos <- this.list$action.infos
    if (nrow(act.infos) > 0) {
      for (j in 1:nrow(act.infos)) {
        this.row <- act.infos[j, ]
        rownames(this.row) <- NULL
        if (this.row["actionid"] == 1) {  # for undirected one, give two directed edge and special symbol representing those
          this.act.result <- append(this.act.result, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "undirected"))
          if (!is.directional) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "undirected"))
          }
        } else {
          if (this.row["actionid"] < 2 || this.row["actionid"] > 7) {
            stop(paste0("Undefined actionid from @param onepair.gmoc$actions.detailed[[", i, "]]!"))
          }
          if (this.row["actionid"] == 2) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "positive"))
          } 
          if (this.row["actionid"] == 3 && !is.directional) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "positive"))
          }
          if (this.row["actionid"] == 4) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "negative"))
          }
          if (this.row["actionid"] == 5 && !is.directional) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "negative"))
          }
          if (this.row["actionid"] == 6) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "unspecified"))
          }
          if (this.row["actionid"] == 7 && !is.directional) {
            this.act.result <- append(this.act.result, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "unspecified"))
          }
        }
      }
    }
  }
  edges.all.infos <- bind_rows(this.act.result)
  ### select target edges.part.infos and vertices.part.infos
  ## check if valid, sel.mode.val, sel.action.effect.val
  predefined.mode.list <- pred.mode
  predefined.action.effect.list <- pred.action.effect
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
#' \code{GenerateVEinfos()}.
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
      this.res <- ScatterSimple.Plot(this.points, this.ctp, 12, coords.xy.colnames = c("gx", "gy"))
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




#' Draw Circos plot for interaction pairs
#'
#' @description
#' This function draws circos plot for interaction pairs of several target clusters, and returns
#' the detailed result tables about these.
#'
#' @inheritParams Inside.DummyInteractPairsActed
#' @inheritParams Inside.DummyFgenes
#' @inheritParams Inside.DummyActionsRefDB
#' @param return.table.vals Logic. If TRUE, the table of interaction pairs will be given in return value, use \code{<-} to save the result.
#' @param select.interacts Character. Interacts names, e.g. \code{paste0("ClusterC", kClustersSplit, "ClusterD")}, and 
#' the direction will be from \code{"ClusterC"} to \code{"ClusterD"}.
#' @param clusters.select.auto Character. If parameter \code{select.interacts} is not specified, it gives all related clusters,
#' and all permutation of interacts between these clusters will be used.
#' @param is.directional Logic. It is passed to \code{GenerateVEinfos}. If TRUE, it only uses single direction 
#' of clusters' interaction, otherwise use the bi-directional.
#' @param sel.mode.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.mode}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.mode}.
#' @param sel.action.effect.val Character. If set NULL, it uses all values in global variables \code{CellTalkDB::kpred.action.effect}, or
#' please specify detailed and accurate values in subset of \code{CellTalkDB::kpred.action.effect}.
#' @param show.legend Logic. Show legend(TRUE) or not(FALSE).
#' @param cluster.label.cex Numeric. The size of letters when labeling cluster names.
#' @param cluster.colour Character. Predefined some colours in harmony, generated by colorbrew2, 
#' see \url{http://colorbrewer2.org} for details.
#' @param link.cor.colour Character. Colour of links, whose length can be 1 or the same as another parameter \code{pred.mode}.
#' @param link.cor.border.colour Character. Colour of border of links.
#' @param link.cor.arrow.drawn Logic. It decides if the arrows are drawn upon the links. If TRUE, draw the arrows, otherwise not.
#' @param link.cor.border.colour Character. Colour of link borders, whose value length can be 1 or the same as another parameter \code{pred.mode}.
#' @param link.cor.arrow.type Character. Arrow type of links, whose length must be same as another parameter \code{pred.action.effect}.
#' @param link.cor.arrow.length Numeric. Arrow length of links, whose length must be same as another parameter \code{pred.action.effect}.
#'
#' @details
#' This function draws circos plot to show multiple clusters invovled interaction pairs. The final circos plot
#' contains 2 tracks and 1 inner links gragh. (If you are not familiar with those terms, please see \url{} for details.)
#' \itemize{
#'   \item {base track: } {It represents all clusters that are drawn in different colours. 
#'         The colours are specified by parameter \code{cluster.colour}}
#'   \item {expression level track: } {It gives all genes' expression levels by values of "LogFC".
#'         The up-regulated ones are marked in yellowish, and the down-regulated ones are marked in green.}
#'   \item {links: } {It draws links pair-by-pair, and shows mode and action type of every interaction pair.}
#' }
#'
#'
#'
#' @import circlize
#' @importFrom graphics par plot.new
#' @importFrom gridBase gridOMI
#' @import grid
#' @importFrom ComplexHeatmap Legend packLegend
#'
#' @export
#'
PlotClusters.Circos <- function(
  interact.pairs.acted,
  fgenes.remapped.all,
  actions.ref.db,
  return.table.vals = FALSE,
  select.interacts = NULL,
  clusters.select.auto = NULL,
  is.directional = TRUE,
  sel.mode.val = NULL,
  sel.action.effect.val = NULL,
  show.legend = FALSE,
  cluster.label.cex = 0.6,
  cluster.colour = c("#fb8072", "#80b1d3", "#fdb462", "#ffffb3", "#bebada", "#b3de69", "#fccde5", "#8dd3c7", "#d9d9d9"),
  link.cor.colour = c("grey"),
  link.cor.border.colour = c("grey"),
  link.cor.arrow.drawn = c(TRUE, TRUE, TRUE, FALSE),
  link.cor.arrow.type = c("triangle", "ellipse", "curved", "circle"),
  link.cor.arrow.length = c(0.4, 0.3, 0.3, 0.2)
) {
  # pre-process set
  pred.mode <- kpred.mode
  pred.action.effect <- kpred.action.effect
  # pre-process check
  if (length(link.cor.colour) != 1 && length(link.cor.colour) < length(pred.mode)) {
    warning("Given insufficent mode-specific colour, use unified colour-'grey' instead.!")
    link.cor.colour <- c("grey")
  }
  # interact pairs selection and alignment of given parameter
  this.sel.interacts.names <- character()
  if (!is.null(select.interacts) && !(length(select.interacts) == 0)) {
    # check if these are valid interact.pairs' names
    logic.ifinnames <- select.interacts %in% interact.pairs.acted$name.allpairs
    if (sum(logic.ifinnames) == length(select.interacts)) {
      this.sel.interacts.names <- select.interacts
    } else {
      stop("Given interacts' names not found: ", paste0(select.interacts[which(logic.ifinnames == FALSE)], collapse = ", "), ".")
    }
  } else {
    if (!is.null(clusters.select.auto) && !((length(clusters.select.auto) == 0))) {
      # check if these are valid clusters' names
      logic.ifclusters <- clusters.select.auto %in% interact.pairs.acted$list.clusters
      if (sum(logic.ifclusters) == length(logic.ifclusters)) {
        # generating interacts' names
        for (i.cluster in clusters.select.auto) {
          for (j.cluster in clusters.select.auto) {
            if (i.cluster == j.cluster) {  # remove self-interactions
              next
            }
            this.sel.interacts.names <- append(this.sel.interacts.names, paste0(i.cluster, kClustersSplit, j.cluster))
          }
        }
      } else {
        stop("Given clusters' names not found: ", paste0(clusters.select.auto[which(logic.ifclusters == FALSE)], collapse = ", "), ".")
      }
    } else {  # use all interact pairs
      print("Given params are not available, automatically use all clusters to this process!")
      this.sel.interacts.names <- interact.pairs.acted$name.allpairs
    }
  }
  # get data of all interaction pairs as specified
  this.all.interacts <- list()  # get all pairs specified by interacts.names
  this.upd.interacts.names <- character()  # handle with the two same clusters interacts
  if (is.null(sel.mode.val)) {
    sel.mode.val <- pred.mode
  }
  if (is.null(sel.action.effect.val)) {
    sel.action.effect.val <- pred.action.effect
  }
  for (one.interact in this.sel.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.AB.1p <- ExtractTargetOnepairClusters(interact.pairs.acted, tmp.rcluster.A, tmp.rcluster.B)
    tmp.gmoc <- GenerateMapDetailOnepairClusters(tmp.AB.1p, actions.ref.db)
    tmp.VEinfos <- GenerateVEinfos(tmp.gmoc, fgenes.remapped.all,
                    is.directional = is.directional,
                    sel.mode.val = sel.mode.val, 
                    sel.action.effect.val = sel.action.effect.val)
    if (tmp.rcluster.A == tmp.rcluster.B) {
      tmp.rcluster.B <- paste0(tmp.rcluster.B, ".mirror")
    }
    tmp.new.interact.name <- paste0(tmp.rcluster.A, kClustersSplit, tmp.rcluster.B)
    this.upd.interacts.names <- append(this.upd.interacts.names, tmp.new.interact.name)
    this.all.interacts[[tmp.new.interact.name]] <- tmp.VEinfos
  }  # after GenerateVEinfos, .mirror names will appear
  # get unique clusters related
  this.related.clusters <- unique(unlist(strsplit(this.upd.interacts.names, split = kClustersSplit, fixed = TRUE)))
  ### extract all genes from related clusters
  # create data structure
  this.genes.clusters.based <- list()
  for (i.rcluster in this.related.clusters) {
    this.genes.clusters.based[[i.rcluster]] <- NULL  # data.frame
  }
  # data preparation for circos plot
  for (one.interact in this.upd.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.interact <- this.all.interacts[[one.interact]]
    tmp.vertex.infos <- tmp.interact$vertices.infos
    tmp.vinfos.A <- tmp.vertex.infos[which(tmp.vertex.infos[, "ClusterName"] == tmp.rcluster.A), ]
    this.genes.clusters.based[[tmp.rcluster.A]] <- rbind(this.genes.clusters.based[[tmp.rcluster.A]], tmp.vinfos.A[, c("GeneName", "Location", "LogFC")])
    tmp.vinfos.B <- tmp.vertex.infos[which(tmp.vertex.infos[, "ClusterName"] == tmp.rcluster.B), ]
    this.genes.clusters.based[[tmp.rcluster.B]] <- rbind(this.genes.clusters.based[[tmp.rcluster.B]], tmp.vinfos.B[, c("GeneName", "Location", "LogFC")])
  }
  # doing (1) unique, (2) sorting, (3) give x-axis values
  for (i.r.cluster in this.related.clusters) {
    this.genes.clusters.based[[i.r.cluster]] <- DoPartUnique(this.genes.clusters.based[[i.r.cluster]], c(1,2))
    # here, use sort in alphabet for "Location" to arrange these genes
    tmp.st <- this.genes.clusters.based[[i.r.cluster]]
    this.genes.clusters.based[[i.r.cluster]] <- tmp.st[order(tmp.st[, "Location"]), ]
    tmp.st <- this.genes.clusters.based[[i.r.cluster]]  # give x-axis values
    rownames(tmp.st) <- NULL  # this gives the default rownames, from 1 to nrow(*)
    tmp.st[, "x.val"] <- as.integer(rownames(tmp.st))
    this.genes.clusters.based[[i.r.cluster]] <- tmp.st
  }
  # get packed values
  this.genes.clusters.packed <- NULL  # data.frame
  for (i.r.cluster in this.related.clusters) {
    this.genes.clusters.packed <- rbind(this.genes.clusters.packed, cbind(this.genes.clusters.based[[i.r.cluster]], ClusterName = i.r.cluster))
  }
  #
  ### --- circos plot ---
  # grid the plot
  plot.new()
  kunit.plot = unit(1, "npc")
  if (show.legend) {
    pushViewport(viewport(x = unit(0, "npc"), y = unit(0.5, "npc"), 
      width = 0.8 * kunit.plot, height = kunit.plot, just = c("left", "center")))
  } else {
    pushViewport(viewport(x = unit(0, "npc"), y = unit(0.5, "npc"), 
      width = kunit.plot, height = kunit.plot, just = c("left", "center")))
  }
  # grid the region
  par(omi = gridOMI())
  par(new = TRUE)
  # extra clear to avoid unexpected program errors
  circos.clear()
  ## draw the plot
  circos.par("track.height" = 0.1, "gap.degree" = 3)
  circ.factors <- factor(this.genes.clusters.packed[, "ClusterName"])
  circ.xlim <- NULL
  for (i.rcluster in levels(circ.factors)) {
    tmp.this.xval <- this.genes.clusters.based[[i.rcluster]][, "x.val"]
    circ.xlim <- rbind(circ.xlim, c(min(tmp.this.xval) - 1, max(tmp.this.xval) + 1))
  }
  # initialize the circos
  circos.initialize(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], xlim = circ.xlim)
  # add base track
  kcolset.basetrack <- cluster.colour
  circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], ylim = c(0, 1),
    track.height = 0.05,
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), CELL_META$sector.index)
      circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2],
        col = kcolset.basetrack[(CELL_META$sector.numeric.index-1) %% length(kcolset.basetrack)+1])
  })
  # add gene name track
  circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], ylim = c(0, 32),
    track.height = 0.20, cell.padding = c(0, 0, 0, 0),
    bg.col = NA, bg.border = NA,
    panel.fun = function(x, y) {
      labels.genes <- this.genes.clusters.packed[which(this.genes.clusters.packed[, "ClusterName"] == CELL_META$sector.index), "GeneName"]
      circos.text(x, CELL_META$ycenter, labels.genes, 
        facing = "clockwise", niceFacing = TRUE,
        cex = cluster.label.cex)
  })
  # add or logfc
  kcolor.mval <- "#dcdcdc"
  if (min(this.genes.clusters.packed[, "LogFC"] > 0)) {
    circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], y = this.genes.clusters.packed[, "LogFC"],
      panel.fun = function(x, y) {
        circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = kcolor.mval)
        circos.rect(x - 0.5, CELL_META$cell.ylim[1], x + 0.5, y, col = "#efe91a")
    })
  } else {
    circos.track(factors = circ.factors, x = this.genes.clusters.packed[, "x.val"], y = this.genes.clusters.packed[, "LogFC"],
      panel.fun = function(x, y) {
        circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = kcolor.mval)
        circos.rect(x - 0.5, 0, x + 0.5, y, col = sapply(y, function(y) {if (y < 0) "green" else "#efe91a"}))
    })
  }
  ## draw links for interaction pairs
  thisf.used.mode <- NULL  # collect all used mode
  thisf.used.action.effect <- NULL  # collect all used action effects
  for (one.interact in this.upd.interacts.names) {
    tmp.related.clusters <- strsplit(one.interact, split = kClustersSplit, fixed = TRUE)[[1]]
    tmp.rcluster.A <- tmp.related.clusters[1]
    tmp.rcluster.B <- tmp.related.clusters[2]
    tmp.this.data <- this.all.interacts[[one.interact]]
    tmp.edges <- tmp.this.data$edges.infos
    tmp.vertices <- tmp.this.data$vertices.infos
    # get all links pairs with c(x1.val, x2.val) & (mode, action.effect)
    tmp.links.xvals <- apply(tmp.edges, MARGIN = 1, 
      vertices = tmp.vertices, gcb.ref = this.genes.clusters.based,
      FUN = function(x, vertices, gcb.ref) {
        # for x1
        uid.x1 <- as.integer(x["from"])  # get UID in this interaction pair
        ind.mv.x1 <- which(vertices[, "UID"] == uid.x1)  # get correspond row in vertices.infos
        tmp.cluster.x1 <- vertices[ind.mv.x1, "ClusterName"]
        tgcb.x1 <- gcb.ref[[tmp.cluster.x1]]  # get this cluster infos in gcb.ref
        ind.mgcb.x1 <- which(tgcb.x1[, "GeneName"] == vertices[ind.mv.x1, "GeneName"])  # find the same gene in gcb.ref[[this.cluster]]
        ind.mgcb.loc.x1 <- which(tgcb.x1[, "Location"] == vertices[ind.mv.x1, "Location"])
        val.x1 <- tgcb.x1[intersect(ind.mgcb.x1, ind.mgcb.loc.x1), "x.val"]  # get x.val in unique 1 row
        # for x2
        uid.x2 <- as.integer(x["to"])
        ind.mv.x2 <- which(vertices[, "UID"] == uid.x2)
        tmp.cluster.x2 <- vertices[ind.mv.x2, "ClusterName"]
        tgcb.x2 <- gcb.ref[[tmp.cluster.x2]]
        ind.mgcb.x2 <- which(tgcb.x2[, "GeneName"] == vertices[ind.mv.x2, "GeneName"])
        ind.mgcb.loc.x2 <- which(tgcb.x2[, "Location"] == vertices[ind.mv.x2, "Location"])
        val.x2 <- tgcb.x2[intersect(ind.mgcb.x2, ind.mgcb.loc.x2), "x.val"]  # get x.val in unique 1 row
        # other
        save.mode <- x["mode"]
        save.action.effect <- x["action.effect"]
        # return
        as.character(c(tmp.cluster.x1, val.x1, tmp.cluster.x2, val.x2, save.mode, save.action.effect))
    })
    tmp.links.xvals <- t(tmp.links.xvals)
    if (is.null(nrow(tmp.links.xvals)) || (nrow(tmp.links.xvals) == 0)) {  # go to next loop
      next
    }
    for (i in 1:nrow(tmp.links.xvals)) {
      # map different mode to different colour
      tmp.link.colour <- link.cor.colour  # assume using unified colour
      tmp.link.border <- link.cor.border.colour  # assume using unified colour
      if (length(link.cor.colour) > 1) {
        tmp.link.colour <- link.cor.colour[which(tmp.links.xvals[i, 5] == pred.mode)]
        tmp.link.border <- link.cor.border.colour[which(tmp.links.xvals[i, 5] == pred.mode)]
        thisf.used.mode <- c(thisf.used.mode, tmp.links.xvals[i, 5])
      }
      # arrow type
      thisf.used.action.effect <- c(thisf.used.action.effect, tmp.links.xvals[i, 6])
      ind.tmp.link.type <- which(tmp.links.xvals[i, 6] == pred.action.effect)
      tmp.this.arr.drawn <- link.cor.arrow.drawn[ind.tmp.link.type]
      tmp.this.arr.type <- link.cor.arrow.type[ind.tmp.link.type]
      tmp.this.arr.length <- link.cor.arrow.length[ind.tmp.link.type]
      if (tmp.this.arr.drawn == TRUE) {
        circos.link(tmp.links.xvals[i, 1], c(as.numeric(tmp.links.xvals[i, 2]) - 0.2, as.numeric(tmp.links.xvals[i, 2]) + 0.2), tmp.links.xvals[i, 3], as.numeric(tmp.links.xvals[i, 4]),
            col = tmp.link.colour, border = tmp.link.border, 
            directional = 1, arr.length = tmp.this.arr.length, arr.type = tmp.this.arr.type)
      } else {
        circos.link(tmp.links.xvals[i, 1], c(as.numeric(tmp.links.xvals[i, 2]) - 0.2, as.numeric(tmp.links.xvals[i, 2]) + 0.2), 
            tmp.links.xvals[i, 3], c(as.numeric(tmp.links.xvals[i, 4]) - 0.2, as.numeric(tmp.links.xvals[i, 4]) + 0.2),
            col = tmp.link.colour, border = tmp.link.border, 
            directional = 0, arr.length = tmp.this.arr.length, arr.type = tmp.this.arr.type)
      }
    }
  }
  circos.clear()
  # grid work, move to parent grid frame
  upViewport()
  ## draw the legend
  # legend exprs track
  leg.logfc.max <- max(this.genes.clusters.packed[, "LogFC"])
  leg.logfc.min <- min(this.genes.clusters.packed[, "LogFC"])
  exprs.p.at <- c("down-regulation", "up-regulation")
  exprs.p.col <- c("green", "#efe91a")
  if (leg.logfc.max * leg.logfc.min < 0) {
    # keep default values
  } else {
    if (leg.logfc.min > 0) {
      exprs.p.at <- c("up-regulation")
      exprs.p.col <- c("#efe91a")
    } else {
      exprs.p.at <- c("down-regulation")
      exprs.p.col <- c("green")
    }
  }
  lgd.exprs <- Legend(at = exprs.p.at, type = "points", 
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(col = exprs.p.col),
            title_position = "topleft", title = paste0("Track2", ":LogFC"))
  # legend arrow type
  thisf.used.action.effect <- levels(factor(thisf.used.action.effect))
  inds.lgd.arr.type <- match(thisf.used.action.effect, pred.action.effect)
  thisf.used.arr.type <- link.cor.arrow.type[inds.lgd.arr.type]
  lgd.links.arr.type <- Legend(at = c(toupper(thisf.used.arr.type), thisf.used.action.effect), ncol = 2,
            grid_width = unit(3, "mm"),
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(fill = c(rep("#ffffff", times = length(thisf.used.arr.type)), rep("#eeeeee", times = length(thisf.used.action.effect)))),
            title_position = "topleft", title = "Links-ActionType")
  # packed legend if no link colour is used,
  lgd.packed.list <- packLegend(lgd.exprs, lgd.links.arr.type)
  # conditionally, legend links colour
  if (length(link.cor.colour) > 1) {  # legend for link colour is required
    thisf.used.mode <- levels(factor(thisf.used.mode))
    inds.lgd.links.colour <- match(thisf.used.mode, pred.mode)
    thisf.used.links.colour <- link.cor.colour[inds.lgd.links.colour]
    lgd.links.colour <- Legend(at = thisf.used.mode, type = "lines",
            labels_gp = gpar(fontsize = 8),
            legend_gp = gpar(col = thisf.used.links.colour, lwd = 2),
            title_position = "topleft", title = paste0("Links-Mode"))
    lgd.packed.list <- packLegend(lgd.exprs, lgd.links.arr.type, lgd.links.colour)
  }
  pushViewport(viewport(x = 0.8 * kunit.plot, y = unit(0.5, "npc"), 
    width = 0.12 * kunit.plot, height = kunit.plot, just = c("left", "center")))
  if (show.legend) {
    grid.draw(lgd.packed.list)
  }
  # grid work
  upViewport()
  # end plotting

  # get result tables
  if (return.table.vals) {
    this.result.table.list <- list()
    for (i in 1:length(this.all.interacts)) {
      tmp.interact.name <- this.upd.interacts.names[i]
      tmp.it <- this.all.interacts[[i]]
      tmp.vertices <- tmp.it$vertices.infos
      tmp.edges <- tmp.it$edges.infos
      inds.tmp.from.match <- match(tmp.edges$from, tmp.vertices$UID)
      inds.tmp.to.match <- match(tmp.edges$to, tmp.vertices$UID)
      tmp.edges[, c("from.ClusterName", "from.GeneName")] <- tmp.vertices[inds.tmp.from.match, c("ClusterName", "GeneName")]
      tmp.edges[, c("to.ClusterName", "to.GeneName")] <- tmp.vertices[inds.tmp.to.match, c("ClusterName", "GeneName")]
      tmp.edges <- tmp.edges[, c("from", "from.ClusterName", "from.GeneName", "to", "to.ClusterName", "to.GeneName", "mode", "action.effect")]
      tmp.vertices <- tmp.vertices[, c("UID", "ClusterName", "GeneName", "Location", "LogFC")]
      tmp.list <- list(tmp.vertices, tmp.edges)
      names(tmp.list) <- c(paste0(tmp.interact.name, "-vertices"), paste0(tmp.interact.name, "-edges"))
      this.result.table.list <- c(this.result.table.list, tmp.list)
    }
    #end#
    return(list(table = this.result.table.list))
  } else {
    return(NULL)  # no return value
  }
}



