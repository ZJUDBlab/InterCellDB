


# [inside usage]
# This function is a random points scatter algorithm.
# The return value of this function is the @param data.veinfo plus 2 new columns
# that given by @param coords.xy.colnames
#
ScatterSimple.Plot <- function(
  data.veinfo,
  outside.cut.percent,  # proportion to area.extend.times to get some around area clear of points
  center.xy,  # (x, y)
  ext.len,  # circle, if non-circle is wanted, use coords.trans functions
  radius.gap.factor = 1,
  sample.shift.degree = 0,  # default goes from +x axis to measure degree
  sample.gap.degree = 60,  # default will be 30
  density.half.near = 1 / 4.0,  # default will be 1/4
  coords.xy.colnames = c("gx", "gy")
) {
  # warning or stop for area.extend.times
  ssp.words <- " Please increases the area.extend.times by multiply it by 10."
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
  this.puts.cnt <- ext.len * (1 - outside.cut.percent)  # 0.1 percent mostly able to exclude the edges and the center
  if (this.puts.cnt < 1) {
    stop("Too small ext.len!", ssp.words)
  }
  # sample
  if (sample.gap.degree < 1 || sample.gap.degree > 360) {
    sample.gap.degree <- 30  # the default
  }
  this.deg.splits <- floor(360 / sample.gap.degree)
  # check capacity
  if (this.puts.cnt * this.deg.splits <= nrow(data.veinfo)) {
    stop("Cannot be located inside!", ssp.words)
  }

  # scatter preparation
  radius.near.center <- ceiling(0.5 * this.puts.cnt)  # half radius
  rad.start.away.center <- radius.near.center + 1
  if (rad.start.away.center >= ext.len) {
    stop(paste0("Too small ext.len: ", ext.len, ssp.words))
  }
  tmp.sel.near.c <- floor(nrow(data.veinfo) * density.half.near)
  data.sel.near.c <- seq_len(ifelse(tmp.sel.near.c >= 1, tmp.sel.near.c, 1))
  data.in.near.c <- data.veinfo[data.sel.near.c, ]
  data.in.away.c <- data.veinfo[setdiff(seq_len(nrow(data.veinfo)), data.sel.near.c), ]
  ## scatter capacity
  # in near center
  tmp.used.near.center <- floor(radius.near.center / radius.gap.factor)
  sp.rad.near.c <- seq_len(tmp.used.near.center) * radius.gap.factor
  sp.deg.near.c <- seq_len(this.deg.splits)
  #sp.rad.near.c <- sample(1:tmp.used.near.center, length(1:tmp.used.near.center)) * radius.gap.factor
  #sp.deg.near.c <- sample(1:this.deg.splits, length(1:this.deg.splits))
  sp.total.near.c <- length(sp.rad.near.c) * length(sp.deg.near.c)
  if (sp.total.near.c < nrow(data.in.near.c) || tmp.used.near.center == 0) {
    stop("Capacity error in near center!", ssp.words)
  }
  # in away from center
  tmp.used.away.center <- floor((this.puts.cnt - rad.start.away.center + 1) / radius.gap.factor)
  sp.rad.away.c <- seq_len(tmp.used.away.center) * radius.gap.factor + rad.start.away.center
  sp.deg.away.c <- seq_len(this.deg.splits)
  #sp.rad.away.c <- sample(1:tmp.used.away.center, length(1:tmp.used.away.center)) * radius.gap.factor + rad.start.away.center
  #sp.deg.away.c <- sample(1:this.deg.splits, length(1:this.deg.splits))
  sp.total.away.c <- length(sp.rad.away.c) * length(sp.deg.away.c)
  if (sp.total.away.c < nrow(data.in.away.c) || tmp.used.away.center == 0) {
    stop("Capacity error away from center!", ssp.words)
  }
  ### scatter process (use deg.* as ref and extend it)
  ## get coords
  gen.coords.ssp <- function(data.input, sp.deg.seq, sp.rad.seq, sample.gap.degree, sample.shift.degree, center.xy) {
    # full sequence of deg and rad
    full.sref.deg <- (rep(sp.deg.seq, times = length(sp.rad.seq)) * sample.gap.degree + sample.shift.degree) * (pi / 180)
    full.sref.rad <- rep(sp.rad.seq, each = length(sp.deg.seq))
    # needed count
    tmp.scnt <- sample(seq_along(full.sref.deg), nrow(data.input))
    tmp.sref.deg <- full.sref.deg[tmp.scnt]
    tmp.sref.rad <- full.sref.rad[tmp.scnt]
    coords.res <- data.frame(x.s = cos(tmp.sref.deg) * tmp.sref.rad + center.xy[1],
                 y.s = sin(tmp.sref.deg) * tmp.sref.rad + center.xy[2])
    coords.res
  }
  # in near center
  coords.near.c <- gen.coords.ssp(data.in.near.c, sp.deg.near.c, sp.rad.near.c, sample.gap.degree, sample.shift.degree, center.xy)
  # in away from center
  coords.away.c <- gen.coords.ssp(data.in.away.c, sp.deg.away.c, sp.rad.away.c, sample.gap.degree, sample.shift.degree, center.xy)
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
  enlarge.xy.ref = c(0, 0),
  rotate.xy.ref = c(0, 0)
) {
  this.xvals <- orig.coords[, 1]
  this.yvals <- orig.coords[, 2]
  # enlarge
  this.xvals <- (this.xvals - enlarge.xy.ref[1]) * enlarge.xy.times[1] + enlarge.xy.ref[1]
  this.yvals <- (this.yvals - enlarge.xy.ref[2]) * enlarge.xy.times[2] + enlarge.xy.ref[2]
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





# This function is to generate coordinates for unit size circle.
Inside.GetCoords.ZeroCircle <- function(
  center.x = 0, 
  center.y = 0, 
  radius = 1.0
) {
  # generating coordinates counter-clockwise
  whole.degree <- 360
  gradient.degree <- 3  # now 120 splits, and in the future, probably use smaller value
  n.cnt <- floor(whole.degree / gradient.degree)
  uc.x.part <- uc.y.part <- numeric(n.cnt)
  for (i in seq_len(n.cnt)) {
    this.degree <- gradient.degree * (i - 1)  # start from 0 degree
    uc.x.part[i] <- cos((this.degree / 180) * pi) * radius + center.x
    uc.y.part[i] <- sin((this.degree / 180) * pi) * radius + center.y
  }
  uc.all.coords <- data.frame(x.cc = uc.x.part, y.cc = uc.y.part)
  #end# return 
  uc.all.coords
}





#' Draw CellPlot for one pair of clusters
#' 
#' @description
#' This function analyzes the interaction pairs between two target clusters, and gives one vivid two-cell
#' graph and several result tables in the return values.
#'
#' @inheritParams Inside.DummyVEinfos
#' @param area.extend.times Numeric. Its default value is 10 which can handle most cases. If a warning given like "Cannot be located inside!" or something else, 
#' one should change this paramemter to be larger to get all vertices allocated.
#' @param hide.locations.A Character. It applies extra limitation on the locations of A in gene pairs formatted as A-B.
#' @param hide.types.A Character. It applies extra limitation on the types(molecular functions) of A in gene pairs formatted as A-B.
#' @param hide.locations.B Character. It applies extra limitation on the locations of B in gene pairs formatted as A-B.
#' @param hide.types.B Character. It applies extra limitation on the types(molecular functions) of B in gene pairs formatted as A-B.
#' @param hide.sole.vertices Character. It hides sole vertices which have no available edges.
#' @param expand.gap.radius.list Numeric. It defines the minimum distance between points(genes) for each plotting area.
#' @param expand.shift.degree.list Numeric. It defines the begin degree that points(genes) are to be drawn. The degree is calculated counter-clockwise.
#' @param expand.gap.degree.list Numeric. It defines the way that points(genes) arrange. If it is set 180 and shift degree is set 90, then points will be aligned in vertical line. 
#' If it is set 90 and shift degree is set 0, then points will be put counter-clockwise from horizontal line to vertical and then back to horizontal and finally at vertical place.
#' @param expand.center.density.list Numeric. It defines the density in each plotting area. Higher value means points concentrating more around 
#' the center of the area, and lower value means more sparse.
#' @param nodes.size.range [TODO]
#' @param nodes.fill [TODO]
#' @param nodes.colour Character. Colour of nodes.
#' @param nodes.alpha Numeric. Alpha of nodes.
#' @param nodes.shape Vector. Use shape reprentable IDs(in integer) or directly shape description(in character). 
#' See \pkg{ggplot2} for details
#' @param nodes.stroke Numeric. Stroke size of nodes. See details in \pkg{ggplot2}.
#' @param nodes.label.size Numeric. Size of label on nodes.
#' @param nodes.label.colour Character. Colour of label on nodes.
#' @param link.size Numeric. Size of link width.
#' @param link.colour Character. Colour of links, the length should be same as \code{InterCellDB::kpred.action.effect}.
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
#' @importFrom cowplot draw_label
#'
#' @export
#'
GetResult.PlotOnepairClusters.CellPlot <- function(
  VEinfos,
  area.extend.times = 10, 
  hide.locations.A = NULL,
  hide.types.A = NULL,
  hide.locations.B = NULL,
  hide.types.B = NULL,
  hide.sole.vertices = TRUE,  # if TRUE, remove those edges cannot formed vertices
  expand.gap.radius.list = list(ECM = 2, PM = 3, CTP = 2, NC = 2, OTHER = 2),
  expand.shift.degree.list = list(ECM = 90, PM = 90, CTP = 30, NC = 30, OTHER = 30), 
  expand.gap.degree.list = list(ECM = 180, PM = 180, CTP = 60, NC = 60, OTHER = 60),
  expand.center.density.list = list(ECM = 0.25, PM = 0.25, CTP = 0.25, NC = 0.25, OTHER = 0.25),
  expand.outside.cut.percent.list = list(ECM = 0.03, PM = 0.03, CTP = 0.03, NC = 0.03, OTHER = 0.03), 
  nodes.size.range = c(1, 3),  # [TODO] make it change to LogFC
  nodes.fill = "grey",  # [TODO] make it change to anno.Type. Use NULL as default, inside give the Discrete colour mapping
  nodes.colour = "grey",
  nodes.alpha = 1.0,
  nodes.shape = 21,
  nodes.stroke = 0.7,
  nodes.label.size = 3,
  nodes.label.colour = "black",
  link.size = 0.1,
  link.colour = c("red", "green", "blue", "grey"),
  link.alpha = 0.5,
  link.arrow.angle = 20,
  link.arrow.length = 10,
  link.arrow.type = "open",
  caption  # [TODO] add in right corner [Term annotation] PM: Plasma Membrane, ECR: Extracellular Region, ER: Endoplasmic Reticulum, Golgi: Golgi Apparatus"
) {
  ## precheck
  # check scatter parameters are correctly settled
  check.require.items <- c("ECM", "PM", "CTP", "NC", "OTHER")
  tmp.1.check <- length(which(names(expand.gap.radius.list) %in% check.require.items)) == length(check.require.items)
  tmp.2.check <- length(which(names(expand.shift.degree.list) %in% check.require.items)) == length(check.require.items)
  tmp.3.check <- length(which(names(expand.gap.degree.list) %in% check.require.items)) == length(check.require.items)
  tmp.4.check <- length(which(names(expand.center.density.list) %in% check.require.items)) == length(check.require.items)
  if (!(tmp.1.check && tmp.2.check && tmp.3.check && tmp.4.check)) {
    tmp.error.m <- c("expand.gap.radius.list", "expand.shift.degree.list", "expand.gap.degree.list", "expand.center.density.list")
    tmp.error.m <- tmp.error.m[c(!tmp.1.check, !tmp.2.check, !tmp.3.check, !tmp.4.check)]
    stop("The following paramemter: ", paste0(tmp.error.m, collapse = ", "), ", do not include all required cell area, which is ", 
      paste0(check.require.items, collapse = ", "), ". Please check default paramemter value to fix this problem.")
  }
  # check given nodes size range valid?
  if (length(nodes.size.range) != 2 || 
    (length(nodes.size.range == 2) && (nodes.size.range[1] <= 0 || nodes.size.range[2] <= 0))) {
    stop("Given nodes.size.range is unvalid, and please give 2 numeric values (> 0).")
  }
  #
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
  vertices.apx.type.B <- VEinfos$vertices.apx.type.B

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
  if (length(tmp.inds.hide) != 0) {
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
    vertices.infos$UID <- seq_len(nrow(vertices.infos))
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
  pred.loc.special <- c("Plasma Membrane", "Extracellular Region", "Other", "Nucleus")
  # others
  pred.loc.common.A <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.A.clustername)])), pred.loc.special)
  pred.loc.common.B <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.B.clustername)])), pred.loc.special)

  ## A
  tmp.inds.is.A <- which(vertices.infos$ClusterName == act.A.clustername)
  # locations
  this.inds.A.pmem <- intersect(which(vertices.infos$Location == "Plasma Membrane"), tmp.inds.is.A)
  this.inds.A.exm <- intersect(which(vertices.infos$Location == "Extracellular Region"), tmp.inds.is.A)
  this.inds.A.other <- intersect(which(vertices.infos$Location == "Other"), tmp.inds.is.A)
  this.inds.A.nucleus <- intersect(which(vertices.infos$Location == "Nucleus"), tmp.inds.is.A)
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
  this.inds.B.nucleus <- intersect(which(vertices.infos$Location == "Nucleus"), tmp.inds.is.B)
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
  this.plot.A.nucleus <- Inside.GetCoords.ZeroCircle(center.x = -100 * area.extend.times, center.y = 0 * area.extend.times, 
    radius = 14 * area.extend.times)  # 15-1 to be suitable circle radius

  # comp B in +x axis
  this.plot.B.cell <- data.frame(x = c(55, 55, 175, 175) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.pmem <- data.frame(x = c(45, 45, 55, 55) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.exm <- data.frame(x = c(15, 15, 45, 45) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.other <- data.frame(x = c(175, 175, 195, 195) * area.extend.times, y = c(-45, 45, 45, -45) * area.extend.times)
  this.plot.B.nucleus <- Inside.GetCoords.ZeroCircle(center.x = 100 * area.extend.times, center.y = 0 * area.extend.times, 
    radius = 14 * area.extend.times)  # 15-1 to be suitable circle radius

  ## plot the base
  tmpb.other.fill <- "white"; tmpb.other.colour <- "black"; tmpb.other.linetype <- "dotted"
  tmpb.cellx.fill <- "white"; tmpb.cellx.colour <- "black"; tmpb.cellx.linetype <- "solid"
  tmpb.exm.fill <- "lightgrey"; tmpb.exm.alpha <- 0.3; tmpb.exm.colour <- "lightgrey"; tmpb.exm.linetype <- "dotted"
  tmpb.pmemx.fill <- "white"; tmpb.pmemx.colour <- "black"; tmpb.pmemx.linetype <- "longdash"
  tmpb.nucleus.fill <- "lightgrey"; tmpb.nucleus.alpha <- 0.8; tmpb.nucleus.colour <- "grey"; tmpb.nucleus.linetype <- "dashed"
  this.graph.raw.base <- this.base.graph +
      geom_polygon(data = this.plot.A.other, aes(x, y), 
        fill = tmpb.other.fill, colour = tmpb.other.colour, linetype = tmpb.other.linetype) +
      geom_polygon(data = this.plot.A.cell, aes(x, y), 
        fill = tmpb.cellx.fill, colour = tmpb.cellx.colour, linetype = tmpb.cellx.linetype) +
      geom_polygon(data = this.plot.A.exm, aes(x, y), 
        fill = tmpb.exm.fill, alpha = tmpb.exm.alpha, colour = tmpb.exm.colour, linetype = tmpb.exm.linetype) + 
      geom_polygon(data = this.plot.A.pmem, aes(x, y), 
        fill = tmpb.pmemx.fill, colour = tmpb.pmemx.colour, linetype = tmpb.pmemx.linetype) + 
      geom_polygon(data = this.plot.A.nucleus, aes(x = x.cc, y = y.cc),
        fill = tmpb.nucleus.fill, alpha = tmpb.nucleus.alpha, colour = tmpb.nucleus.colour, linetype = tmpb.nucleus.linetype) + 
      geom_polygon(data = this.plot.B.other, aes(x, y), 
        fill = tmpb.other.fill, colour = tmpb.other.colour, linetype = tmpb.other.linetype) +
      geom_polygon(data = this.plot.B.cell, aes(x, y),
        fill = tmpb.cellx.fill, colour = tmpb.cellx.colour, linetype = tmpb.cellx.linetype) +
      geom_polygon(data = this.plot.B.exm, aes(x, y), 
        fill = tmpb.exm.fill, alpha = tmpb.exm.alpha, colour = tmpb.exm.colour, linetype = tmpb.exm.linetype) + 
      geom_polygon(data = this.plot.B.pmem, aes(x, y), 
        fill = tmpb.pmemx.fill, colour = tmpb.pmemx.colour, linetype = tmpb.pmemx.linetype) + 
      geom_polygon(data = this.plot.B.nucleus, aes(x = x.cc, y = y.cc), 
        fill = tmpb.nucleus.fill, alpha = tmpb.nucleus.alpha, colour = tmpb.nucleus.colour, linetype = tmpb.nucleus.linetype)

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
    r.label = rep(c("Cytoplasm & Nucleus", "PM", "ECR", "Other"), each = 2),  # ECR = Extracellular Region, PM = Plasma Membrane
    r.label.nudge.y = rep(c(0, 0, 0, 0), each = 2) * area.extend.times)
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
  ## plot Plasma Membrane
  # A
  this.pmem.A.ctp.xy <- c(-50, 0) * area.extend.times
  this.pmem.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.pmem, ], expand.outside.cut.percent.list$PM, this.pmem.A.ctp.xy, 5 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$PM, 
    sample.shift.degree = expand.shift.degree.list$PM, 
    sample.gap.degree = expand.gap.degree.list$PM, 
    density.half.near = expand.center.density.list$PM, 
    coords.xy.colnames = c("gx", "gy"))
  this.pmem.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.pmem.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 9), rotate.degree = 0, 
    enlarge.xy.ref = this.pmem.A.ctp.xy, rotate.xy.ref = this.pmem.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.A.vx.ext)
  # B
  this.pmem.B.ctp.xy <- c(50, 0) * area.extend.times
  this.pmem.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.pmem, ], expand.outside.cut.percent.list$PM, this.pmem.B.ctp.xy, 5 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$PM, 
    sample.shift.degree = expand.shift.degree.list$PM, 
    sample.gap.degree = expand.gap.degree.list$PM, 
    density.half.near = expand.center.density.list$PM, 
    coords.xy.colnames = c("gx", "gy"))
  this.pmem.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.pmem.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 9), rotate.degree = 0, 
    enlarge.xy.ref = this.pmem.B.ctp.xy, rotate.xy.ref = this.pmem.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.B.vx.ext)

  ## plot Extracellular Region
  # A
  this.exm.A.ctp.xy <- c(-30, 0) * area.extend.times
  this.exm.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.exm, ], expand.outside.cut.percent.list$ECM, this.exm.A.ctp.xy, 10 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$ECM, 
    sample.shift.degree = expand.shift.degree.list$ECM, 
    sample.gap.degree = expand.gap.degree.list$ECM, 
    density.half.near = expand.center.density.list$ECM, 
    coords.xy.colnames = c("gx", "gy"))
  this.exm.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, 
    enlarge.xy.ref = this.exm.A.ctp.xy, rotate.xy.ref = this.exm.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.A.vx.ext)
  # B
  this.exm.B.ctp.xy <- c(30, 0) * area.extend.times
  this.exm.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.exm, ], expand.outside.cut.percent.list$ECM, this.exm.B.ctp.xy, 10 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$ECM, 
    sample.shift.degree = expand.shift.degree.list$ECM, 
    sample.gap.degree = expand.gap.degree.list$ECM, 
    density.half.near = expand.center.density.list$ECM,
    coords.xy.colnames = c("gx", "gy"))
  this.exm.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, 
    enlarge.xy.ref = this.exm.B.ctp.xy, rotate.xy.ref = this.exm.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.B.vx.ext)

  ## plot Nucleus
  # A
  this.nucleus.A.ctp.xy <- c(-100, 0) * area.extend.times
  this.nucleus.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.nucleus, ], expand.outside.cut.percent.list$NC, this.nucleus.A.ctp.xy, 14 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$NC, 
    sample.shift.degree = expand.shift.degree.list$NC, 
    sample.gap.degree = expand.gap.degree.list$NC, 
    density.half.near = expand.center.density.list$NC, 
    coords.xy.colnames = c("gx", "gy"))
  this.nucleus.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.nucleus.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 1), rotate.degree = 0, 
    enlarge.xy.ref = this.nucleus.A.ctp.xy, rotate.xy.ref = this.nucleus.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.nucleus.A.vx.ext)
  # B
  this.nucleus.B.ctp.xy <- c(100, 0) * area.extend.times
  this.nucleus.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.nucleus, ], expand.outside.cut.percent.list$NC, this.nucleus.B.ctp.xy, 14 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$NC, 
    sample.shift.degree = expand.shift.degree.list$NC, 
    sample.gap.degree = expand.gap.degree.list$NC, 
    density.half.near = expand.center.density.list$NC, 
    coords.xy.colnames = c("gx", "gy"))
  this.nucleus.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.nucleus.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 1), rotate.degree = 0, 
    enlarge.xy.ref = this.nucleus.B.ctp.xy, rotate.xy.ref = this.nucleus.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.nucleus.B.vx.ext)
  
  ## plot Other
  # A
  this.other.A.ctp.xy <- c(-185, 0) * area.extend.times
  this.other.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.other, ], expand.outside.cut.percent.list$OTHER, this.other.A.ctp.xy, 10 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$OTHER, 
    sample.shift.degree = expand.shift.degree.list$OTHER, 
    sample.gap.degree = expand.gap.degree.list$OTHER, 
    density.half.near = expand.center.density.list$OTHER, 
    coords.xy.colnames = c("gx", "gy"))
  this.other.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, 
    enlarge.xy.ref = this.other.A.ctp.xy, rotate.xy.ref = this.other.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.A.vx.ext)
  # B
  this.other.B.ctp.xy <- c(185, 0) * area.extend.times
  this.other.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.other, ], expand.outside.cut.percent.list$OTHER, this.other.B.ctp.xy, 10 * area.extend.times, 
    radius.gap.factor = expand.gap.radius.list$OTHER, 
    sample.shift.degree = expand.shift.degree.list$OTHER, 
    sample.gap.degree = expand.gap.degree.list$OTHER, 
    density.half.near = expand.center.density.list$OTHER, 
    coords.xy.colnames = c("gx", "gy"))
  this.other.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = c(1, 4), rotate.degree = 0, 
    enlarge.xy.ref = this.other.B.ctp.xy, rotate.xy.ref = this.other.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.B.vx.ext)

  ## plot common
  gen.coords.put.ps <- function(x, inds.ref, ctp.ref, vertices.infos, radius.gap.factor, sample.shift.degree, sample.gap.degree, density.half.near, outside.cut.percent) {
    this.res <- NULL
    this.inds <- inds.ref[[x]]
    if (length(this.inds) > 0) {
      tmp.ind.m <- match(x, ctp.ref$Map.Items)
      this.ctp <- as.numeric(ctp.ref[tmp.ind.m, c("tpx", "tpy")])
      this.points <- vertices.infos[this.inds, ]
      this.res <- ScatterSimple.Plot(this.points, outside.cut.percent, this.ctp, 14 * area.extend.times,  # 14 is special by allocating from 120 total cytoplasm area width&height 
        radius.gap.factor = radius.gap.factor, coords.xy.colnames = c("gx", "gy"),
        sample.shift.degree = sample.shift.degree, sample.gap.degree = sample.gap.degree, 
        density.half.near = density.half.near)
    }
    this.res
  }

  # A common
  this.common.A.ctp.xy <- data.frame(tpx = rep(c(-70, -100, -130,  -160), each = 3) * area.extend.times,
                     tpy = rep(c(30, 0, -30), times = 4) * area.extend.times)
  # as (-100, 0) being taken by Nucleus, remove it from common ctp.xy
  tmp.rest.inds.A <- setdiff(seq_len(nrow(this.common.A.ctp.xy)), 
    intersect(which(this.common.A.ctp.xy$tpx == -100 * area.extend.times), which(this.common.A.ctp.xy$tpy == 0 * area.extend.times)))
  this.common.A.ctp.xy <- this.common.A.ctp.xy[tmp.rest.inds.A, ]
  # tmp-ly shortcuts  [TODO]
  tmp.plc <- pred.loc.common.A
  length(tmp.plc) <- 12 - 1  # 1 removed as Nucleus seperated itself
  this.common.A.ctp.xy$Map.Items <- tmp.plc
  #
  this.common.A.vx.inlist <- lapply(names(this.inds.A.common),
    inds.ref = this.inds.A.common,
    ctp.ref = this.common.A.ctp.xy,
    vertices.infos = vertices.infos,
    radius.gap.factor = expand.gap.radius.list$CTP, 
    sample.shift.degree = expand.shift.degree.list$CTP, 
    sample.gap.degree = expand.gap.degree.list$CTP, 
    density.half.near = expand.center.density.list$CTP, 
    outside.cut.percent = expand.outside.cut.percent.list$CTP, 
    FUN = gen.coords.put.ps
    )
  this.vx.ext.infos <- rbind(this.vx.ext.infos, bind_rows(this.common.A.vx.inlist))

  # B common
  this.common.B.ctp.xy <- data.frame(tpx = rep(c(70, 100, 130, 160), each = 3) * area.extend.times,
                     tpy = rep(c(30, 0, -30), times = 4) * area.extend.times)
  # as (100, 0) being taken by Nucleus, remove it from common ctp.xy
  tmp.rest.inds.B <- setdiff(seq_len(nrow(this.common.B.ctp.xy)), 
    intersect(which(this.common.B.ctp.xy$tpx == 100 * area.extend.times), which(this.common.B.ctp.xy$tpy == 0 * area.extend.times)))
  this.common.B.ctp.xy <- this.common.B.ctp.xy[tmp.rest.inds.B, ]
  # tmp-ly shortcuts  [TODO]
  tmp.plc <- pred.loc.common.B
  length(tmp.plc) <- 12 - 1  # 1 removed as Nucleus seperated itself
  this.common.B.ctp.xy$Map.Items <- tmp.plc
  #
  this.common.B.vx.inlist <- lapply(names(this.inds.B.common),
    inds.ref = this.inds.B.common,
    ctp.ref = this.common.B.ctp.xy,
    vertices.infos = vertices.infos,
    radius.gap.factor = expand.gap.radius.list$CTP, 
    sample.shift.degree = expand.shift.degree.list$CTP, 
    sample.gap.degree = expand.gap.degree.list$CTP, 
    density.half.near = expand.center.density.list$CTP, 
    outside.cut.percent = expand.outside.cut.percent.list$CTP, 
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
      const.x.dec <- 14 * area.extend.times
      const.y.dec <- 14.5 * area.extend.times
      this.ploygon.xy <- data.frame(x = rep(c(this.cx - const.x.dec, this.cx + const.x.dec), each = 2), 
        y = c(this.cy - const.y.dec, this.cy + const.y.dec, this.cy + const.y.dec, this.cy - const.y.dec))
      tmp.common.polygons <- append(tmp.common.polygons, 
        geom_polygon(data = this.ploygon.xy, aes(x, y),
          fill = "white", colour = "orange", linetype = "dashed"
        ))
    }
    this.graph.add.ps <- this.graph.add.ps + tmp.common.polygons

    tmp.common.ctp.xy[, "label.pos.tpy"] <- tmp.common.ctp.xy[, "tpy"] + 13.5 * area.extend.times
    # use short names

    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Endoplasmic Reticulum"), "Map.Items"] <- "ER"
    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Golgi Apparatus"), "Map.Items"] <- "Golgi"
    tmp.common.ctp.xy[which(tmp.common.ctp.xy[, "Map.Items"] == "Cytoplasm"), "Map.Items"] <- "OtherInPlasm"
    #
    this.graph.add.ps <- this.graph.add.ps + 
      geom_text(data = tmp.common.ctp.xy, 
        mapping = aes(x = tpx, y = label.pos.tpy, label = Map.Items),
        colour = "orange", size = 2)
  }
  # after all regions and their names being settled, get caption drawn
  # [TODO] caption here
  

  ## set points size range by LogFC
  # check validity of given range 
  tmp.linear.size <- (nodes.size.range[2] - nodes.size.range[1]) / (max(this.vx.ext.infos$LogFC) - min(this.vx.ext.infos$LogFC))
  this.vx.ext.infos$nodes.size <- (this.vx.ext.infos$LogFC - min(this.vx.ext.infos$LogFC)) * tmp.linear.size + nodes.size.range[1]

  # -------------
  this.graph.gplot <- this.graph.add.ps
  ## plot vertices
  this.graph.gplot <- this.graph.gplot +
    geom_point(data = this.vx.ext.infos,
      mapping = aes(x = gx, y = gy, size = nodes.size),
      colour = nodes.colour,
      alpha = nodes.alpha,
      shape = nodes.shape,
      stroke = nodes.stroke) + 
    scale_size_identity(name = "LogFC", guide = "legend")
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

  # table exported columns
  tmp.vertices.colnames <- setdiff(colnames(this.vx.ext.infos), c("gx", "gy", "nodes.size"))
  tmp.edges.colnames <- setdiff(colnames(edges.infos), c("from.gx", "from.gy", "to.gx", "to.gy"))
  ## res
  list(plot = this.graph.gplot, 
    table = list(vertices.infos = this.vx.ext.infos[, tmp.vertices.colnames], 
                 edges.infos = edges.infos[, tmp.edges.colnames]))
}




