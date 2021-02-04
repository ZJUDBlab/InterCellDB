


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
  ssp.words <- " Please increases the area.extend.times by multiply or increment it by 10."
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
# This function is a random points scatter algorithm used for some area like circle ring
# The return value of this function is the @param data.veinfo plus 2 new columns
# that given by @param coords.xy.colnames
#
ScatterRing.Plot <- function(
  vertices.infos,
  center.xy,  # (x, y)
  ext.len,  # circle, if non-circle is wanted, use coords.trans functions
  ring.radius.range.percent = c(0, 1),  # [1] means the inner radius percent, [2] means the outer one
  radius.gap.factor = 1,
  sample.shift.degree = 0,  # default goes from +x axis to measure degree
  sample.gap.degree = 60,  # default will be 30
  coords.xy.colnames = c("gx", "gy")
) {
  # pre-check
  if (nrow(vertices.infos) < 1) {
    if (!is.null(colnames(vertices.infos))) {
      tmp.ret <- data.frame(ttx = numeric(0), tty = numeric(0))
      colnames(tmp.ret) <- coords.xy.colnames
      return(cbind(vertices.infos, tmp.ret))
    } else {
      return(NULL)
    }
  }
  #
  if (length(coords.xy.colnames) != 2) {
    stop("Unexpected colnames whose length <!=> 2!")
  }
  # get area radius both the near center one and the away one
  radius.near.center <- ceiling(ring.radius.range.percent[1] * ext.len)  # use ceiling to ensure the minimum gets larger
  radius.away.center <- floor(ring.radius.range.percent[2] * ext.len)  # use floor to make sure the maximum gets smaller

  this.puts.cnt <- radius.away.center - radius.near.center + 1 - 2  # count the available points position
  if (this.puts.cnt < 0) {
    stop("Area is too small to locate even one point!")
  }
  # sample
  if (sample.gap.degree < 1 || sample.gap.degree > 360) {
    sample.gap.degree <- 30  # the default
  }
  this.deg.splits <- floor(360 / sample.gap.degree)

  # scatter preparation
  tmp.ring.unit.locs <- floor(this.puts.cnt / radius.gap.factor)  # add 1 to make the number count rather than number substraction
  sp.ring.it <- seq_len(tmp.ring.unit.locs) * radius.gap.factor + radius.near.center
  sp.deg.it <- seq_len(this.deg.splits) - 1
  sp.total.it <- length(sp.ring.it) * length(sp.deg.it)
  # check capacity 
  if (sp.total.it < nrow(vertices.infos) || tmp.ring.unit.locs == 0) {
    stop("Capacity error: Cannot allocate all points!")
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
  # ring coords
  coords.ring.it <- gen.coords.ssp(vertices.infos, sp.deg.it, sp.ring.it, sample.gap.degree, sample.shift.degree, center.xy)
  ## merge back to df
  vertices.infos[, coords.xy.colnames] <- coords.ring.it
  # return
  vertices.infos
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



# This function is to generate coordinates for unit size circle under default parameters. 
Inside.GetCoords.ZeroCircle <- function(
  center.x = 0, 
  center.y = 0, 
  radius = 1.0,
  gradient.degree = 3  # now 120 splits, and in the future, probably use smaller value
) {
  # generating coordinates counter-clockwise
  whole.degree <- 360
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



# This function is to generate coordinates for sector of circle.
# It is used after getting result of \code{Inside.GetCoords.ZeroCircle}.
Inside.GetCoords.PartialCircle <- function(
  uc.all.coords,
  start.degree = 0,
  end.degree = 360,
  gradient.degree = 3  # get it set as it is used in \code{Inside.GetCoords.ZeroCircle}
) {
  # align degree with the closest(floored) neighbor that matches the gradient degree
  start.deg <- floor(start.degree / gradient.degree) * gradient.degree
  end.deg <- floor(end.degree / gradient.degree) * gradient.degree
  while (start.deg < 0) start.deg <- start.deg + 360
  while (start.deg > 360) start.deg <- start.deg - 360
  while (end.deg < 0) end.deg <- end.deg + 360
  while (end.deg > 360) end.deg <- end.deg - 360
  # 
  # grain.level <- as.integer(nrow(uc.all.coords) / floor(360 / gradient.degree))  # if gradient.degree = 1, then get 360 data rows. 
  if (start.deg <= end.deg) {
    ind.start <- start.deg / gradient.degree + 1
    ind.end <- ifelse(end.deg < 360, end.deg + gradient.degree, end.deg) / gradient.degree  # here add ifelse to avoid inds exceeding range error
    vc.part.coords <- uc.all.coords[ind.start:ind.end, ]
  } else {  # it always go counter-clockwise, so if start from larger degree, it goes across the +x axis
    ind.start.1 <- start.deg / gradient.degree + 1
    ind.end.dummy.1 <- 360 / gradient.degree
    vc.part.1.c <- uc.all.coords[ind.start.1:ind.end.dummy.1, ]
    ind.end.2 <- ifelse(end.deg < 360, end.deg + gradient.degree, end.deg) / gradient.degree  # here add ifelse to avoid inds exceeding range error
    vc.part.2.c <- uc.all.coords[1:ind.end.2, ]
    vc.part.coords <- rbind(vc.part.1.c, vc.part.2.c)
  }
  #end# return
  vc.part.coords
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
#' @param hide.sole.vertices Character. It hides sole vertices which cannot form available edges anymore.
#' @param expand.gap.radius.list Numeric. It defines the minimum distance between points(genes) for each plotting area.
#' @param expand.shift.degree.list Numeric. It defines the begin degree that points(genes) are to be drawn. The degree is calculated counter-clockwise.
#' @param expand.gap.degree.list Numeric. It defines the way that points(genes) arrange. If it is set 180 and shift degree is set 90, then points will be aligned in vertical line. 
#' If it is set 90 and shift degree is set 0, then points will be put counter-clockwise from horizontal line to vertical and then back to horizontal and finally at vertical place.
#' @param expand.center.density.list Numeric. It defines the density in each plotting area. Higher value means points concentrating more around 
#' the center of the area, and lower value means more sparse.
#' @param expand.outside.cut.percent.list Numeric. It is used to restrict plotting area range against center point in each plotting area. The value of this parameter 
#' defines the width of outside not-drawing area, and if is set 0, means plotting over the total area, or if is set 0.5, means plotting inside circle of half radius.
#' @param expand.PM.gap.len Numeric. It defines the gap of each points(genes) plotted in area:Plasma Membrane. The larger value means sparser the points being plotted, while smaller denser.
#' @param locate.PM.method Character. It defines the way to allocate points(genes) plotted in area:Plasma Membrane. It supports 2 methods for now: 
#' "random" and "uniform".
#' @param nodes.size.range Numeric of length 2. The former gives the minimum size, while the latter gives the maximum. Node(genes) sizes are reflecting the actual LogFC value of every gene. 
#' @param nodes.size.gap Numeric. It is used along with \code{nodes.size.range}. It defines the resolution of changing point sizes. For example, defines \code{nodes.size.range = c(3,6)} and 
#'  \code{nodes.size.gap = 1}, then point sizes are \code{c(3,4,5,6)}, and will be shown the same in graph legend.
#' @param nodes.fill.updn Character of length 2. It defines the 2 colours used to distinguish the up-regulated and down-regulated nodes(genes). 
#' @param nodes.colour Character. Colour of nodes.
#' @param nodes.alpha Numeric. Alpha of nodes.
#' @param nodes.shape Vector. Use shape reprentable IDs(in integer) or directly shape description(in character). 
#' See \pkg{ggplot2} for details
#' @param nodes.stroke Numeric. Stroke size of nodes. See details in \pkg{ggplot2}.
#' @param label.size.nodes Numeric. It defines the size of label text on nodes(genes). If get length 2, the former defines size of nodes in the left part of graph, and 
#' the latter corresponds to the right part of graph. 
#' @param label.colour.nodes Character. It defines the colour of label.
#' @param label.vjust Numeric. It defines the vertical alignment value for labels.
#' @param label.hjust Numeric. It defines the horizontal alignment value for labels. 
#' @param label.nudge.x Numeric. It defines the slight adjustment movement along x-axis when plotting.
#' @param label.nudge.y Numeric. It defines the slight adjustment movement along y-axis when plotting.
#' @param label.padding.itself Numeric. It defines the padding size of label.
#' @param label.size.itself Numeric. It defines the size of label itself(includes padding, etc).
#' @param link.size Numeric. Size of link width.
#' @param link.colour Character. Colour of links, the length should be same as \code{InterCellDB::kpred.mode}. Otherwise, named vector could be provided 
#' with the desired colours corresponding to some of \code{InterCellDB::kpred.mode}.
#' @param link.alpha Numeric. Alpha of link.
#' @param link.arrow.angle Numeric. Angle of link arrow.
#' @param link.arrow.length Numeric. Length of link arrow.
#' @param link.arrow.type Character. Type of link arrow, either \code{open} or \code{closed}.
#' @param legend.show.fill.updn.label Character of length 2. It gives the labels putting in legend:fill.
#' @param legend.show.fill.override.point.size Numeric. It changes the size of template points used in legend:fill.
#' @param legend.show.size.override.colour Character. It changes the colour of template points used in legend:colour. 
#' @param legend.show.size.override.stroke Numeric. It changes the stroke of template points used in legend:colour.
#' @param legend.show.size.override.size.proportion Numeric. It changes the size of template points used in legend:size, to avoid 
#' unbalanced large or small size shown in legend. 
#' @param legend.manual.left.spacing Unit. It spcifies the left spacing of the manual legend to \pkg{ggplot2} automatically-generated legends.
#' @param legend.manual.internal.spacing Unit. It defines the spacing between 2 manual legend about action mode and action effect.
#'
#'
#'
#' @details
#' All parameters starts with 'label.' can be specified of 2 values or only 1 value. As the graph contains 2 cells(drawn like it), with 1 cell put in the 
#' left side and the other put in right side, there are some situations that number of nodes(genes) is quite unbalanced between them. The label may hard to 
#' be placed properly. To deal with this problem, this function gives the way to specify different label patterns for each part. That's why parameters started with 
#' 'label.' can be set of 2 values as well as 1 value.
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item grid.plot: graphical output generated from \pkg{grid}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#'
#'
#' @import grid
#' @import gtable
#' @import ggplot2
#' @importFrom dplyr bind_rows
#'
#' @export
#'
GetResult.PlotOnepairClusters.CellPlot.SmallData <- function(
  VEinfos,
  area.extend.times = 10, 
  hide.locations.A = c("Other"),
  hide.types.A = NULL,
  hide.locations.B = c("Other"),
  hide.types.B = NULL,
  hide.sole.vertices = TRUE, 
  expand.gap.radius.list = list(ECM = 2, CTP = 2, NC = 2, OTHER = 2),  
  expand.shift.degree.list = list(ECM = 90, CTP = 30, NC = 30, OTHER = 30), 
  expand.gap.degree.list = list(ECM = 180, CTP = 60, NC = 60, OTHER = 60),
  expand.center.density.list = list(ECM = 0.25, NC = 0.25, OTHER = 0.25),
  expand.outside.cut.percent.list = list(ECM = 0.03, NC = 0.03, OTHER = 0.03), 
  expand.PM.gap.len = 2, 
  locate.PM.method = "uniform",  # or "random"
  nodes.size.range = c(6, 12), 
  nodes.size.gap = 1, 
  nodes.fill.updn = c("#D07F86", "#7AAF7A"), 
  nodes.colour = c("lightgrey"), 
  nodes.alpha = 1.0,
  nodes.shape = 21,
  nodes.stroke = 0,
  label.size.nodes = c(3, 3), 
  label.colour.nodes = c("black", "black"),
  label.vjust = c(0, 0),
  label.hjust = c(1, 0),
  label.nudge.x = c(0, 0), 
  label.nudge.y = c(2, 2), 
  label.padding.itself = list(unit(0.25, "lines"), unit(0.25, "lines")),
  label.size.itself = c(0.25, 0.25),
  link.size = 0.6,
  link.colour = c("#D70051", "#00913A", "#1296D4", "#956134", "#C8DC32", "#B5B5B6", "#0A0AFF"), 
  link.alpha = 1, 
  link.linetype = c("solid", "solid", "solid", "44"), 
  link.arrow.angle = c(20, 90, 60, 0), 
  link.arrow.length = c(12, 6, 8, 0), 
  link.arrow.type = c("closed", "open", "open", "open"), 
  legend.show.fill.updn.label = c("UP", "DN"), 
  legend.show.fill.override.point.size = 3, 
  legend.show.size.override.colour = "black", 
  legend.show.size.override.stroke = 0.7, 
  legend.show.size.override.size.proportion = 1,
  legend.manual.left.spacing = unit(0.2, "cm"),
  legend.manual.internal.spacing = unit(0.6, "cm") 
) {
  ## precheck
  # check scatter parameters are correctly settled
  check.require.items <- c("ECM", "CTP", "NC", "OTHER")
  check.require.items.ex <- c("ECM", "NC", "OTHER")
  tmp.1.check <- length(which(names(expand.gap.radius.list) %in% check.require.items)) == length(check.require.items)
  tmp.2.check <- length(which(names(expand.shift.degree.list) %in% check.require.items)) == length(check.require.items)
  tmp.3.check <- length(which(names(expand.gap.degree.list) %in% check.require.items)) == length(check.require.items)
  tmp.4.check <- length(which(names(expand.center.density.list) %in% check.require.items.ex)) == length(check.require.items.ex)
  tmp.5.check <- length(which(names(expand.outside.cut.percent.list) %in% check.require.items.ex)) == length(check.require.items.ex)
  if (!(tmp.1.check && tmp.2.check && tmp.3.check && tmp.4.check)) {
    tmp.error.m <- c("expand.gap.radius.list", "expand.shift.degree.list", "expand.gap.degree.list", "expand.center.density.list", "expand.outside.cut.percent.list")
    tmp.error.m <- tmp.error.m[c(!tmp.1.check, !tmp.2.check, !tmp.3.check, !tmp.4.check)]
    stop("The following paramemter: ", paste0(tmp.error.m, collapse = ", "), ", do not include all required cell area, which is ", 
      paste0(check.require.items, collapse = ", "), ". Please check default paramemter value to fix this problem.")
  }
  # check given nodes size range valid?
  if (length(nodes.size.range) != 2 || 
    (length(nodes.size.range == 2) && (nodes.size.range[1] <= 0 || nodes.size.range[2] <= 0))) {
    stop("Given nodes.size.range is unvalid, and please give 2 numeric values (> 0).")
  }
  ## check for label.* to make it length 2
  if (length(label.size.nodes) < 1 || is.null(label.size.nodes)) stop("Paramemter `label.size.nodes` is invalid!")
  label.size.nodes <- if (length(label.size.nodes) == 1) rep(label.size.nodes, times = 2) else label.size.nodes[1:2]
  if (length(label.colour.nodes) < 1 || is.null(label.colour.nodes)) stop("Paramemter `label.colour.nodes` is invalid!")
  label.colour.nodes <- if (length(label.colour.nodes) == 1) rep(label.colour.nodes, times = 2) else label.colour.nodes[1:2]
  if (length(label.vjust) < 1 || is.null(label.vjust)) stop("Paramemter `label.vjust` is invalid!")
  label.vjust <- if (length(label.vjust) == 1) rep(label.vjust, times = 2) else label.vjust[1:2]
  if (length(label.hjust) < 1 || is.null(label.hjust)) stop("Paramemter `label.hjust` is invalid!")
  label.hjust <- if (length(label.hjust) == 1) rep(label.hjust, times = 2) else label.hjust[1:2]
  if (length(label.nudge.x) < 1 || is.null(label.nudge.x)) stop("Paramemter `label.nudge.x` is invalid!")
  label.nudge.x <- if (length(label.nudge.x) == 1) rep(label.nudge.x, times = 2) else label.nudge.x[1:2]
  if (length(label.nudge.y) < 1 || is.null(label.nudge.y)) stop("Paramemter `label.nudge.y` is invalid!")
  label.nudge.y <- if (length(label.nudge.y) == 1) rep(label.nudge.y, times = 2) else label.nudge.y[1:2]
  if (length(label.padding.itself) < 1 || is.null(label.padding.itself)) stop("Paramemter `label.padding.itself` is invalid!")
  label.padding.itself <- if (length(label.padding.itself) == 1) rep(label.padding.itself, times = 2) else label.padding.itself[1:2]
  if (length(label.size.itself) < 1 || is.null(label.size.itself)) stop("Paramemter `label.size.itself` is invalid!")
  label.size.itself <- if (length(label.size.itself) == 1) rep(label.size.itself, times = 2) else label.size.itself[1:2]
  # check and set link colour
  if (is.null(names(link.colour))) {
    if (length(link.colour) != length(kpred.mode)) {
      length(link.colour) <- length(kpred.mode)
      link.colour[which(is.na(link.colour))] <- "grey"  # use grey to all other
      warning("Given link colour are shorter than expected, and are automatically extended.")
    }
  } else {
    template.link.colour <- c("#D70051", "#00913A", "#1296D4", "#956134", "#C8DC32", "#B5B5B6", "#0A0AFF")
    tmp.inds.link.col <- match(names(link.colour), kpred.mode)
    tmp.link.colour <- link.colour[!is.na(tmp.inds.link.col)]  # get the valid ones
    warning("Named link colour has some unvalid names: ", paste0(link.colour[is.na(tmp.inds.link.col)], collapse = ", "), 
      ", which will be automatically replaced with default colour values.")
    template.link.colour[match(names(tmp.link.colour), kpred.mode)] <- tmp.link.colour
    link.colour <- template.link.colour
  }
  names(link.colour) <- kpred.mode  # set the names right
  
  # check given nodes.size.gap
  if (length(nodes.size.gap) != 1 || nodes.size.gap > abs(nodes.size.range[2] - nodes.size.range[1])) {
    stop("Given `nodes.size.gap` is invalid or too large!")
  }

  #
  act.A.clustername <- VEinfos$cluster.name.A
  act.B.clustername <- VEinfos$cluster.name.B
  edges.infos <- VEinfos$edges.infos
  vertices.infos <- VEinfos$vertices.infos  # infos included: location
  vertices.apx.type.A <- VEinfos$vertices.apx.type.A
  vertices.apx.type.B <- VEinfos$vertices.apx.type.B

  ## defining location properly
  pred.loc.special <- c("Plasma Membrane", "Extracellular Region", "Other", "Nucleus")
  # others
  pred.loc.common.A <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.A.clustername)])), pred.loc.special)
  pred.loc.common.B <- setdiff(levels(factor(vertices.infos$Location[which(vertices.infos$ClusterName == act.B.clustername)])), pred.loc.special)

  ## hide some vertices & merge all cytoplasm locations
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
  # so then the vertices
  vertices.infos <- vertices.infos[setdiff(seq_len(nrow(vertices.infos)), tmp.inds.hide), ]
  # hide sole vertices
  if (hide.sole.vertices == TRUE) {
    tmp.edge.points <- unique(c(edges.infos[, "from"], edges.infos[, "to"]))
    tmp.UID.keep <- intersect(vertices.infos$UID, tmp.edge.points)
    vertices.infos <- vertices.infos[which(vertices.infos[, "UID"] %in% tmp.UID.keep), ]
  }
  ## merge cytoplasm locations
  # A
  tmp.inds.A.cyto.locs <- intersect(which(vertices.infos$ClusterName == act.A.clustername), 
    which(vertices.infos$Location %in% pred.loc.common.A))
  tmp.cor.A.cyto.gnames <- vertices.infos[tmp.inds.A.cyto.locs, "GeneName"]
  tmp.inds.A.cyto.list <- tapply(tmp.inds.A.cyto.locs, tmp.cor.A.cyto.gnames, function(x) { x })
  tmp.inds.A.cyto.getall <- as.integer(unlist(tmp.inds.A.cyto.list))
  tmp.inds.A.cyto.reserve <- as.integer(unlist(lapply(tmp.inds.A.cyto.list, function(x) { x[[1]] })))
  vertices.infos$Location[tmp.inds.A.cyto.reserve] <- rep("Cytoplasm", times = length(tmp.inds.A.cyto.reserve))
  vertices.infos <- vertices.infos[setdiff(seq_len(nrow(vertices.infos)), setdiff(tmp.inds.A.cyto.getall, tmp.inds.A.cyto.reserve)), ]
  # B
  tmp.inds.B.cyto.locs <- intersect(which(vertices.infos$ClusterName == act.B.clustername), 
    which(vertices.infos$Location %in% pred.loc.common.B))
  tmp.cor.B.cyto.gnames <- vertices.infos[tmp.inds.B.cyto.locs, "GeneName"]
  tmp.inds.B.cyto.list <- tapply(tmp.inds.B.cyto.locs, tmp.cor.B.cyto.gnames, function(x) { x })
  tmp.inds.B.cyto.getall <- as.integer(unlist(tmp.inds.B.cyto.list))
  tmp.inds.B.cyto.reserve <- as.integer(unlist(lapply(tmp.inds.B.cyto.list, function(x) { x[[1]] })))
  vertices.infos$Location[tmp.inds.B.cyto.reserve] <- rep("Cytoplasm", times = length(tmp.inds.B.cyto.reserve))
  vertices.infos <- vertices.infos[setdiff(seq_len(nrow(vertices.infos)), setdiff(tmp.inds.B.cyto.getall, tmp.inds.B.cyto.reserve)), ]
  # check if no edges left
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

  ## specifiy locations when ploting
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
  this.inds.A.common <- as.integer(unlist(this.inds.A.common))  # all are in Cytoplasm
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
  this.inds.B.common <- as.integer(unlist(this.inds.B.common))  # all are in Cytoplasm

  ## calculate every area with its plotting parameters
  # basic metric
  # nucleus R: 5, cell-partial R: 10, co-centeralized
  # cell-extra(near 0,0) L: 2, PM L-x-shift(between 2 membrane): 1, PM L-x-len: 2, PM L-y-len: 20
  # outside-ECM L-x: 8
  # outside-OTHER (full align with cell in x axis) L-x-len: 22, L-y-len: 5
  # [TODO] if use program to automatically calculate the most suitable base.gen.it (like area.extend.times)
  # 
  # base.gen.it <- 1
  # 
  tmp.start.point <- c(0, 0)
  basep.ecm.x <- 8; basep.PM.x.shift <- 1; basep.PM.x.len <- 2; basep.cell.extra <- 2; 
  basep.cellp.radius <- 10; basep.nucleus.radius <- 4; 
  basep.other.cellp.dist.in.y <- 1; basep.other.y.len <- 3;
  basep.label.cluster.extra.move.y <- 0.5
  # apply area.extend.times
  basep.ecm.x <- basep.ecm.x * area.extend.times; basep.PM.x.shift <- basep.PM.x.shift * area.extend.times;
  basep.PM.x.len <- basep.PM.x.len * area.extend.times; basep.cell.extra <- basep.cell.extra * area.extend.times;
  #
  basep.cellp.radius <- basep.cellp.radius * area.extend.times; basep.nucleus.radius <- basep.nucleus.radius * area.extend.times;
  basep.other.cellp.dist.in.y <- basep.other.cellp.dist.in.y * area.extend.times; basep.other.y.len <- basep.other.y.len * area.extend.times
  #
  basep.label.cluster.extra.move.y <- basep.label.cluster.extra.move.y * area.extend.times
  ### use base to calc every stick point and center
  ## stick point (in +X direction, it always goes from left-corner, left-top, right-top, right-corner
  ## or goes clockwise)
  # ECM
  this.ecm.stick.x <- c(0, 0, basep.ecm.x, basep.ecm.x)
  this.ecm.stick.y <- basep.cellp.radius * c(-1, 1, 1, -1)
  # PM (use circle coords but has 4 stick points like rect)
  this.pm.stick.x <- basep.ecm.x + basep.PM.x.len + c(rep(0, times = 2), rep(basep.PM.x.shift, times = 2))
  this.pm.stick.y <- basep.cellp.radius * c(-1, 1, 1, -1)
  # --- some sum ---
  tmp.start.cell.p.x <- basep.ecm.x + basep.PM.x.len + basep.PM.x.shift
  # CELL [Exception:goes from left-top, clockwise]
  # (in +X direction, it leaves out the left line of the rect)
  this.cell.path.x <- c(tmp.start.cell.p.x, tmp.start.cell.p.x + 2 * basep.cellp.radius, tmp.start.cell.p.x + 2 * basep.cellp.radius, tmp.start.cell.p.x)
  this.cell.path.y <- basep.cellp.radius * c(1, 1, -1, -1)
  # NUCLEUS
  # Nucleus is circle, no stick x y
  # OTHER
  this.other.stick.x <- this.cell.path.x  # align with cell, so the stick.y gets slightly different to the rest
  this.other.stick.y <- 0 - basep.cellp.radius - basep.other.cellp.dist.in.y - basep.other.y.len / 2 + basep.other.y.len / 2 * c(1, 1, -1, -1)

  ## center
  this.ecm.cxy <- c(basep.ecm.x / 2, 0)
  # PM needs no center to scatter
  this.cell.cxy <- c(tmp.start.cell.p.x + basep.cell.extra + basep.cellp.radius - basep.cell.extra, 0)
  this.nucleus.cxy <- this.cell.cxy
  this.other.cxy <- c(this.cell.cxy[1] - basep.cell.extra / 2, 0 - basep.cellp.radius - basep.other.cellp.dist.in.y - basep.other.y.len / 2)

  ## scatter radius
  this.ecm.sct.rad <- min(basep.ecm.x / 2, basep.cellp.radius)
  # PM doesn't scatter
  this.cell.sct.rad.outer <- basep.cellp.radius
  this.cell.sct.rad.inner <- basep.nucleus.radius  # as use ring scatter method
  this.nucleus.sct.rad <- basep.nucleus.radius
  this.other.sct.rad <- min(basep.cell.extra + 2 * basep.cellp.radius, basep.other.y.len) / 2

  ## scatter transform parameter
  this.ecm.trans <- c(basep.ecm.x / 2, basep.cellp.radius) / this.ecm.sct.rad
  # PM doesn't scatter
  # CELL use ring method, no need to tranform
  this.nucleus.trans <- c(basep.nucleus.radius, basep.nucleus.radius) / this.nucleus.sct.rad
  this.other.trans <- c(basep.cell.extra + 2 * basep.cellp.radius, basep.other.y.len) / 2 / this.other.sct.rad

  ## rest things
  # PM curves
  tmp.ref.near.c <- c(basep.ecm.x + basep.PM.x.len, 0)
  tmp.ref.away.c <- c(basep.ecm.x + basep.PM.x.len + basep.PM.x.shift, 0)
  tmp.uc.near.c <- Inside.GetCoords.ZeroCircle(center.x = tmp.ref.near.c[1], center.y = tmp.ref.near.c[2], radius = basep.PM.x.len, gradient.degree = 1)
  tmp.uc.away.c <- Inside.GetCoords.ZeroCircle(center.x = tmp.ref.away.c[1], center.y = tmp.ref.away.c[2], radius = basep.PM.x.len, gradient.degree = 1)
  #- near.c
  this.PM.curve.near.c <- Inside.GetCoords.PartialCircle(tmp.uc.near.c, start.degree = 90, end.degree = 270, gradient.degree = 1)
  this.PM.curve.near.c <- Inside.TransCoords.Enlarge.Rotate(this.PM.curve.near.c, c(basep.PM.x.len, basep.cellp.radius) / basep.PM.x.len, 0, 
    enlarge.xy.ref = tmp.ref.near.c)
  this.PM.curve.near.c[, 1] <- rev(this.PM.curve.near.c[, 1])  # x in near.c
  this.PM.curve.near.c[, 2] <- rev(this.PM.curve.near.c[, 2])  # y in near.c
  #- away.c
  this.PM.curve.away.c <- Inside.GetCoords.PartialCircle(tmp.uc.away.c, start.degree = 90, end.degree = 270, gradient.degree = 1)
  this.PM.curve.away.c <- Inside.TransCoords.Enlarge.Rotate(this.PM.curve.away.c, c(basep.PM.x.len, basep.cellp.radius) / basep.PM.x.len, 0, 
    enlarge.xy.ref = tmp.ref.away.c)
  # label (on cluster)
  this.label.cluster.dxy <- c(this.cell.cxy[1], basep.cellp.radius + basep.label.cluster.extra.move.y)
  # vertical split line
  this.vertical.line.extra.mul.to.cellp.radius <- 0.12
  tmp.vertical.y.val <- basep.cellp.radius * (1 + this.vertical.line.extra.mul.to.cellp.radius)
  this.vertical.plot.data <- data.frame(x = rep(0, times = 2), y = tmp.vertical.y.val * c(-1, 1))

  ### BASE plot
  this.base.graph <- ggplot()
  # coords of plotting compts (must have)
  # comp A in -x axis
  # CELL path
  this.plot.A.cell.path <- data.frame(x = (-1) * this.cell.path.x, y = this.cell.path.y)
  # CELL base
  this.plot.A.cell.base <- data.frame(x = (-1) * c(this.pm.stick.x[4], this.PM.curve.away.c$x.cc, this.pm.stick.x[3], this.cell.path.x[2:3]), 
    y = c(this.pm.stick.y[4], this.PM.curve.away.c$y.cc, this.pm.stick.y[3], this.cell.path.y[2:3]))
  # NUCLEUS
  this.plot.A.nucleus <- Inside.GetCoords.ZeroCircle(center.x = (-1) * this.nucleus.cxy[1], center.y = this.nucleus.cxy[2], radius = this.nucleus.sct.rad)
  # PM get special for add 2 curve between stick points
  this.plot.A.pmem <- data.frame(x = (-1) * c(this.pm.stick.x[1], this.PM.curve.near.c$x.cc, this.pm.stick.x[2:3], this.PM.curve.away.c$x.cc, this.pm.stick.x[4]),
    y = c(this.pm.stick.y[1], this.PM.curve.near.c$y.cc, this.pm.stick.y[2:3], this.PM.curve.away.c$y.cc, this.pm.stick.y[4]))
  # ECM not drawn
  
  # comp B in +x axis
  # CELL path
  this.plot.B.cell.path <- data.frame(x = this.cell.path.x, y = this.cell.path.y)
  # CELL base
  this.plot.B.cell.base <- data.frame(x = c(this.pm.stick.x[4], this.PM.curve.away.c$x.cc, this.pm.stick.x[3], this.cell.path.x[2:3]), 
    y = c(this.pm.stick.y[4], this.PM.curve.away.c$y.cc, this.pm.stick.y[3], this.cell.path.y[2:3]))
  # NUCLEUS
  this.plot.B.nucleus <- Inside.GetCoords.ZeroCircle(center.x = this.nucleus.cxy[1], center.y = this.nucleus.cxy[2], radius = this.nucleus.sct.rad)
  # PM get special for add 2 curve between stick points
  this.plot.B.pmem <- data.frame(x = c(this.pm.stick.x[1], this.PM.curve.near.c$x.cc, this.pm.stick.x[2:3], this.PM.curve.away.c$x.cc, this.pm.stick.x[4]), 
    y = c(this.pm.stick.y[1], this.PM.curve.near.c$y.cc, this.pm.stick.y[2:3], this.PM.curve.away.c$y.cc, this.pm.stick.y[4]))
  # ECM not drawn

  # [HIDE] area OTHER is plotted optionally
  # this.plot.A.other <- data.frame(x = (-1) * this.other.stick.x, y = this.other.stick.y)
  # this.plot.B.other <- data.frame(x = this.other.stick.x, y = this.other.stick.y)

  ## plot the base
  #tmpb.other.fill <- "white"; tmpb.other.colour <- "black"; tmpb.other.linetype <- "dotted"
  #tmpb.cellx.colour <- "black"; tmpb.cellx.linetype <- "solid"
  tmpb.cellbase.fill <- "#F8F5FA"; tmpb.cellbase.alpha <- 1; tmpb.cellbase.colour <- "#F8F5FA"; tmpb.cellbase.linetype <- "solid"
  tmpb.pmemx.fill <- "#B3AFC5"; tmpb.pmemx.alpha <- 1; tmpb.pmemx.colour <- "#B3AFC5"; tmpb.pmemx.linetype <- "solid"
  tmpb.nucleus.fill <- "#9A92AE"; tmpb.nucleus.alpha <- 1; tmpb.nucleus.colour <- "#9A92AE"; tmpb.nucleus.linetype <- "solid"
  # extra: vertical split line
  tmpex.vline.colour <- "lightgrey"; tmpex.vline.linetype <- "dashed"
  this.graph.raw.base <- this.base.graph +
      # [HIDE] other
      #geom_polygon(data = this.plot.A.other, aes(x, y), 
      #  fill = tmpb.other.fill, colour = tmpb.other.colour, linetype = tmpb.other.linetype) +
      geom_polygon(data = this.plot.A.cell.base, aes(x, y),
        fill = tmpb.cellbase.fill, alpha = tmpb.cellbase.alpha, colour = tmpb.cellbase.colour, linetype = tmpb.cellbase.linetype) + 
      #geom_path(data = this.plot.A.cell.path, aes(x, y), 
      #  colour = tmpb.cellx.colour, linetype = tmpb.cellx.linetype) +
      geom_polygon(data = this.plot.A.pmem, aes(x, y), 
        fill = tmpb.pmemx.fill, alpha = tmpb.pmemx.alpha, colour = tmpb.pmemx.colour, linetype = tmpb.pmemx.linetype) + 
      geom_polygon(data = this.plot.A.nucleus, aes(x = x.cc, y = y.cc),
        fill = tmpb.nucleus.fill, alpha = tmpb.nucleus.alpha, colour = tmpb.nucleus.colour, linetype = tmpb.nucleus.linetype) + 
      # [HIDE] other
      #geom_polygon(data = this.plot.B.other, aes(x, y), 
      #  fill = tmpb.other.fill, colour = tmpb.other.colour, linetype = tmpb.other.linetype) +
      geom_polygon(data = this.plot.B.cell.base, aes(x, y),
        fill = tmpb.cellbase.fill, alpha = tmpb.cellbase.alpha, colour = tmpb.cellbase.colour, linetype = tmpb.cellbase.linetype) + 
      #geom_path(data = this.plot.B.cell.path, aes(x, y),
      #  colour = tmpb.cellx.colour, linetype = tmpb.cellx.linetype) +
      geom_polygon(data = this.plot.B.pmem, aes(x, y), 
        fill = tmpb.pmemx.fill, alpha = tmpb.pmemx.alpha, colour = tmpb.pmemx.colour, linetype = tmpb.pmemx.linetype) + 
      geom_polygon(data = this.plot.B.nucleus, aes(x = x.cc, y = y.cc), 
        fill = tmpb.nucleus.fill, alpha = tmpb.nucleus.alpha, colour = tmpb.nucleus.colour, linetype = tmpb.nucleus.linetype) + 
      # extra: vertical split line
      geom_path(data = this.vertical.plot.data, aes(x, y), 
        colour = tmpex.vline.colour, linetype = tmpex.vline.linetype)

  # label the cluster
  this.label.clusters <- data.frame(
    x.lc = c(-1 * this.label.cluster.dxy[1], this.label.cluster.dxy[1]), 
    y.lc = rep(this.label.cluster.dxy[2], times = 2),
    r.label.c = c(act.A.clustername, act.B.clustername))
  this.graph.raw.base <- this.graph.raw.base + 
      geom_text(data = this.label.clusters, 
        mapping = aes(x = x.lc, y = y.lc, label = r.label.c),
        colour = "black", size = 4, 
        vjust = 0, hjust = 0.5)


  ### plot points
  this.graph.add.ps <- this.graph.raw.base
  this.vx.ext.infos <- data.frame()
  ## plot Plasma Membrane
  # fetch available coords of points
  tmp.pmem.y.lim <- c(this.pm.stick.y[1], this.pm.stick.y[2])
  tmp.pmem.seeker <- this.pm.stick.y[1]
  tmp.pmem.inds.seek <- integer()
  for (ip in seq_along(this.PM.curve.near.c$y.cc)) {
    tmp.it <- this.PM.curve.near.c$y.cc[ip]
    if (tmp.it < tmp.pmem.seeker) {
      next
    } else {
      tmp.pmem.inds.seek <- c(tmp.pmem.inds.seek, ip)
      tmp.pmem.seeker <- tmp.pmem.seeker + expand.PM.gap.len
    }
  }
  # remove the first and last few percentage to avoid edges effect, here remove 1 more in last direction as the following algorithm 
  # uses biased 1:max, which got 0 position lacked, so 1 more removed from last ones
  tmp.pmem.inds.seek <- tmp.pmem.inds.seek[setdiff(seq_along(tmp.pmem.inds.seek), c(1, (length(tmp.pmem.inds.seek) - 1):length(tmp.pmem.inds.seek)))]
  if (length(tmp.pmem.inds.seek) == 0) {
    stop("Error: could not allocate even 1 gene as setting too big gap for placing genes in PM.")
  }
  if (length(tmp.pmem.inds.seek) < max(this.inds.A.pmem, this.inds.B.pmem)) {
    stop("Cannot be located inside PM! Please increases the area.extend.times by multiply or increment it by 10.")
  }
  this.pmem.all.avp <- data.frame(gx = this.PM.curve.near.c$x.cc[tmp.pmem.inds.seek], gy = this.PM.curve.near.c$y.cc[tmp.pmem.inds.seek])# available in +x
  ## to plot
  # A 
  this.pmem.A.vx.ext <- vertices.infos[this.inds.A.pmem, ]
  this.pmem.A.vx.ext[, c("gx", "gy")] <- switch(locate.PM.method, 
    "random" = this.pmem.all.avp[sample(seq_len(nrow(this.pmem.all.avp)), length(this.inds.A.pmem)), c("gx", "gy")],
    "uniform" = this.pmem.all.avp[floor((nrow(this.pmem.all.avp) / length(this.inds.A.pmem)) * seq_along(this.inds.A.pmem)), c("gx", "gy")], 
    stop("Undefined method, only `random` and `uniform` are allowed!"))
  this.pmem.A.vx.ext$gx <- (-1) * this.pmem.A.vx.ext$gx
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.A.vx.ext)
  # B
  this.pmem.B.vx.ext <- vertices.infos[this.inds.B.pmem, ]
  this.pmem.B.vx.ext[, c("gx", "gy")] <- switch(locate.PM.method, 
    "random" = this.pmem.all.avp[sample(seq_len(nrow(this.pmem.all.avp)), length(this.inds.B.pmem)), c("gx", "gy")],
    "uniform" = this.pmem.all.avp[floor((nrow(this.pmem.all.avp) / length(this.inds.B.pmem)) * seq_along(this.inds.B.pmem)), c("gx", "gy")], 
    stop("Undefined method, only `random` and `uniform` are allowed!"))
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.pmem.B.vx.ext)

  ## plot Extracellular Region
  # A
  this.exm.A.ctp.xy <- c(-1 * this.ecm.cxy[1], this.ecm.cxy[2])
  this.exm.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.exm, ], 
    expand.outside.cut.percent.list$ECM, this.exm.A.ctp.xy, this.ecm.sct.rad, 
    radius.gap.factor = expand.gap.radius.list$ECM, 
    sample.shift.degree = expand.shift.degree.list$ECM, 
    sample.gap.degree = expand.gap.degree.list$ECM, 
    density.half.near = expand.center.density.list$ECM, 
    coords.xy.colnames = c("gx", "gy"))
  this.exm.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = this.ecm.trans, rotate.degree = 0, 
    enlarge.xy.ref = this.exm.A.ctp.xy, rotate.xy.ref = this.exm.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.A.vx.ext)
  # B
  this.exm.B.ctp.xy <- this.ecm.cxy
  this.exm.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.exm, ], 
    expand.outside.cut.percent.list$ECM, this.exm.B.ctp.xy, this.ecm.sct.rad, 
    radius.gap.factor = expand.gap.radius.list$ECM, 
    sample.shift.degree = expand.shift.degree.list$ECM, 
    sample.gap.degree = expand.gap.degree.list$ECM, 
    density.half.near = expand.center.density.list$ECM,
    coords.xy.colnames = c("gx", "gy"))
  this.exm.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.exm.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = this.ecm.trans, rotate.degree = 0, 
    enlarge.xy.ref = this.exm.B.ctp.xy, rotate.xy.ref = this.exm.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.exm.B.vx.ext)

  ## plot Nucleus
  # A
  this.nucleus.A.ctp.xy <- c(-1 * this.nucleus.cxy[1], this.nucleus.cxy[2])
  this.nucleus.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.nucleus, ], 
    expand.outside.cut.percent.list$NC, this.nucleus.A.ctp.xy, this.nucleus.sct.rad, 
    radius.gap.factor = expand.gap.radius.list$NC, 
    sample.shift.degree = expand.shift.degree.list$NC, 
    sample.gap.degree = expand.gap.degree.list$NC, 
    density.half.near = expand.center.density.list$NC, 
    coords.xy.colnames = c("gx", "gy"))
  this.nucleus.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.nucleus.A.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = this.nucleus.trans, rotate.degree = 0, 
    enlarge.xy.ref = this.nucleus.A.ctp.xy, rotate.xy.ref = this.nucleus.A.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.nucleus.A.vx.ext)
  # B
  this.nucleus.B.ctp.xy <- this.nucleus.cxy
  this.nucleus.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.nucleus, ], 
    expand.outside.cut.percent.list$NC, this.nucleus.B.ctp.xy, this.nucleus.sct.rad, 
    radius.gap.factor = expand.gap.radius.list$NC, 
    sample.shift.degree = expand.shift.degree.list$NC, 
    sample.gap.degree = expand.gap.degree.list$NC, 
    density.half.near = expand.center.density.list$NC, 
    coords.xy.colnames = c("gx", "gy"))
  this.nucleus.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.nucleus.B.vx.ext[, c("gx", "gy")], 
    enlarge.xy.times = this.nucleus.trans, rotate.degree = 0, 
    enlarge.xy.ref = this.nucleus.B.ctp.xy, rotate.xy.ref = this.nucleus.B.ctp.xy)
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.nucleus.B.vx.ext)
  
  ## plot Other [HIDE]
  ## A
  #this.other.A.ctp.xy <- c(-1 * this.other.cxy[1], this.other.cxy[2])
  #this.other.A.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.A.other, ], 
  #  expand.outside.cut.percent.list$OTHER, this.other.A.ctp.xy, this.other.sct.rad, 
  #  radius.gap.factor = expand.gap.radius.list$OTHER, 
  #  sample.shift.degree = expand.shift.degree.list$OTHER, 
  #  sample.gap.degree = expand.gap.degree.list$OTHER, 
  #  density.half.near = expand.center.density.list$OTHER, 
  #  coords.xy.colnames = c("gx", "gy"))
  #this.other.A.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.A.vx.ext[, c("gx", "gy")], 
  #  enlarge.xy.times = this.other.trans, rotate.degree = 0, 
  #  enlarge.xy.ref = this.other.A.ctp.xy, rotate.xy.ref = this.other.A.ctp.xy)
  #this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.A.vx.ext)
  ## B
  #this.other.B.ctp.xy <- this.other.cxy
  #this.other.B.vx.ext <- ScatterSimple.Plot(vertices.infos[this.inds.B.other, ], 
  #  expand.outside.cut.percent.list$OTHER, this.other.B.ctp.xy, this.other.sct.rad, 
  #  radius.gap.factor = expand.gap.radius.list$OTHER, 
  #  sample.shift.degree = expand.shift.degree.list$OTHER, 
  #  sample.gap.degree = expand.gap.degree.list$OTHER, 
  #  density.half.near = expand.center.density.list$OTHER, 
  #  coords.xy.colnames = c("gx", "gy"))
  #this.other.B.vx.ext[, c("gx", "gy")] <- Inside.TransCoords.Enlarge.Rotate(this.other.B.vx.ext[, c("gx", "gy")], 
  #  enlarge.xy.times = this.other.trans, rotate.degree = 0, 
  #  enlarge.xy.ref = this.other.B.ctp.xy, rotate.xy.ref = this.other.B.ctp.xy)
  #this.vx.ext.infos <- rbind(this.vx.ext.infos, this.other.B.vx.ext)

  ## plot Cytoplasm
  # A
  this.cyto.A.ctp.xy <- c(-1 * this.cell.cxy[1], this.cell.cxy[2])
  this.cyto.A.vx.ext <- ScatterRing.Plot(vertices.infos[this.inds.A.common, ], 
    this.cyto.A.ctp.xy, this.cell.sct.rad.outer, 
    ring.radius.range.percent = c(this.cell.sct.rad.inner, this.cell.sct.rad.outer) / this.cell.sct.rad.outer,
    radius.gap.factor = expand.gap.radius.list$CTP, 
    sample.shift.degree = expand.shift.degree.list$CTP, 
    sample.gap.degree = expand.gap.degree.list$CTP, 
    coords.xy.colnames = c("gx", "gy"))
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.cyto.A.vx.ext)
  # B
  this.cyto.B.ctp.xy <- this.cell.cxy
  this.cyto.B.vx.ext <- ScatterRing.Plot(vertices.infos[this.inds.B.common, ], 
    this.cyto.B.ctp.xy, this.cell.sct.rad.outer, 
    ring.radius.range.percent = c(this.cell.sct.rad.inner, this.cell.sct.rad.outer) / this.cell.sct.rad.outer, 
    radius.gap.factor = expand.gap.radius.list$CTP, 
    sample.shift.degree = expand.shift.degree.list$CTP, 
    sample.gap.degree = expand.gap.degree.list$CTP, 
    coords.xy.colnames = c("gx", "gy"))
  this.vx.ext.infos <- rbind(this.vx.ext.infos, this.cyto.B.vx.ext)  


  ## set points size range by LogFC
  tmp.linear.size <- (nodes.size.range[2] - nodes.size.range[1]) / (max(abs(this.vx.ext.infos$LogFC)) - min(abs(this.vx.ext.infos$LogFC)))
  this.vx.ext.infos$nodes.size <- (abs(this.vx.ext.infos$LogFC) - min(abs(this.vx.ext.infos$LogFC))) * tmp.linear.size + nodes.size.range[1]
  
  ## add nodes fill by LogFC UP (>0) or DN (<0)
  this.vx.ext.infos$UPDN <- as.character(unlist(lapply(this.vx.ext.infos$LogFC, function(x) {
    ifelse(x >= 0, legend.show.fill.updn.label[1], legend.show.fill.updn.label[2])
    })))
  this.updn.fill <- nodes.fill.updn
  names(this.updn.fill) <- legend.show.fill.updn.label[1:2]

  # -------------
  this.graph.gplot <- this.graph.add.ps
  ## plot vertices
  this.graph.gplot <- this.graph.gplot +
    geom_point(data = this.vx.ext.infos,
      mapping = aes(x = gx, y = gy, fill = UPDN, size = nodes.size),
      colour = nodes.colour,
      alpha = nodes.alpha,
      shape = nodes.shape,
      stroke = nodes.stroke) + 
    scale_fill_manual(name = "Exprs Change", values = this.updn.fill, breaks = names(this.updn.fill), aesthetics = "fill", 
      guide = guide_legend(override.aes = list(size = legend.show.fill.override.point.size),
        order = 1)) + 
    scale_size_continuous(name = "LogFC", range = c(nodes.size.range[1], nodes.size.range[2]), 
      breaks = seq(from = nodes.size.range[1], to = nodes.size.range[2], by = nodes.size.gap), 
      guide = guide_legend(override.aes = list(colour = legend.show.size.override.colour, 
        stroke = legend.show.size.override.stroke,
        size = legend.show.size.override.size.proportion * seq(from = nodes.size.range[1], to = nodes.size.range[2], by = nodes.size.gap),
        order = 2)))

  ## split the label function to be cluster specific
  # A
  tmp.plot.label.A <- geom_label(data = this.vx.ext.infos[which(this.vx.ext.infos$ClusterName == act.A.clustername), ],
      mapping = aes(x = gx, y = gy, label = GeneName),
      size = label.size.nodes[[1]],
      colour = label.colour.nodes[[1]],
      vjust = label.vjust[[1]], hjust = label.hjust[[1]], 
      nudge_x = label.nudge.x[[1]], nudge_y = label.nudge.y[[1]], 
      label.padding = label.padding.itself[[1]],
      label.size = label.size.itself[[1]])
  # B
  tmp.plot.label.B <- geom_label(data = this.vx.ext.infos[which(this.vx.ext.infos$ClusterName == act.B.clustername), ],
      mapping = aes(x = gx, y = gy, label = GeneName),
      size = label.size.nodes[[2]],
      colour = label.colour.nodes[[2]],
      vjust = label.vjust[[2]], hjust = label.hjust[[2]], 
      nudge_x = label.nudge.x[[2]], nudge_y = label.nudge.y[[2]], 
      label.padding = label.padding.itself[[2]],
      label.size = label.size.itself[[2]]) 
  #
  this.graph.gplot <- this.graph.gplot + list(tmp.plot.label.A, tmp.plot.label.B)

  # -------------
  ## plot edges
  inds.e.from.match <- match(edges.infos$from, this.vx.ext.infos$UID)
  inds.e.to.match <- match(edges.infos$to, this.vx.ext.infos$UID)
  edges.infos[, c("from.gx", "from.gy")] <- this.vx.ext.infos[inds.e.from.match, c("gx", "gy")]
  edges.infos[, c("to.gx", "to.gy")] <- this.vx.ext.infos[inds.e.to.match, c("gx", "gy")]
  #
  show.mode.val <- levels(factor(edges.infos$mode))
  show.action.effect.val <- levels(factor(edges.infos$action.effect))
  # check if global variables cover all action modes and effects
  if (sum(show.mode.val %in% kpred.mode) != length(show.mode.val)) {
    warning("Database changes with new action mode: ", paste0(setdiff(show.mode.val, kpred.mode), collapse = ", "), "!")
  }
  if (sum(show.action.effect.val %in% kpred.action.effect) != length(show.action.effect.val)) {
    warning("Database changes with new action effect: ", paste0(setdiff(show.action.effect.val, kpred.action.effect), collapse = ", "), "!")
  }

  # draw edges
  # use reverse order to plot positive & negative above the unspecified or undirected
  for (j in rev(1:length(kpred.action.effect))) {
    if (!(kpred.action.effect[j] %in% show.action.effect.val)) {
      next
    }
    this.edges.infos <- edges.infos[which(edges.infos[, "action.effect"] == kpred.action.effect[j]), ]
    this.graph.gplot <- this.graph.gplot + geom_segment(data = this.edges.infos,
      mapping = aes(x = from.gx, y = from.gy, xend = to.gx, yend = to.gy, colour = mode),
      alpha = link.alpha, size = link.size, linetype = link.linetype[j], 
      arrow = arrow(angle = link.arrow.angle[j], length = unit(link.arrow.length[j], "pt"), type = link.arrow.type[j]))
  }
  this.graph.gplot <- this.graph.gplot + 
    scale_colour_manual(guide = FALSE, values = link.colour, breaks = names(link.colour), aesthetics = "colour")

  # other settings
  this.graph.gplot <- this.graph.gplot + 
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme(panel.background = element_blank())

  # --- then goes to grid, no more ggplot ---
  this.gGrob <- ggplotGrob(this.graph.gplot)
  ### create manual legend for action mode and action effect
  ## for action mode, for at least 7 types valid given in database, and add some extra types manually
  ## action mode only need one line and one colour, so can be easily generated
  plot.lwd.L.mode <- 1
  plot.grey.bg.L.mode <- polygonGrob(x = c(0, 1, 1, 0),
      y = c(1, 1, 0, 0),
      gp = gpar(
        col = NA, 
        fill = "lightgrey",
        alpha = 0.4)
    )
  tmp.legend.mode <- kpred.mode[which(kpred.mode %in% show.mode.val)]
  LT.m.mode.list <- lapply(seq_along(tmp.legend.mode), 
    tmp.legend.mode = tmp.legend.mode, 
    link.colour.mode = link.colour, link.alpha.mode = link.alpha, 
    plot.lwd.L.mode = plot.lwd.L.mode, plot.grey.bg.L.mode = plot.grey.bg.L.mode, 
    function(x, tmp.legend.mode, link.colour.mode, link.alpha.mode, plot.lwd.L.mode, plot.grey.bg.L.mode) {
        this.L <- grobTree(
          plot.grey.bg.L.mode, 
          polylineGrob(x = c(.1, .9),
            y = c(.5, .5),
            id = c(1, 1)),
          gp = gpar(
            col = link.colour.mode[x],
            alpha = link.alpha.mode, 
            lwd = plot.lwd.L.mode)
        )
        this.T <- textGrob(tmp.legend.mode[x], x = .1, y = .5, just = "left", 
          gp = gpar(fontsize = 10))
        list(L = this.L, T = this.T)
      })
  # legend key and label
  L.mode.to.use <- lapply(LT.m.mode.list, function(x) { x$L })
  T.mode.to.use <- lapply(LT.m.mode.list, function(x) { x$T })
  # legend title
  title.mode <- textGrob("Action Mode", x = .0, y = .5, just = "left")

  ## for action effect, only 4 types
  # legend key
  plot.lwd.L.act.eff <- 1
  plot.grey.bg.L.act.eff <- polygonGrob(x = c(0, 1, 1, 0),
      y = c(1, 1, 0, 0),
      gp = gpar(
        col = NA, 
        fill = "lightgrey",
        alpha = 0.4)
    )
  L.act.eff.1 <- grobTree(
      plot.grey.bg.L.act.eff, 
      polylineGrob(x = c(.1, .9), 
        y = c(.5, .5),
        id = c(1, 1)),
      polygonGrob(x = c(.9, .6, .6),
        y = c(.5, .4, .6)),
      gp = gpar(
        #col = "#D70051",  # for positive
        col = "black", 
        # fill = "#D70051",  # for positive
        fill = "black", 
        alpha = link.alpha, 
        lwd = plot.lwd.L.act.eff)
    )
  L.act.eff.2 <- grobTree(
      plot.grey.bg.L.act.eff,
      polylineGrob(x = c(.1, .85, .85, .85, .85, .85),
        y = c(.5, .5, .5, .35, .5, .65),
        id = c(1, 1, 2, 2, 3, 3)),
      gp = gpar(
        #col = "#00913A",  # for negative
        col = "black", 
        alpha = link.alpha, 
        lwd = plot.lwd.L.act.eff)
    )
  L.act.eff.3 <- grobTree(
      plot.grey.bg.L.act.eff, 
      polylineGrob(x = c(.1, .9, .9, .8, .9, .8),
        y = c(0.5, 0.5, .5, .3, .5, .7),
        id = c(1, 1, 2, 2, 3, 3)),  # [NOTE] point in 0:25 is 75% size as given in size (measure in diameter)
      #pointsGrob(x = 0.8, y = 0.5, size = unit(0.266, "npc"), pch = 16, default.units = "npc"),
      gp = gpar(
        #col = "#956134",  # for unspecified
        col = "black", 
        alpha = link.alpha, 
        lwd = plot.lwd.L.act.eff)
    )
  L.act.eff.4 <- grobTree(
      plot.grey.bg.L.act.eff, 
      linesGrob(x = c(.1, .9),
        y = c(.5, .5)),
      gp = gpar(
        #col = "#B5B5B6",  # for undirected
        col = "black", 
        alpha = link.alpha, 
        lwd = plot.lwd.L.act.eff,
        lty = "22")
    )
  # legend label
  T.act.eff.1 <- textGrob("positive", x = .1, y = .5, just = "left", gp = gpar(fontsize = 10))
  T.act.eff.2 <- textGrob("negative", x = .1, y = .5, just = "left", gp = gpar(fontsize = 10))
  T.act.eff.3 <- textGrob("unspecified", x = .1, y = .5, just = "left", gp = gpar(fontsize = 10))
  T.act.eff.4 <- textGrob("undirected", x = .1, y = .5, just = "left", gp = gpar(fontsize = 10))
  # legend title
  title.act.eff <- textGrob("Action Effect", x = .0, y = .5, just = "left")
  # pack all action effect elements up (except the title)
  L.act.eff.list <- list(L.act.eff.1, L.act.eff.2, L.act.eff.3, L.act.eff.4)
  T.act.eff.list <- list(T.act.eff.1, T.act.eff.2, T.act.eff.3, T.act.eff.4)
  # choose the existing ones to be used in legend
  exist.inds.act.eff <- which(kpred.action.effect %in% show.action.effect.val)
  L.act.eff.to.use <- L.act.eff.list[exist.inds.act.eff]
  T.act.eff.to.use <- T.act.eff.list[exist.inds.act.eff]

  ## create gtable for action legend
  # get fit width for text
  table.fit.width.for.text <- max(unlist(lapply(T.act.eff.to.use, function(x) {
      convertWidth(grobWidth(x), "cm", TRUE)
    })),
    unlist(lapply(T.mode.to.use, function(x) {
      convertWidth(grobWidth(x), "cm", TRUE)
    }))
  )
  # create gtable
  table.act.width  <- c(.6, table.fit.width.for.text / 0.9)
  table.act.height.mode <- unit(c(.6, rep(1, times = length(L.mode.to.use)) * 0.6), "cm")
  table.act.gap <- legend.manual.internal.spacing
  if (class(legend.manual.internal.spacing) != "unit") {
    if (!is.numeric(legend.manual.internal.spacing)) {
      stop("Non-numeric value given to specify spacing! Please check paramemter 'legend.manual.internal.spacing'!")
    } else {
      warning("Given spacing is not class:unit. Auto change to unit in 'cm'.")
      table.act.gap <- unit(legend.manual.internal.spacing, "cm")
    }
  }
  table.act.height.act.eff <- unit(c(.6, rep(1, times = length(L.act.eff.to.use)) * 0.6), "cm")
  leg.m.act <- gtable(width = unit(table.act.width, "cm"), 
    height = unit.c(table.act.height.mode, table.act.gap, table.act.height.act.eff))
  ## add all grobs into it
  tmp.pos.at.mode <- 1
  # add act mode legend title
  leg.m.act <- gtable_add_grob(leg.m.act, title.mode, t = tmp.pos.at.mode, l = 1, r = 2)
  # add act mode rest grobs
  for (i in seq_along(L.mode.to.use)) {
    leg.m.act <- gtable_add_grob(leg.m.act, L.mode.to.use[[i]], t = tmp.pos.at.mode + i, l = 1)
    leg.m.act <- gtable_add_grob(leg.m.act, T.mode.to.use[[i]], t = tmp.pos.at.mode + i, l = 2)
  }
  # action effect
  tmp.pos.at.act.eff <- 1 + length(L.mode.to.use) + 2
  # add act effect legend title
  leg.m.act <- gtable_add_grob(leg.m.act, title.act.eff, t = tmp.pos.at.act.eff, l = 1, r = 2)
  # add act effect rest grobs
  for (i in seq_along(L.act.eff.to.use)) {
    leg.m.act <- gtable_add_grob(leg.m.act, L.act.eff.to.use[[i]], t = tmp.pos.at.act.eff + i, l = 1)
    leg.m.act <- gtable_add_grob(leg.m.act, T.act.eff.to.use[[i]], t = tmp.pos.at.act.eff + i, l = 2)
  }


  # merge ggplot2 and manual legend
  mleg.start.pos <- this.gGrob$layout[which(this.gGrob$layout$name == "guide-box"), c("t", "l")]
  this.gGrob <- gtable_add_cols(this.gGrob, sum(leg.m.act$widths), mleg.start.pos$l)
  this.gGrob <- gtable_add_grob(this.gGrob, leg.m.act, t = mleg.start.pos$t, l = mleg.start.pos$l + 1)
  this.gGrob <- gtable_add_cols(this.gGrob, legend.manual.left.spacing, mleg.start.pos$l)

  # table exported columns
  tmp.vertices.colnames <- setdiff(colnames(this.vx.ext.infos), c("gx", "gy", "nodes.size"))
  tmp.edges.colnames <- setdiff(colnames(edges.infos), c("from.gx", "from.gy", "to.gx", "to.gy"))

  ## res
  list(plot = NULL, grid.plot = this.gGrob, 
    table = list(vertices.infos = this.vx.ext.infos[, tmp.vertices.colnames], 
                 edges.infos = edges.infos[, tmp.edges.colnames]))
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
#' @param expand.outside.cut.percent.list Numeric. It is used to restrict plotting area range against center point in each plotting area. The value of this parameter 
#' defines the width of outside not-drawing area, and if is set 0, means plotting over the total area, or if is set 0.5, means plotting inside circle of half radius.
#' @param nodes.size.range Numeric of length 2. The former gives the minimum size, while the latter gives the maximum. Node(genes) sizes are reflecting the actual LogFC value of every gene. 
#' @param nodes.fill Character. Fill of nodes.
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
#'
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
GetResult.PlotOnepairClusters.CellPlot.LargeData <- function(
  VEinfos,
  area.extend.times = 10, 
  hide.locations.A = NULL,
  hide.types.A = NULL,
  hide.locations.B = NULL,
  hide.types.B = NULL,
  hide.sole.vertices = TRUE,
  expand.gap.radius.list = list(ECM = 2, PM = 3, CTP = 2, NC = 2, OTHER = 2),
  expand.shift.degree.list = list(ECM = 90, PM = 90, CTP = 30, NC = 30, OTHER = 30), 
  expand.gap.degree.list = list(ECM = 180, PM = 180, CTP = 60, NC = 60, OTHER = 60),
  expand.center.density.list = list(ECM = 0.25, PM = 0.25, CTP = 0.25, NC = 0.25, OTHER = 0.25),
  expand.outside.cut.percent.list = list(ECM = 0.03, PM = 0.03, CTP = 0.03, NC = 0.03, OTHER = 0.03), 
  nodes.size.range = c(1, 3), 
  nodes.fill = "grey", 
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
  link.arrow.type = "open"
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
  # [TODO] caption here, not to added for now
  # add in right corner [Term annotation] PM: Plasma Membrane, ECR: Extracellular Region, ER: Endoplasmic Reticulum, Golgi: Golgi Apparatus"
  

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
    scale_colour_manual(name = "Action Effect", values = link.colour, breaks = names(link.colour), aesthetics = "colour")
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





