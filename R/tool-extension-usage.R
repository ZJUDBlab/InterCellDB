


#' Add user-defined dataset
#' 
#' @description
#' This function gives a standlized procedure that package users can easily add its 
#' definitions about genes and their functions or types or locations etc.
#'
#' @param user.def.db Data.frame. It should be a data.frame with at least 2 columns, whose 1st column is gene names, and the 2nd is 
#' definitions in users' personal will, and colums from the 3rd will be inherited intactly.
#' @inheritParams Inside.DummyGenesRefDB
#' @param warning.given Character. It applies extra limitation on the types(molecular functions) of A in gene pairs formatted as A-B.
#'
#' @return Data.frame.
#' | <col-1> <col-2> <col-3> ...|
#'
#' | GeneID Gene.name user.type ...|
#'
#' | <data> <data> <data> ...|
#'
#'
#'
#'
#' @importFrom dplyr left_join
#'
#' @export
#'
Tool.AddUserRestrictDB <- function(
	user.def.db,
	genes.ref.db,
	warning.given = "genes"
) {
	colnames(user.def.db)[1] <- "gene"
	user.def.db <- DataPrep.RemapClustersMarkers(user.def.db, genes.ref.db, warning.given)
	colnames(user.def.db)[1] <- "Gene.name"
	# add GeneID
	entrez.db <- genes.ref.db$gene.ncbi.db
	user.def.db <- left_join(user.def.db, entrez.db[, c("GeneID", "Symbol_from_nomenclature_authority")], 
		by = c("Gene.name" = "Symbol_from_nomenclature_authority"))
	tmp.ncol <- length(colnames(user.def.db))
	user.def.db <- user.def.db[, c(tmp.ncol, 1:(tmp.ncol-1))]
	# return
	user.def.db  # GeneID Gene.name user.type *1 *2 *3 ...
}






#' Find genes annotated in specific GO terms
#' 
#' @description
#' This function uses GO terms or GO IDs to get a specific list of genes.
#'
#'
#' @param go.todolist Character. Several GO_terms or GO_IDs or mixed, which will 
#' be used to get subsets of feature genes. 
#' @inheritParams Inside.DummyGenesRefDB
#' @inheritParams Inside.DummyGORefDB
#' 
#' @return Character. A gene list of given GO IDs or terms.
#'
#'
#'
#' @export
#'
Tool.FindGenesFromGO <- function(
	go.todolist,
	genes.ref.db,
	go.ref.db
) {
	# pre-process
	entrez.ref.db <- genes.ref.db$gene.ncbi.db
	### --- doing all checks ---
	# 1. check if GO_IDs or GO_terms given are available
	# 2. auto transform of all go.todolist to be GO_ID
	# 3. get spliting
	inds.ID.given <- grep("^GO:", go.todolist)
	if (length(inds.ID.given) > 0) {
		go.ID.given.list <- go.todolist[inds.ID.given]
		go.term.given.list <- go.todolist[-inds.ID.given]
	} else {
		go.ID.given.list <- NULL
		go.term.given.list <- go.todolist
	}
	## Giving error report
	# ID matching
	go.ID.given.nonexist <- character()
	go.ID.given.exist <- character()
	if (length(go.ID.given.list) > 0) {
		inds.mID.s <- match(go.ID.given.list, go.ref.db$GO_ID)
		go.ID.given.nonexist <- go.ID.given.list[which(is.na(inds.mID.s))]
		go.ID.given.exist <- go.ID.given.list[which(!is.na(inds.mID.s))]
		if (length(go.ID.given.nonexist) > 0) {
			warning("The following GO_IDs are not found: \n ", paste0(go.ID.given.nonexist, collapse = ", "), ".")
		}
	}
  res.go.id.match.list <- go.ref.db[which(go.ref.db$GO_ID %in% go.ID.given.exist), "GeneName"]
	# term matching
	go.term.given.nonexist <- character()
	go.term.given.exist <- character()
	if (length(go.term.given.list) > 0) { 
		inds.mTerm.s <- match(go.term.given.list, go.ref.db$GO_term)
		go.term.given.nonexist <- go.term.given.list[which(is.na(inds.mTerm.s))]
		go.term.given.exist <- go.term.given.list[which(!is.na(inds.mTerm.s))]
		if (length(go.term.given.nonexist) > 0) {
			warning("The following GO_terms are not found: \n ", paste0(go.term.given.nonexist, collapse = ", "), ".")
		}
	}
  res.go.term.match.list <- go.ref.db[which(go.ref.db$GO_term %in% go.term.given.exist), "GeneName"]
	# finish transformation
	res.go.rel.genes <- unique(as.character(c(res.go.id.match.list, res.go.term.match.list)))
	res.go.rel.genes <- res.go.rel.genes[order(res.go.rel.genes)]
  # return
	res.go.rel.genes
}





#' Draw Circos plot for multiple interaction pairs
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
#' @param if.ignore.annos Logic. It is passed to \code{GenerateVEinfos}. If TRUE, genes with different locations or types documented will
#' be treated as the same, and only one row information will be reserved.
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
Tool.PlotInteractsInMultipleClusters <- function(
  interact.pairs.acted,
  fgenes.remapped.all,
  actions.ref.db,
  return.table.vals = FALSE,
  select.interacts = NULL,
  clusters.select.auto = NULL,
  is.directional = TRUE,
  if.ignore.annos = TRUE,
  sel.mode.val = NULL,
  sel.action.effect.val = NULL,
  show.legend = FALSE,
  cluster.label.cex = 0.6,
  cluster.colour = c("#fb8072", "#80b1d3", "#fdb462", "#ffffb3", "#bebada", "#b3de69", "#fccde5", "#8dd3c7", "#d9d9d9"),
  link.cor.colour = c("grey"),
  link.cor.border.colour = c("grey"),
  link.cor.arrow.drawn = c(TRUE, TRUE, TRUE, FALSE),
  link.cor.arrow.type = c("triangle", "ellipse", "curved", "circle"),
  link.cor.arrow.length = c(0.4, 0.3, 0.3, 0.2),
  ...
) {
  # pre-process set
  pred.mode <- kpred.mode
  pred.action.effect <- kpred.action.effect
  # pre-process check
  if ((length(link.cor.colour) != 1) && (length(link.cor.colour) < length(pred.mode))) {
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
      logic.ifclusters <- clusters.select.auto %in% as.character(unique(unlist(interact.pairs.acted$list.clusters)))
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
    tmp.VEinfos <- GenerateVEinfos(tmp.gmoc, fgenes.remapped.all, is.directional = is.directional, if.ignore.annos = if.ignore.annos)
    tmp.VEinfos <- TrimVEinfos(tmp.VEinfos, sel.mode.val = sel.mode.val, sel.action.effect.val = sel.action.effect.val)
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

