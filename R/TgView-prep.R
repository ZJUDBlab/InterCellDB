

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
	if (!(cluster.name.C %in% as.character(unlist(interact.pairs.acted$list.clusters)) &&
			cluster.name.D %in% as.character(unlist(interact.pairs.acted$list.clusters)))) {  # check if given cluster.name.* are in clusters.name list
		stop("Error: Target one-pair clusters not found, as giving undefined name of clusters.")
	}
	this.clusters.name <- c(cluster.name.C, cluster.name.D)
	this.sel.name <- paste0(cluster.name.C, kClustersSplit, cluster.name.D) 
	## slim - locations
	this.loc.slim.tg.cols <- c("Gene.name", "GO.Term.target")
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
	act.B.genename,
	put.colnames
) {
	res.mode <- onerow.info["mode"]
	res.actionid <- 1
	ifisdir <- ifelse(onerow.info["is_directional"] == 't', TRUE, FALSE)
	if (!ifisdir) {
		if (onerow.info["a_is_acting"] == 't') {
			stop("Exception: database is broken by mannually modification, as is_directional is f but a_is_acting is t!")
		} # else res.actionid <- 1 (stay default)
	} else {
		ifconv <- ifelse(onerow.info[put.colnames] == "conv", 0, 1)
		# in database, under is_directional = t, 
		# if a_is_acting = t, confirm that it is A-act-upon-B
		# if a_is_acting = f, in oppsite way, B-act-upon-A
		ifa.act <- ifelse(onerow.info["a_is_acting"] == 't', 0, 1)
		# get final offset
		if.op <- (ifconv + ifa.act) %% 2
		if (onerow.info["action"] == "") {
			res.actionid <- match("A--oB", kpred.ext.action.effect) + if.op
		} else {
			if (onerow.info["action"] == "activation") {  # the database record positive actions in "activation"
				res.actionid <- match("A-->B", kpred.ext.action.effect) + if.op
			} else {
				if (onerow.info["action"] == "inhibition") {  # the database record negative actions in "inhibition"
					res.actionid <- match("A--|B", kpred.ext.action.effect) + if.op
				} else {
					stop("Exception: database is broken by mannually modification, as undefined values appear in column('action')!")
				}
			}
		}
	}
	#end# return
	c(as.character(res.mode), as.character(res.actionid), as.character(onerow.info["score"]))
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
	
	# 
	this.act.conv.pairs <- left_join(bt.pairs, actions.ref.db[, setdiff(colnames(actions.ref.db), c("inter.GeneID.A", "inter.GeneID.B"))], 
		by = c("inter.GeneName.A", "inter.GeneName.B"))
	# use rev actions to join again
	this.act.rev.pairs <- left_join(bt.pairs, actions.ref.db[, setdiff(colnames(actions.ref.db), c("inter.GeneID.A", "inter.GeneID.B"))], 
		by = c("inter.GeneName.A" = "inter.GeneName.B", "inter.GeneName.B" = "inter.GeneName.A"))
	# additional infos for action_id
	tmp.tried.indicator.colname <- "__Indicator__"
	for (tried.indicator.i in 1:100) {
		if ((tmp.tried.indicator.colname %in% colnames(this.act.conv.pairs)) == TRUE) {
			tmp.tried.indicator.colname <- paste0(tmp.tried.indicator.colname, i)
		} else {
			break
		}
		if (tried.indicator.i == 100) {
			stop("GenerateMapDetailOnepairClusters has no way to generate different column names that are not in the original data!")
		}
	}
	this.act.conv.pairs[, tmp.tried.indicator.colname] <- "conv"
	this.act.rev.pairs[, tmp.tried.indicator.colname] <- "rev"
	#
	this.act.all.pairs <- rbind(this.act.conv.pairs, this.act.rev.pairs)
	# get those with mode defined collected
	this.act.put.pairs <- this.act.all.pairs[which(!is.na(this.act.all.pairs[, "mode"])), ]
	# those without mode specified, it will has replicates in conv and rev
	this.act.pv.pairs <- this.act.all.pairs[which(is.na(this.act.all.pairs[, "mode"])), ]
	this.act.pv.pairs <- DoPartUnique(this.act.pv.pairs, which(colnames(this.act.pv.pairs) %in% c("inter.GeneName.A", "inter.GeneName.B")))
	# rough selected put & preserved pairs
	this.put.pairs <- this.act.put.pairs
	this.pv.pairs <- this.act.pv.pairs

	# further selection upon [put pairs] by action_id pattern
	this.put.sc.ind.list <- list()
	this.pv.sc.ind.list <- list()
	this.put.short.cut <- paste(this.put.pairs[, "inter.GeneID.A"], this.put.pairs[, "inter.GeneID.B"], sep = ">")
	this.put.sc.ind.list <- tapply(seq_len(length(this.put.short.cut)), this.put.short.cut, function(x) {x})
	# further selection upon [pv pairs] by overlap with [put pairs]
	this.pv.short.cut <- paste(this.pv.pairs[, "inter.GeneID.A"], this.pv.pairs[, "inter.GeneID.B"], sep = ">")
	this.pv.sc.ind.list <- tapply(seq_len(length(this.pv.short.cut)), this.pv.short.cut, function(x) {x})
	# collect overlap ones and further remove these
	tmp.overlap <- intersect(names(this.put.sc.ind.list), names(this.pv.sc.ind.list))
	tmp.inds.overlap <- as.integer(unlist(lapply(tmp.overlap, total.list = this.pv.sc.ind.list, function(x, total.list){
		total.list[[x]]
		})))
	this.pv.pairs <- this.pv.pairs[setdiff(seq_len(nrow(this.pv.pairs)), tmp.inds.overlap), ]

	
	#
	prog.bar.gmoc <- progress::progress_bar$new(total = length(this.put.sc.ind.list))
	prog.bar.gmoc$tick(0)
	this.put.act.detailed <- lapply(this.put.sc.ind.list, 
		tmp.prog = prog.bar.gmoc, this.put.pairs = this.put.pairs, put.colnames = tmp.tried.indicator.colname, 
		function(x, tmp.prog, this.put.pairs, put.colnames) {
			tmp.put.pairs <- this.put.pairs[x, ]
			tmp.put.act.infos <- apply(tmp.put.pairs, MARGIN = 1,
				act.A.genename = tmp.put.pairs[1, "inter.GeneName.A"],
				act.B.genename = tmp.put.pairs[1, "inter.GeneName.B"],
				put.colnames = put.colnames, 
				function(x, act.A.genename, act.B.genename, put.colnames) {
					Inside.CollectActionMapping(x, act.A.genename, act.B.genename, put.colnames)
					})
			tmp.put.act.infos <- t(tmp.put.act.infos)
			## data trimming
			# deliminate those 1s when higher IDs exist, which means if directional is determined, leave out old ambiguous ones.
			# see human database F3->F7, then this situation will be understood
			tmp.put.act.infos.df <- data.frame(
				mode = as.character(tmp.put.act.infos[, 1]), 
				actionid = as.integer(tmp.put.act.infos[, 2]),
				score = as.integer(tmp.put.act.infos[, 3]), 
				stringsAsFactors = FALSE
			)
			# to examine if more specific mode has been recorded, that is the same mode but more detailed action type
			inds.del.1 <- which(tmp.put.act.infos.df$actionid == 1)
			inds.rest <- which(tmp.put.act.infos.df$actionid != 1)
			logic.del.1 <- tmp.put.act.infos.df[inds.del.1, "mode"] %in% tmp.put.act.infos.df[inds.rest, "mode"]
			tmp.put.act.infos.df.exm <- rbind(tmp.put.act.infos.df[inds.del.1[which(logic.del.1 == FALSE)], ],
										 tmp.put.act.infos.df[inds.rest, ]
									 )
			tmp.put.act.infos.df.exm <- unique(tmp.put.act.infos.df.exm)
			# tick
			tmp.prog$tick()
			# return
			list(
				act.A.genename = tmp.put.pairs[1, "inter.GeneName.A"],  # use first row as it has at least 1 row
				act.B.genename = tmp.put.pairs[1, "inter.GeneName.B"],
				act.A.logfc = tmp.put.pairs[1, "inter.LogFC.A"],
				act.B.logfc = tmp.put.pairs[1, "inter.LogFC.B"],
				action.infos = tmp.put.act.infos.df.exm
			)
	})

	# further pack up upon [pv pairs]
	this.pv.act.detailed <- list()
	prog.bar.sub.pv <- progress::progress_bar$new(total = nrow(this.pv.pairs))
	prog.bar.sub.pv$tick(0)
	this.pv.act.detailed <- lapply(seq_len(nrow(this.pv.pairs)), this.pv.pairs = this.pv.pairs, 
		function(x, this.pv.pairs) {
		prog.bar.sub.pv$tick()
		list(act.A.genename = this.pv.pairs[x, "inter.GeneName.A"],
				 act.B.genename = this.pv.pairs[x, "inter.GeneName.B"],
				 act.A.logfc = this.pv.pairs[x, "inter.LogFC.A"],
				 act.B.logfc = this.pv.pairs[x, "inter.LogFC.B"],
				 action.infos = data.frame(mode = "other", actionid = 1, score = NA, stringsAsFactors = FALSE)
				)
		})  

	bt.pairs.result$actions.detailed <- c(this.put.act.detailed, this.pv.act.detailed)
	#end# return
	bt.pairs.result
}





#' Generate data about vertices and edges
#'
#' @description
#' This function uses detailed informations about one interaction pair(return value of 
#' \code{GenerateMapDetailOnepairClusters()}), to generate data for drawing relation plot.
#'
#' @param onepair.gmoc List. Return value of \code{\link{GenerateMapDetailOnepairClusters}}.
#' @inheritParams Inside.DummyFgenes 
#' @param direction.X.to.Y [TODO]
#' @param if.ignore.location Logic. Logic. It is passed to \code{GenerateVEinfos}. If TRUE, genes with different locations or types documented will
#' be treated as the same, and only one row information will be reserved.
#'
#' @details
#' This function uses actions that are recorded in STRING act database, but only a small part of 
#' actions are thoroughly difined in the database.
#' This function is used to generate formatted data structure(with vertices and edges).
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
	direction.X.to.Y = NULL,
	if.ignore.location = FALSE
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
	# do unique if locations and types are not cared
	if (if.ignore.location == TRUE) {
		vertices.all.infos <- DoPartUnique(vertices.all.infos, cols.select = match(c("GeneName", "ClusterName"), colnames(vertices.all.infos)))
	}
	vertices.all.infos$UID <- 1:nrow(vertices.all.infos)
	rownames(vertices.all.infos) <- NULL
	# change colnames in vertices.all
	tmp.cols.change <- match(c("GO.Term.target", "LogFC", "PVal"), colnames(vertices.all.infos))
	colnames(vertices.all.infos)[tmp.cols.change] <- c("Location", "LogFC", "PVal")
	tmp.cols.first6 <- c("UID", "ClusterName", "GeneName", "Location", "LogFC", "PVal")
	vertices.all.infos <- vertices.all.infos[, c(tmp.cols.first6, setdiff(colnames(vertices.all.infos), tmp.cols.first6))]  # rearrange the columns
	# change colnames in apx.*
	colnames(vertices.A.apx.types) <- colnames(vertices.B.apx.types) <- c("GeneName", "Type")

	## --- edges ---
	# predefined function
	gen.edges.vei.inside <- function(act.part1.UID, act.part2.UID, action.mode, action.effect, action.score) {
		# this function is to generate all permutation of act.part1.UID ~ act.part2.UID, e.g. A*2 B*3 will get 2*3 results
		tmp.all.pert <- lapply(act.part1.UID,
			act.part2.UID = act.part2.UID, action.mode = action.mode, action.effect = action.effect, action.score = action.score, 
			FUN = function(x, act.part2.UID, action.mode, action.effect, action.score) {
				data.frame(from = x, to = act.part2.UID, 
					action.mode = action.mode, action.effect = action.effect, action.score = action.score, 
					stringsAsFactors = FALSE)
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
						if (is.null(direction.X.to.Y) || direction.X.to.Y == TRUE) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "undirected", this.row["score"]))  
						}
						if (is.null(direction.X.to.Y) || !direction.X.to.Y) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "undirected", this.row["score"]))
						}
					} else {
						if (this.row["actionid"] < 2 || this.row["actionid"] > 7) {
							stop(paste0("Undefined actionid from @param onepair.gmoc$actions.detailed[[", i, "]]!"))
						}
						if (this.row["actionid"] == 2 && (is.null(direction.X.to.Y) || direction.X.to.Y == TRUE)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "positive", this.row["score"]))
						} 
						if (this.row["actionid"] == 3 && (is.null(direction.X.to.Y) || !direction.X.to.Y)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "positive", this.row["score"]))
						}
						if (this.row["actionid"] == 4 && (is.null(direction.X.to.Y) || direction.X.to.Y == TRUE)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "negative", this.row["score"]))
						}
						if (this.row["actionid"] == 5 && (is.null(direction.X.to.Y) || !direction.X.to.Y)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "negative", this.row["score"]))
						}
						if (this.row["actionid"] == 6 && (is.null(direction.X.to.Y) || direction.X.to.Y == TRUE)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.A.UID, act.B.UID, this.row["mode"], "unspecified", this.row["score"]))
						}
						if (this.row["actionid"] == 7 && (is.null(direction.X.to.Y) || !direction.X.to.Y)) {
							tmp.act.res <- c(tmp.act.res, gen.edges.vei.inside(act.B.UID, act.A.UID, this.row["mode"], "unspecified", this.row["score"]))
						}
					}
				}
			}
			prog.bar.edge.collect$tick()
			bind_rows(tmp.act.res)
	})
	edges.all.infos <- as.data.frame(bind_rows(this.act.result))
	#end# return
	VEinfos <- list(cluster.name.A = afterV.A.clustername, cluster.name.B = afterV.B.clustername,
		edges.infos = edges.all.infos, 
		vertices.infos = vertices.all.infos,
		vertices.apx.type.A = vertices.A.apx.types,
		vertices.apx.type.B = vertices.B.apx.types
		)
	return(VEinfos)
}


# as .mirror for A-A cluster group, sometimes it is need to removed it 
getOrigClusterNameTgVEInfo <- function(
	object
) {
	this.veinfos <- object@tg.veinfo
	cluster.name.B <- this.veinfos$cluster.name.B
	ind.match <- grep("mirror$", cluster.name.B)
	if (length(ind.match) > 0) {
		cluster.name.B <- strsplit(cluster.name.B, split = ".mirror", fixed = TRUE)[[1]][1]
	}
	list(cluster.name.A = this.veinfos$cluster.name.A, cluster.name.B = cluster.name.B)
}





# This function is to fetch interactions in given 2-cell groups
FetchInterOI <- function(
	object,
	cluster.x,
	cluster.y,
	if.ignore.location = FALSE
) {
	# check cluster.x cluster.y are embeded in `ExtractTargetOnepairClusters()`
	# process
	tg.inter <- ExtractTargetOnepairClusters(getFullViewResult(object), cluster.x, cluster.y)
	object <- setTgActionPairs(object, GenerateMapDetailOnepairClusters(tg.inter, object@database@actions.db))
	object <- setTgVEInfo(object, 
		GenerateVEinfos(getTgActionPairs(object), object@fgenes, 
			direction.X.to.Y = NULL,  # keep all in default setting
			if.ignore.location)
	)

	# return
	return(object)
}





#' Select part of gene pairs
#'
#' @description
#' This function is to further select subset of interactions in given 2-cell group.
#'
#' @param object [TODO]
#' @param direction.X.to.Y [TODO]
#' @param sel.exprs.change Character. It selects the expression change status that gene pairs can be. It has total 4 options:
#' "Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn", which are defined in global variables \code{kexprs.change}.
#' @param sel.some.genes.X Character. It selects some genes in cluster A(A is defined in parameter \code{VEinfos} and use \code{VEinfos$cluster.name.A} to fetch that).
#' @param sel.some.genes.Y Character. It selects some genes in cluster B(A is defined in parameter \code{VEinfos} and use \code{VEinfos$cluster.name.B} to fetch that).
#' @param sel.genes.option Character. Its allowed values are "intersect" and "union". It defines the way that merges the result of subset selected by some given genes.  
#' The default way is to intersect the subsets, and the other option is union.
#' @param sel.gene.pairs Data.frame. It is at-least-2-column data.frame, which records gene pairs with each column settling 
#' one of the participated genes. The 2 required column names need to be specified by another parameter \code{sel.some.gene.pairs.colnames}.
#' @param sel.some.gene.pairs.colnames Character of length 2. It strictly specifies the column names that records the genes of given gene pairs, 
#' and it also implies the direction goes from the first column to the second (AtoB). So, make sure putting genes in their proper positions.
#' @param sel.mode.val Character. If set NULL, it uses all values in global variables \code{InterCellDB::kpred.mode}, or
#' please specify detailed and accurate values in subset of \code{InterCellDB::kpred.mode}.
#' @param sel.action.effect.val Character. If set NULL, it uses all values in global variables \code{InterCellDB::kpred.action.effect}, or
#' please specify detailed and accurate values in subset of \code{InterCellDB::kpred.action.effect}.
#' @param sel.mode.action.option Character. Its allowed values are "intersect" and "union". It defines the way that merges the result of subset selected by mode 
#' and subset selected by action.effect. The default way is to intersect the subsets, and the other option is union.
#'
#'
#' @details
#' The whole list of mode or action.effect is defined as global variable that is given within the package.
#' kpred.mode defines 9 modes, while kpred.action.effect defines 4 action effects.
#' As the action database given now is not completed well, it is recommended to select upon
#' action.effect, while leaving all different modes preserved.
#'
#'
#'
#' @importFrom dplyr left_join 
#'
#' @export
#'
SelectInterSubset <- function(
	object, 
	direction.X.to.Y = NULL, 
	sel.exprs.change = c("Xup.Yup", "Xup.Ydn", "Xdn.Yup", "Xdn.Ydn"), 
	sel.some.genes.X = NULL, 
	sel.some.genes.Y = NULL, 
	sel.genes.option = "intersect", 
	sel.gene.pairs = NULL, 
	sel.some.gene.pairs.colnames = c("inter.GeneName.A", "inter.GeneName.B"), 
	sel.mode.val = NULL, 
	sel.action.effect.val = NULL, 
	sel.mode.action.option = "intersect" 
) {
	VEinfos <- getTgVEInfo(object)
	#
	afterV.A.clustername <- VEinfos$cluster.name.A
	afterV.B.clustername <- VEinfos$cluster.name.B
	vertices.all.infos <- VEinfos$vertices.infos
	edges.all.infos <- VEinfos$edges.infos
	vertices.A.apx.types <- VEinfos$vertices.apx.type.A
	vertices.B.apx.types <- VEinfos$vertices.apx.type.B
	### select target edges.part.infos and vertices.part.infos by some genes
	if (!is.null(sel.some.genes.X) || !is.null(sel.some.genes.Y)) {
		# --- merge option ---
		if(!(sel.genes.option[[1]] %in% c("union", "intersect"))) {  # only use [1]
			stop("Sel-some-genes merge option error: undefined merge options! only 'union' and 'intersect' are allowed.")
		}
		uid.somegenes.sel.A <- uid.somegenes.sel.B <- integer()
		if (is.null(sel.some.genes.X)) {
			uid.somegenes.sel.A <- vertices.all.infos[which(vertices.all.infos$ClusterName == afterV.A.clustername), "UID"]
		} else {
			uid.somegenes.sel.A <- vertices.all.infos[intersect(which(vertices.all.infos$ClusterName == afterV.A.clustername), which(vertices.all.infos$GeneName %in% sel.some.genes.X)), "UID"]    
		}
		if (is.null(sel.some.genes.Y)) {
			uid.somegenes.sel.B <- vertices.all.infos[which(vertices.all.infos$ClusterName == afterV.B.clustername), "UID"]
		} else {
			uid.somegenes.sel.B <- vertices.all.infos[intersect(which(vertices.all.infos$ClusterName == afterV.B.clustername), which(vertices.all.infos$GeneName %in% sel.some.genes.Y)), "UID"]    
		}
		inds.somegenes.result <- integer()
		if (sel.genes.option == "intersect") {
			# conv
			inds.somegenes.conv <- intersect(which(edges.all.infos$from %in% uid.somegenes.sel.A), which(edges.all.infos$to %in% uid.somegenes.sel.B))
			# rev
			inds.somegenes.rev <- intersect(which(edges.all.infos$from %in% uid.somegenes.sel.B), which(edges.all.infos$to %in% uid.somegenes.sel.A))
			inds.somegenes.result <- unique(c(inds.somegenes.conv, inds.somegenes.rev))
		} else {
			if (sel.genes.option == "union") {
				# conv
				inds.somegenes.conv <- union(which(edges.all.infos$from %in% uid.somegenes.sel.A), which(edges.all.infos$to %in% uid.somegenes.sel.B))
				# rev
				inds.somegenes.rev <- union(which(edges.all.infos$from %in% uid.somegenes.sel.B), which(edges.all.infos$to %in% uid.somegenes.sel.A))
				inds.somegenes.result <- unique(c(inds.somegenes.conv, inds.somegenes.rev))
			}
		}
		edges.sel0.infos <- edges.all.infos[inds.somegenes.result, ]
	} else {
		edges.sel0.infos <- edges.all.infos
	}

	### select target edges.part.infos and vertices.part.infos by some gene pairs
	## As the process in selecting mode & action.effect using 'edges' as reference, So here, only selecting 'edges' is enough
	if (!is.null(sel.gene.pairs)) {
		if (class(sel.gene.pairs) == "data.frame" && nrow(sel.gene.pairs) > 0 && ncol(sel.gene.pairs) >= 2) {
			tmp.x <- sum(sel.some.gene.pairs.colnames %in% colnames(sel.gene.pairs))
			if (tmp.x != length(sel.some.gene.pairs.colnames) || tmp.x < 2) {
				stop("Given colnames in `sel.some.gene.pairs.colnames` are unvalid or less than 2 (unable to define gene pairs)! They are ",
					paste0(sel.some.gene.pairs.colnames, collapse = ", "))
			}
			tmp.sel.gp.df <- sel.gene.pairs[, sel.some.gene.pairs.colnames]
			tmp.p1.by.vec <- "GeneName"
			names(tmp.p1.by.vec) <- sel.some.gene.pairs.colnames[1]
			tmp.p2.by.vec <- "GeneName"
			names(tmp.p2.by.vec) <- sel.some.gene.pairs.colnames[2]
			# join is considering to be cluster-gene perfectly matched
			tob.res.sel.gp.df <- left_join(sel.gene.pairs, vertices.all.infos[, c("UID", "GeneName")], by = tmp.p1.by.vec)
			colnames(tob.res.sel.gp.df)[which(colnames(tob.res.sel.gp.df) %in% c("UID"))] <- "part1.UID"
			tob.res.sel.gp.df <- left_join(tob.res.sel.gp.df, vertices.all.infos[, c("UID", "GeneName")], by = tmp.p2.by.vec)
			colnames(tob.res.sel.gp.df)[which(colnames(tob.res.sel.gp.df) %in% c("UID"))] <- "part2.UID"
			# get subset, using short-inter-pair to match
			tob.res.short.interacts.conv <- paste(tob.res.sel.gp.df[, "part1.UID"], tob.res.sel.gp.df[, "part2.UID"], sep = "->")
			tob.res.short.interacts.rev <- paste(tob.res.sel.gp.df[, "part2.UID"], tob.res.sel.gp.df[, "part1.UID"], sep = "->")
			tmp.edges.short.interacts <- paste(edges.sel0.infos[, "from"], edges.sel0.infos[, "to"], sep = "->")
			edges.sel1.infos <- edges.sel0.infos[which(tmp.edges.short.interacts %in% unique(c(tob.res.short.interacts.conv, tob.res.short.interacts.rev))), ]
		} else {
			stop("Given parameter `sel.gene.pairs` should be data.frame that has 2 columns.")      
		}
	} else {
		edges.sel1.infos <- edges.sel0.infos
	}

	### select target edges.part.infos and vertices.part.infos by exprs changes of genes
	## As the process in selecting mode & action.effect using 'edges' as reference, So here, only selecting 'edges' is enough
	# check if exprs change selection is valid
	if (sum(sel.exprs.change %in% kexprs.change) != length(sel.exprs.change)) {
		stop("Given parameter `sel.exprs.change` has some undefined values: ", 
			paste0(setdiff(sel.exprs.change, kexprs.change), collapse = ", "), "!")
	}
	## select by expression changes
	# cluster belongs inds
	tmp.inds.cluster.va <- which(vertices.all.infos$ClusterName == afterV.A.clustername)
	tmp.inds.cluster.vb <- which(vertices.all.infos$ClusterName == afterV.B.clustername) 
	# exprs change belongs inds
	tmp.inds.all.upreg <- which(vertices.all.infos$LogFC > 0)
	tmp.inds.all.dnreg <- which(vertices.all.infos$LogFC <= 0)
	# saved exprs change result
	inds.exprs.va <- integer()
	inds.exprs.vb <- integer()
	# up.up
	if ("Xup.Yup" %in% sel.exprs.change) {
		inds.exprs.va <- c(inds.exprs.va, intersect(tmp.inds.cluster.va, tmp.inds.all.upreg))
		inds.exprs.vb <- c(inds.exprs.vb, intersect(tmp.inds.cluster.vb, tmp.inds.all.upreg))
	}
	# up.dn
	if ("Xup.Ydn" %in% sel.exprs.change) {
		inds.exprs.va <- c(inds.exprs.va, intersect(tmp.inds.cluster.va, tmp.inds.all.upreg))
		inds.exprs.vb <- c(inds.exprs.vb, intersect(tmp.inds.cluster.vb, tmp.inds.all.dnreg))
	}
	# dn.up
	if ("Xdn.Yup" %in% sel.exprs.change) {
		inds.exprs.va <- c(inds.exprs.va, intersect(tmp.inds.cluster.va, tmp.inds.all.dnreg))
		inds.exprs.vb <- c(inds.exprs.vb, intersect(tmp.inds.cluster.vb, tmp.inds.all.upreg))
	}
	# dn.dn
	if ("Xdn.Ydn" %in% sel.exprs.change) {
		inds.exprs.va <- c(inds.exprs.va, intersect(tmp.inds.cluster.va, tmp.inds.all.dnreg))
		inds.exprs.vb <- c(inds.exprs.vb, intersect(tmp.inds.cluster.vb, tmp.inds.all.dnreg))
	}
	tmp.exprs.valid.UIDs <- vertices.all.infos[unique(c(inds.exprs.va, inds.exprs.vb)), "UID"]
	edges.sel2.infos <- edges.sel1.infos[intersect(which(edges.sel1.infos$from %in% tmp.exprs.valid.UIDs), which(edges.sel1.infos$to %in% tmp.exprs.valid.UIDs)), ]

	## select by direction
	# cluster belongs inds
	tmp.inds.cluster.vA <- which(vertices.all.infos$ClusterName == afterV.A.clustername)
	tmp.inds.cluster.vB <- which(vertices.all.infos$ClusterName == afterV.B.clustername) 
	edges.sel3.infos <- edges.sel2.infos
	if (!is.null(direction.X.to.Y) && direction.X.to.Y == TRUE) {
		tmp.from.matches <- which(edges.sel3.infos$from %in% vertices.all.infos[tmp.inds.cluster.vA, "UID"])
		tmp.to.matches <- which(edges.sel3.infos$to %in% vertices.all.infos[tmp.inds.cluster.vB, "UID"])
		edges.sel3.infos <- edges.sel3.infos[intersect(tmp.from.matches, tmp.to.matches), ]
	}
	if (!is.null(direction.X.to.Y) && !direction.X.to.Y) {
		tmp.from.matches <- which(edges.sel3.infos$from %in% vertices.all.infos[tmp.inds.cluster.vB, "UID"])
		tmp.to.matches <- which(edges.sel3.infos$to %in% vertices.all.infos[tmp.inds.cluster.vA, "UID"])
		edges.sel3.infos <- edges.sel3.infos[intersect(tmp.from.matches, tmp.to.matches), ]
	}

	# as it plot either directed or undirected graphs, new definition of action effects are given as below
	# for "A---B",              given type: "undirected"  --- kpred.ext.action.effect[1]
	# for "A-->B" or "A<--B",   given type: "positive"  --- kpred.ext.action.effect[c(2,3)]
	# for "A--|B" or "A|--B",   given type: "negative"  --- kpred.ext.action.effect[c(4,5)]
	# for "A--oB" or "Ao--B",   given type: "unspecified" --- kpred.ext.action.effect[c(6,7)]
	#
	### select target edges.part.infos and vertices.part.infos by mode & action.effect
	## check if valid, sel.mode.val, sel.action.effect.val
	predefined.mode.list <- kpred.mode
	predefined.action.effect.list <- kpred.action.effect
	if ((sum(sel.mode.val %in% predefined.mode.list) == length(sel.mode.val) ||
		 is.null(sel.mode.val)) &&
		(sum(sel.action.effect.val %in% predefined.action.effect.list) == length(sel.action.effect.val) ||
		 is.null(sel.action.effect.val))) {
		inds.full.a1 <- seq_len(nrow(edges.sel3.infos))
		# --- mode ---
		if (is.null(sel.mode.val)) {
			inds.mode.sel <- inds.full.a1
		} else {
			inds.mode.sel <- which(edges.sel3.infos[, "mode"] %in% sel.mode.val)
		}
		# --- action.effect ---
		if (is.null(sel.action.effect.val)) {
			inds.actf.sel <- inds.full.a1
		} else {
			inds.actf.sel <- which(edges.sel3.infos[, "action.effect"] %in% sel.action.effect.val)
		}
		# --- merge option ---
		if (sel.mode.action.option == "union") {
			inds.mode.actf.sel <- union(inds.mode.sel, inds.actf.sel)
		} else {
			if (sel.mode.action.option == "intersect") {
				inds.mode.actf.sel <- intersect(inds.mode.sel, inds.actf.sel)
			} else {
				stop("Mode-ActionEffect merge option error: undefined merge options! only 'union' and 'intersect' are allowed.")
			}
		}
		edges.part.infos <- edges.sel3.infos[inds.mode.actf.sel, ]

		# recheck if nrow() > 0
		if (nrow(edges.part.infos) == 0) {
			tmp.warn.names <- getOrigClusterNameTgVEInfo(object)
			stop("No given subset of interactions between cluster: ", tmp.warn.names$cluster.name.A, " and cluster: ", tmp.warn.names$cluster.name.B, "!")
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
		edges.infos = as.data.frame(edges.part.infos), 
		vertices.infos = as.data.frame(vertices.part.infos),
		vertices.apx.type.A = vertices.A.apx.types,
		vertices.apx.type.B = vertices.B.apx.types
		)
	object <- setTgVEInfo(object, VEinfos)
	return(object)
}


