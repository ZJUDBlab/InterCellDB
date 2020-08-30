

# [inside usage] [TODO-after-know math]
# This function evaluates the power of gene-interact pairs by the formula:
#   SUM(abs(LogFC[geneA]) * abs(LogFC[geneB])) / COUNT(<pairs>)
Inside.EvaluateByABS <- function(
  pairs.ref,
  colname.eval = c("inter.LogFC.A", "inter.LogFC.B")
) {
  eval.res <- 0.0
  if (nrow(pairs.ref) > 0) {
    for (i in 1:nrow(pairs.ref)) {
      eval.res <- eval.res + abs(pairs.ref[i, colname.eval[1]]) * abs(pairs.ref[i, colname.eval[2]])
    }
    eval.res <- eval.res / nrow(pairs.ref)
  }
  eval.res
}
EvaluateByFunc <- Inside.EvaluateByABS




# [inside usage]
# This function analyzes count and power of interaction pairs among all clusters.
#
# @param markers.remapped.all SEE @param fgenes.remapped.all in AnalyzeClustersInteracts().
# @param pairs.ref SEE @param pairs.ref in AnalyzeClustersInteracts().
# @param anno.location.ref SEE @param anno.location.ref in AnalyzeClustersInteracts().
# @param anno.type.ref SEE @param anno.type.ref in AnalyzeClustersInteracts().
# @param subgroup.options Settings that can select a subset of interaction pairs for analysis.
# @param user.type.database SEE @param user.type.database in AnalyzeClustersInteracts().
# @param sub.sel.user.type.colname SEE @param sub.sel.user.type.colname in AnalyzeClustersInteracts().
# @param restricted.some.genes SEE @param restricted.some.genes in AnalyzeClustersInteracts().
# @param restricted.gene.pairs SEE @param restricted.gene.pairs in AnalyzeClustersInteracts().
# @param ind.colname.end.dual SEE @param ind.colname.end.dual in AnalyzeClustersInteracts().
#
#
# import dplyr::left_join
# import dplyr::bind_rows
#
Inside.AnalyzeClustersInteracts <- function(
  markers.remapped.all,
  pairs.ref,
  anno.location.ref,
  anno.type.ref,
  subgroup.options,
  user.type.database = NULL,
  sub.sel.user.type.colname = NULL,
  restricted.some.genes = NULL,
  restricted.gene.pairs = NULL,
  ind.colname.end.dual = 4
) {
  # check if the column named "cluster" exists, so as "gene" and "avg_logFC"
  tmp.precheck.cols <- c("cluster", "gene", "avg_logFC")
  tmp.ind.valid.cols <- which(tmp.precheck.cols %in% colnames(markers.remapped.all))
  if (length(tmp.ind.valid.cols) != length(tmp.precheck.cols)) {
    stop("columns: ", paste0(tmp.precheck.cols[-tmp.ind.valid.cols], collapse = ", "), " are NOT AVAILABLE.")
  }
  # find all cluster, and factor() it
  fac.clusters <- levels(as.factor(markers.remapped.all[, "cluster"]))
  # check cluster name validity
  inds.ifvalid <- grep(kClustersSplit, fac.clusters, fixed = TRUE)
  if (length(inds.ifvalid) > 0) {
    stop("Name of clusters contains the specified `kClustersSplit`: \"", kClustersSplit, 
      "\", please change this global value by hand (like kClustersSplit <- '?') to avoid this error.")
  }
  # find all interactions between every two clusters
  interact.pairs.all <- list(
    list.clusters = fac.clusters,
    data.allpairs = list(),
    anno.allpairs = list(location.A = list(), location.B = list(),
                         type.A = list(), type.B = list()),
    name.allpairs = character(),
    cnt.allpairs = integer(),
    strength.allpairs = single()
  )
  prog.bar.p.p <- progress::progress_bar$new(total = length(fac.clusters))
  prog.bar.p.p$tick(0)
  for (i in 1:length(fac.clusters)) {
    for (k in 1:length(fac.clusters)) {
      interact.name <- paste0(fac.clusters[i], kClustersSplit, fac.clusters[k])  # naming interaction
      markers.rall.i <- markers.remapped.all[which(markers.remapped.all[, "cluster"] == fac.clusters[i]), ]  # extract i & k markers.remapped.* list
      markers.rall.k <- markers.remapped.all[which(markers.remapped.all[, "cluster"] == fac.clusters[k]), ]
      genes.list.i <- markers.rall.i[, "gene"]
      genes.list.k <- markers.rall.k[, "gene"]
      ## find pairs
      #1 A-B
      inds.ia.match <- which(pairs.ref[, "inter.GeneName.A"] %in% genes.list.i)
      inds.kb.match <- which(pairs.ref[, "inter.GeneName.B"] %in% genes.list.k)
      inds.ia.kb <- intersect(inds.ia.match, inds.kb.match)
      pairs.ia.kb <- pairs.ref[inds.ia.kb, ]
      #2 B-A
      inds.ib.match <- which(pairs.ref[, "inter.GeneName.B"] %in% genes.list.i)
      inds.ka.match <- which(pairs.ref[, "inter.GeneName.A"] %in% genes.list.k)
      inds.ib.ka <- intersect(inds.ib.match, inds.ka.match)
      pairs.ib.ka <- pairs.ref[inds.ib.ka, c(ReverseOddEvenCols(ind.colname.end.dual), (ind.colname.end.dual+1):ncol(pairs.ref))]  # reverse A-B to be matching B-A
      colnames(pairs.ib.ka) <- colnames(pairs.ref)
      # all
      pairs.sp.ikab <- rbind(pairs.ia.kb, pairs.ib.ka)
      pairs.sp.ikab <- DoPartUnique(pairs.sp.ikab, c(1,2))  # do part unique
      if (!is.null(restricted.gene.pairs)) {  # using restricted gene pairs
        if (is.data.frame(restricted.gene.pairs)) {  # repack pairs
          # check validity if given data.frame
          if (length(restricted.gene.pairs) != 2 || nrow(restricted.gene.pairs) == 0) {
            stop("Given data.frame of gene pairs have no records in details. Please check given data again.")
          }
          # process with data.frame
          tmp.df.a <- as.character(restricted.gene.pairs[, 1])
          tmp.df.b <- as.character(restricted.gene.pairs[, 2])
          tmp.mat.ab <- matrix(data = c(tmp.df.a, tmp.df.b), nrow = 2, byrow = TRUE)
          restricted.gene.pairs <- as.character(tmp.mat.ab)
        }
        # check validity
        if (length(restricted.gene.pairs) < 1 && length(restricted.gene.pairs) %% 2 != 0) {
          stop("Given wrong list of gene pairs with non-even positive length!")
        }
        # priority of restricted.gene.pairs, overwrite restricted.some.genes
        restricted.some.genes <- restricted.gene.pairs  # workflow merging, use with care
      }
      if (!is.null(restricted.some.genes)) {  # using retricted genes
        df.ref.somegenes <- data.frame(part.genes = restricted.some.genes, indicator.A = 1, indicator.B = 1, stringsAsFactors = FALSE)
        pairs.sp.ikab <- left_join(pairs.sp.ikab, df.ref.somegenes[, c("part.genes", "indicator.A")], by = c("inter.GeneName.A" = "part.genes"))
        pairs.sp.ikab <- pairs.sp.ikab[which(!is.na(pairs.sp.ikab[, "indicator.A"])), ]
        pairs.sp.ikab <- left_join(pairs.sp.ikab, df.ref.somegenes[, c("part.genes", "indicator.B")], by = c("inter.GeneName.B" = "part.genes"))
        pairs.sp.ikab <- pairs.sp.ikab[which(!is.na(pairs.sp.ikab[, "indicator.B"])), ]
        inds.rm.some.genes <- which(colnames(pairs.sp.ikab) %in% c("indicator.A", "indicator.B"))
        pairs.sp.ikab <- pairs.sp.ikab[, -inds.rm.some.genes]
      }
      if (!is.null(restricted.gene.pairs)) {  # using restricted gene pairs, finish rest process
        tmp.rgp.a <- restricted.gene.pairs[seq(1, length(restricted.gene.pairs) - 1, 2)]
        tmp.rgp.b <- restricted.gene.pairs[seq(2, length(restricted.gene.pairs), 2)]
        tmp.rpg.paste.conv <- paste(tmp.rgp.a, tmp.rgp.b, sep = "}")  # special sep, to avoid the same letter appeared in gene names
        tmp.rpg.paste.rev <- paste(tmp.rgp.b, tmp.rgp.a, sep = "}")   # and gene pairs are of no direction, the reverse of given gene pairs are generated
        tmp.ikab.paste <- paste(pairs.sp.ikab[, "inter.GeneName.A"], pairs.sp.ikab[, "inter.GeneName.B"], sep = "}")
        inds.ikab.pmatch.conv <- which(tmp.ikab.paste %in% tmp.rpg.paste.conv)
        inds.ikab.pmatch.rev <- which(tmp.ikab.paste %in% tmp.rpg.paste.rev)
        inds.ikab.pmatch <- union(inds.ikab.pmatch.conv, inds.ikab.pmatch.rev)
        pairs.sp.ikab <- pairs.sp.ikab[inds.ikab.pmatch, ]
      }
      # assign Exprs,LogFC to each pair
      #pairs.sp.ikab <- left_join(pairs.sp.ikab, markers.rall.i[, c("gene", "exprs")], by = c("inter.GeneName.A" = "gene"))
      #colnames(pairs.sp.ikab)[ncol(pairs.sp.ikab)] <- "inter.Exprs.A"
      #pairs.sp.ikab <- left_join(pairs.sp.ikab, markers.rall.k[, c("gene", "exprs")], by = c("inter.GeneName.B" = "gene"))
      #colnames(pairs.sp.ikab)[ncol(pairs.sp.ikab)] <- "inter.Exprs.B"
      pairs.sp.ikab <- left_join(pairs.sp.ikab, markers.rall.i[, c("gene", "avg_logFC")], by = c("inter.GeneName.A" = "gene"))
      colnames(pairs.sp.ikab)[ncol(pairs.sp.ikab)] <- "inter.LogFC.A"
      pairs.sp.ikab <- left_join(pairs.sp.ikab, markers.rall.k[, c("gene", "avg_logFC")], by = c("inter.GeneName.B" = "gene"))
      colnames(pairs.sp.ikab)[ncol(pairs.sp.ikab)] <- "inter.LogFC.B"
      # rearrange the columns. Make dual-matched columns being in the front, and single-matched columns being in the back.
      pairs.sp.ikab <- pairs.sp.ikab[, c(1:ind.colname.end.dual, (ncol(pairs.sp.ikab)-1):(ncol(pairs.sp.ikab)), (ind.colname.end.dual+1):(ncol(pairs.sp.ikab)-2))]
      ## do subgroup based on @param subgroup.options
      #1 - exprs change (based on logFC)
      pairs.subg.logfc <- list(xup.yup = data.frame(), xup.ydn = data.frame(), xdn.yup = data.frame(), xdn.ydn = data.frame())
      ind.A.up <- which(pairs.sp.ikab[, "inter.LogFC.A"] > 0)
      ind.A.dn <- which(pairs.sp.ikab[, "inter.LogFC.A"] < 0)
      ind.B.up <- which(pairs.sp.ikab[, "inter.LogFC.B"] > 0)
      ind.B.dn <- which(pairs.sp.ikab[, "inter.LogFC.B"] < 0)
      if ("Xup.Yup" %in% subgroup.options$exprs.logfc) {
        pairs.subg.logfc$xup.yup <- pairs.sp.ikab[intersect(ind.A.up, ind.B.up), ]
      }
      if ("Xup.Ydown" %in% subgroup.options$exprs.logfc) {
        pairs.subg.logfc$xup.ydn <- pairs.sp.ikab[intersect(ind.A.up, ind.B.dn), ]
      }
      if ("Xdown.Yup" %in% subgroup.options$exprs.logfc) {
        pairs.subg.logfc$xdn.yup <- pairs.sp.ikab[intersect(ind.A.dn, ind.B.up), ]
      }
      if ("Xdown.Ydown" %in% subgroup.options$exprs.logfc) {
        pairs.subg.logfc$xdn.ydn <- pairs.sp.ikab[intersect(ind.A.dn, ind.B.dn), ]
      }
      pairs.v1.after.logfc <- rbind(rbind(pairs.subg.logfc$xup.yup, pairs.subg.logfc$xup.ydn), rbind(pairs.subg.logfc$xdn.yup, pairs.subg.logfc$xdn.ydn))
      ## 2 - Location
      # Location - A
      if (is.null(subgroup.options[["X.Location"]])) {
        this.A.locations <- unique(pairs.v1.after.logfc[, c("inter.GeneID.A", "inter.GeneName.A")])
        colnames(this.A.locations) <- c("GeneID", "Gene.name")
        this.A.locations <- left_join(this.A.locations, anno.location.ref[, c("GeneID", "Gene.name", "GO.Term.target", "Source", "Evidence", "score")], by = c("GeneID", "Gene.name"))
        pairs.v2.pre.A <- pairs.v1.after.logfc
      } else {
        tmp.inds.A.locations <- which(anno.location.ref$GeneID %in% levels(factor(pairs.v1.after.logfc[, "inter.GeneID.A"])))  # fetch location
        tmp.A.locations <- anno.location.ref[tmp.inds.A.locations, c("GeneID", "Gene.name", "GO.Term.target", "Source", "Evidence", "score")]
        this.A.locations <- tmp.A.locations[which(tmp.A.locations[, "GO.Term.target"] %in% subgroup.options[["X.Location"]]), ]
        # 2.1 slim
        pairs.v2.pre.A <- pairs.v1.after.logfc[which(pairs.v1.after.logfc[, "inter.GeneID.A"] %in% this.A.locations[, "GeneID"]), ]
      }
      # Location - B
      if (is.null(subgroup.options[["Y.Location"]])) {
        this.B.locations <- unique(pairs.v2.pre.A[, c("inter.GeneID.B", "inter.GeneName.B")])
        colnames(this.B.locations) <- c("GeneID", "Gene.name")
        this.B.locations <- left_join(this.B.locations, anno.location.ref[, c("GeneID", "Gene.name", "GO.Term.target", "Source", "Evidence", "score")], by = c("GeneID", "Gene.name"))
        pairs.v2.aft.loc <- pairs.v2.pre.A
      } else {
        tmp.inds.B.locations <- which(anno.location.ref$GeneID %in% levels(factor(pairs.v2.pre.A[, "inter.GeneID.B"])))  # fetch location
        tmp.B.locations <- anno.location.ref[tmp.inds.B.locations, c("GeneID", "Gene.name", "GO.Term.target", "Source", "Evidence", "score")]
        this.B.locations <- tmp.B.locations[which(tmp.B.locations[, "GO.Term.target"] %in% subgroup.options[["Y.Location"]]), ]
        # 2.2 slim
        pairs.v2.aft.loc <- pairs.v2.pre.A[which(pairs.v2.pre.A[, "inter.GeneID.B"] %in% this.B.locations[, "GeneID"]), ]
      }
      ## 3 - Type
      # Type - A
      if (is.null(subgroup.options[["X.Type"]])) {
        this.A.types <- unique(pairs.v2.aft.loc[, c("inter.GeneID.A", "inter.GeneName.A")])
        colnames(this.A.types) <- c("GeneID", "Gene.name")
        this.A.types <- left_join(this.A.types, anno.type.ref[, c("GeneID", "Gene.name", "Keyword.Name")], by = c("GeneID", "Gene.name"))
        pairs.v3.spre.A <- pairs.v2.aft.loc
      } else {
        tmp.inds.A.types <- which(anno.type.ref$GeneID %in% levels(factor(pairs.v2.aft.loc[, "inter.GeneID.A"])))
        tmp.A.types <- anno.type.ref[tmp.inds.A.types, c("GeneID", "Gene.name", "Keyword.Name")]
        this.A.types <- tmp.A.types[which(tmp.A.types[, "Keyword.Name"] %in% subgroup.options$X.Type), ]
        # 3.1 slim
        pairs.v3.spre.A <- pairs.v2.aft.loc[which(pairs.v2.aft.loc[, "inter.GeneID.A"] %in% this.A.types[, "GeneID"]), ]
      }
      # Type - B
      if (is.null(subgroup.options[["Y.Type"]])) {
        this.B.types <- unique(pairs.v3.spre.A[, c("inter.GeneID.B", "inter.GeneName.B")])
        colnames(this.B.types) <- c("GeneID", "Gene.name")
        this.B.types <- left_join(this.B.types, anno.type.ref[, c("GeneID", "Gene.name", "Keyword.Name")], by = c("GeneID", "Gene.name"))
        pairs.v3.aft.type <- pairs.v3.spre.A
      } else {
        tmp.inds.B.types <- which(anno.type.ref$GeneID %in% levels(factor(pairs.v3.spre.A[, "inter.GeneID.B"])))
        tmp.B.types <- anno.type.ref[tmp.inds.B.types, c("GeneID", "Gene.name", "Keyword.Name")]
        this.B.types <- tmp.B.types[which(tmp.B.types[, "Keyword.Name"] %in% subgroup.options$Y.Type), ]
        # 3.2 slim
        pairs.v3.aft.type <- pairs.v3.spre.A[which(pairs.v3.spre.A[, "inter.GeneID.B"] %in% this.B.types[, "GeneID"]), ]
      }
      ## 4 - user.type
      if (!is.null(user.type.database) && !is.null(sub.sel.user.type.colname)) {
        ## as user.type.database usually doesn't cover all genes, so here use different strategy
        # user.type - A
        if (is.null(subgroup.options[["X.user.type"]])) {
          pairs.v4.upre.A <- pairs.v3.aft.type
        } else {
          # slim the dataset
          tmp.inds.A.user.types <- which(user.type.database$GeneID %in% levels(factor(pairs.v3.aft.type[, "inter.GeneID.A"])))
          tmp.A.user.types <- user.type.database[tmp.inds.A.user.types, c("GeneID", "Gene.name", sub.sel.user.type.colname)]
          this.A.user.types <- tmp.A.user.types[which(tmp.A.user.types[, sub.sel.user.type.colname] %in% subgroup.options$X.user.type), ]
          pairs.v4.upre.A <- pairs.v3.aft.type[which(pairs.v3.aft.type[, "inter.GeneID.A"] %in% this.A.user.types[, "GeneID"]), ]
        }
        # user.type - B
        if (is.null(subgroup.options[["Y.user.type"]])) {
          pairs.v4.aft.user.type <- pairs.v4.upre.A
        } else {
          # slim the dataset
          tmp.inds.B.user.types <- which(user.type.database$GeneID %in% levels(factor(pairs.v4.upre.A[, "inter.GeneID.B"])))
          tmp.B.user.types <- user.type.database[tmp.inds.B.user.types, c("GeneID", "Gene.name", sub.sel.user.type.colname)]
          this.B.user.types <- tmp.B.user.types[which(tmp.B.user.types[, sub.sel.user.type.colname] %in% subgroup.options$Y.user.type), ]
          pairs.v4.aft.user.type <- pairs.v4.upre.A[which(pairs.v4.upre.A[, "inter.GeneID.B"] %in% this.B.user.types[, "GeneID"]), ]  
        }
      } else {
        pairs.v4.aft.user.type <- pairs.v3.aft.type
      }
      ## after doing subgroup
      pairs.subg.result <- pairs.v4.aft.user.type
      ## re-slim with Location
      func.location.score.inside <- function(data.loc, option) {
        ret.val <- data.loc
        if (is.character(option)) {  # use pre-defined strategies
          if (option == "the most confident") {  # the only strategy supported yet
            tmp.gene.ids <- levels(factor(data.loc[, "GeneID"]))
            tmp.ret.dfs <- lapply(tmp.gene.ids, ref.data = data.loc, function(x, ref.data) {
                this.gene.loc <- ref.data[which(ref.data[, "GeneID"] == x), ]
                this.max.score <- max(this.gene.loc[, "score"])
                this.gene.loc[which(this.gene.loc[, "score"] == this.max.score), ]
              })
            ret.val <- bind_rows(tmp.ret.dfs)
            if (nrow(ret.val) == 0) {
              ret.val <- data.loc[1, ][-1, ]  # to save the colnames
            }
          }  # else, do nothing
        } else {  # use score
          if (is.numeric(option)) {
            option <- as.integer(option)  # coerce to integer score
            ret.val <- data.loc[which(data.loc[, "score"] %in% option), ]
          } else {
            warning("Score limit get neither NUMERIC nor PRE-DEFINED options, and will not do anything.")
          }
        }
        ret.val
      }
      # - Location
      this.A.locations <- this.A.locations[which(this.A.locations[, "GeneID"] %in% pairs.subg.result[, "inter.GeneID.A"]), ]
      this.A.locations <- func.location.score.inside(this.A.locations, subgroup.options$X.Location.score.limit)
      this.B.locations <- this.B.locations[which(this.B.locations[, "GeneID"] %in% pairs.subg.result[, "inter.GeneID.B"]), ]
      this.B.locations <- func.location.score.inside(this.B.locations, subgroup.options$Y.Location.score.limit)
      # re-slim rescue here
      pairs.subg.result <- pairs.subg.result[which(pairs.subg.result[, "inter.GeneID.A"] %in% this.A.locations[, "GeneID"]), ]
      pairs.subg.result <- pairs.subg.result[which(pairs.subg.result[, "inter.GeneID.B"] %in% this.B.locations[, "GeneID"]), ]
      # - Type
      this.A.types <- this.A.types[which(this.A.types[, "GeneID"] %in% pairs.subg.result[, "inter.GeneID.A"]), ]
      this.B.types <- this.B.types[which(this.B.types[, "GeneID"] %in% pairs.subg.result[, "inter.GeneID.B"]), ]
      ## add to result
      interact.pairs.all$data.allpairs[[interact.name]] <- pairs.subg.result
      interact.pairs.all$anno.allpairs$location.A[[interact.name]] <- this.A.locations
      interact.pairs.all$anno.allpairs$location.B[[interact.name]] <- this.B.locations
      interact.pairs.all$anno.allpairs$type.A[[interact.name]] <- this.A.types
      interact.pairs.all$anno.allpairs$type.B[[interact.name]] <- this.B.types     
      interact.pairs.all$name.allpairs <- append(interact.pairs.all$name.allpairs, interact.name)
      interact.pairs.all$cnt.allpairs <- append(interact.pairs.all$cnt.allpairs, nrow(pairs.subg.result))
      interact.pairs.all$strength.allpairs <- append(interact.pairs.all$strength.allpairs,
        EvaluateByFunc(pairs.subg.result, c("inter.LogFC.A", "inter.LogFC.B")))
    }
    prog.bar.p.p$tick()
  }
  #end# return
  interact.pairs.all
}




#' Analyze interaction pairs based on clusters
#'
#' @description
#' This function analyzes count and power of interaction pairs among all clusters.
#'
#' @inheritParams Inside.DummyFgenes 
#' @inheritParams Inside.DummyPairsRef
#' @inheritParams Inside.DummyAnnoLocationRefDB  
#' @inheritParams Inside.DummyAnnoTypeRefDB
#' @param user.type.database Data.frame. It specifies one special data.frame that contains user-defined informations.
#' For generating such database, see \code{Tool.AddUserRestrictDB} for help. 
#' @param restricted.some.genes Character. Analysis will be restricted in interaction pairs that contain
#' at least one genes given in this parameter.
#' @param restricted.gene.pairs Character or data.frame. Analysis will be restricted in given gene pairs(i.e. interaction pairs). 
#' The given format is explained in the below, see details for help.
#' @param sub.sel.exprs.changes Character. Use subset of \code{c("Xup.Yup", "Xup.Ydown", "Xdown.Yup", "Xdown.Ydown")}.
#' @param sub.sel.X.Location Character. Its value depends on the database used, see details for help.
#' @param sub.sel.X.Location.score.limit Character or Integer. The one in \code{character()} will be treated as predefined strategy, while
#' the one in \code{integer()} will be treated as score limit range. See details for help.
#' @param sub.sel.Y.Location Like \code{sub.sel.X.Location}.
#' @param sub.sel.Y.Location.score.limit Like \code{sub.sel.X.Location.score.limit}.
#' @param sub.sel.X.Type Character. Its value depends on the database used, see details for help.
#' @param sub.sel.Y.Type Like \code{sub.sel.X.Type}.
#' @param sub.sel.user.type.colname Character. It gives the column name that will be used for user-defined purposes.
#' @param sub.sel.X.user.type Its mode depends on what datatype user has defined. The function will accept anything given in this parameter
#' without any additional check.
#' @param sub.sel.Y.user.type Like \code{sub.sel.X.user.type}.
#' @param ind.colname.end.dual Integer. Use default value provided only when the pairs.ref database is modified by users.
#'
#' @details
#' This function gives as much as possible freedom for users to choose to analyze a subset of interaction pairs.
#' \itemize{
#'   \item \code{restricted.some.genes}: (1) Use GO terms or IDs to get some genes.
#'         (2) Directly specify some genes(use remap.ref.<taxonomy>.dblist to get authorized genenames).
#'   \item \code{restricted.gene.pairs}: If it is given, the parameter \code{restricted.some.genes} will be ignored. For this parameter, the allowed
#'         data formats are:
#'
#'         (1) Character vector with genes in every pair listed in consecutive order, e.g. if \code{c("geneA", "geneB", "geneC", "geneD")} is given, 
#'         it will be interpreted to be 2 gene pairs: geneA - geneB, geneC - geneD. In this format, the length of the given vector is implied to be 
#'         even number.
#'
#'         (2) Data.frame with 2 columns, which records the gene pairs in every row.
#'
#'   \item \code{sub.sel.exprs.changes}: For interaction pair X-Y, take "Xup.Ydown" for example, and it means the subset that
#'         gene X is up-regulated and gene Y is down-regulated is selected.
#'   \item \code{sub.sel.X.Location} & \code{sub.sel.Y.Location}: Location of genes, i.e. Extracellular Region, Plasma Membrane, Cytosol, etc. 
#'         To fetch available values from current used database, \bold{e.g.}
#'
#'         In human, use database \code{anno.location.human.ref.db}, and use the code \code{levels(factor(anno.location.human.ref.db$GO.Term.target))} to 
#'         see all available values.        
#'   \item \code{sub.sel.X.Location.score.limit} & \code{sub.sel.Y.Location.score.limit}: Further restrition applied on location of genes, which can be given 
#'         in \code{character()} or \code{integer()}.
#'
#'         In \code{character()}, this gives the predefined strategies, which support only 1 strategy in current setting. The "the most confident" strategy, 
#'         which selects all locations of the highest score for one specific gene, and it means multiply locations can be selected together if their scores are 
#'         the same.
#'
#'         In \code{integer()}, this gives the range of score. The full range of score is 1~5, and other integers out of this range will not be recognized, 
#'         \bold{e.g.} \code{sub.sel.X.Location.score.limit = c(4,5)} is fine, and \code{sub.sel.X.Location.score.limit = c(4,5,6)} has the same meaning as the former 
#'         one when 6 is out of range and will not be used.
#'         \bold{TO NOTE}: BOTH strategies will have some effect on final runout data, e.g if use \code{.score.limit = c(5)}, and get gene A with Plasma Membrane in score 4, 
#'         then gene A will not be seen as it can located in Plasma Membrane.  
#'   \item \code{sub.sel.X.Type} & \code{sub.sel.Y.Type}: More specific type definition of genes, i.e. Cytokine, G-protein Coupled Receptor, etc.
#          To fetch available values, the method is like the above one but uses the database for types, \bold{e.g.} 
#'
#'         In mouse, use database \code{anno.type.mouse.ref.db}, and use the code \code{levels(factor(anno.type.mouse.ref.db$Keyword.Name))}.
#'   \item \code{sub.sel.X.user.type} & \code{sub.sel.Y.user.type} + \code{sub.sel.user.type.colname}: Those give a more flexible way for subset selection.
#'         It uses the database given in another parameter \code{user.type.database}. \code{sub.sel.user.type.colname} selects which column to be used. 
#'         \code{sub.sel.X.user.type} & \code{sub.sel.Y.user.type} select the specific items in that column, and if any of them were set \code{NULL}, the  
#'         corresponding part will not add restriction on this user type.
#'   \item \code{ind.colname.end.dual}: It is the index when column of pairs.ref begin to provide information about the pair itself
#'          but not for the 2 genes in that pair. In almost all circumstances, this parameter is required not to be modified.
#' }
#'
#' @return A list.
#' \itemize{
#'   \item list.clusters: names of clusters.
#'   \item data.allpairs: a list of all interaction pairs. e.g. cluster: Astrocyte -> cluster: Microglia, if kClustersSplit is "%",
#'         then, the interaction pair will be named "Astrocyte%Microglia", and use $data.allpairs[["Astrocyte%Microglia"]] to get details.
#'   \item anno.allpairs: a list of lists. The sublists are 
#'         \itemize{
#'           \item location.A: it records locations of A in gene pairs formatted as A-B.
#'           \item location.B: it records locations of B in gene pairs formatted as A-B.
#'           \item type.A: it records molecular functions of A in gene pairs formatted as A-B.
#'           \item type.B: it records molecular functions of B in gene pairs formatted as A-B.
#'         }
#'   \item name.allpairs: the names of $data.allpairs
#'   \item cnt.allpairs: count of each interaction pair in $data.allpairs
#'   \item strength.allpairs: power of each interaction pair in $data.allpairs
#' }
#'
#' @examples
#' \dontrun{
#' # restrict to some subset
#' tmp.gene.pairs <- data.frame(geneA = c("Fas", "Tgfbr2", "Rac1"),
#'                              geneB = c("Casp8", "Tgfbr1", "Ctnnb1"))
#' AnalyzeClustersInteracts(fgenes.remapped.all,
#'   pairs.ref = pairs.mouse.experiments.db,
#'   restricted.gene.pairs = tmp.gene.pairs,
#'   sub.sel.X.Location = c("Plasma Membrane"),
#'   sub.sel.Y.Location = c("Plasma Membrane", "Cytosol")
#' )
#' }
#'
#'
#'
#' @importFrom dplyr bind_rows left_join
#' @import progress
#'
#' @export
#'
AnalyzeClustersInteracts <- function(
  fgenes.remapped.all,
  pairs.ref,
  anno.location.ref,
  anno.type.ref,
  user.type.database = NULL,
  restricted.some.genes = NULL,
  restricted.gene.pairs = NULL,
  sub.sel.exprs.changes = c("Xup.Yup", "Xup.Ydown", "Xdown.Yup", "Xdown.Ydown"),
  sub.sel.X.Location = character(),
  sub.sel.X.Location.score.limit = c("the most confident"),
  sub.sel.Y.Location = character(),
  sub.sel.Y.Location.score.limit = c("the most confident"),
  sub.sel.X.Type = character(),
  sub.sel.Y.Type = character(),
  sub.sel.user.type.colname = NULL,
  sub.sel.X.user.type = NULL,
  sub.sel.Y.user.type = NULL,
  ind.colname.end.dual = 4
) {
  # generate default settings
  user.settings <- list(
    exprs.logfc = c("Xup.Yup", "Xup.Ydown", "Xdown.Yup", "Xdown.Ydown"),
    X.Location.score.limit = sub.sel.X.Location.score.limit,
    Y.Location.score.limit = sub.sel.Y.Location.score.limit)
  ### then the user parameters
  ## fold change
  if (!is.null(sub.sel.exprs.changes) && length(sub.sel.exprs.changes) != 0) {
    # rather matching
    rather.m.exprs.changes <- character()
    for (i.item in sub.sel.exprs.changes) {
      tmp.sl.i <- strsplit(i.item, split = ".", fixed = TRUE)[[1]]
      if (length(tmp.sl.i) >= 2) {  # rather matching need (>= 2) splits. So, format "???up.???down" will be correctly recognized.
        tif.Xup   <- if (length(grep("up", tolower(tmp.sl.i[1]), fixed = TRUE))   > 0) TRUE else FALSE
        tif.Xdown <- if (length(grep("down", tolower(tmp.sl.i[1]), fixed = TRUE)) > 0) TRUE else FALSE
        tif.Yup   <- if (length(grep("up", tolower(tmp.sl.i[2]), fixed = TRUE))   > 0) TRUE else FALSE
        tif.Ydown <- if (length(grep("down", tolower(tmp.sl.i[2]), fixed = TRUE)) > 0) TRUE else FALSE
        if (tif.Xup && tif.Yup)   rather.m.exprs.changes <- append(rather.m.exprs.changes, "Xup.Yup")
        if (tif.Xup && tif.Ydown)   rather.m.exprs.changes <- append(rather.m.exprs.changes, "Xup.Ydown")
        if (tif.Xdown && tif.Yup) rather.m.exprs.changes <- append(rather.m.exprs.changes, "Xdown.Yup")
        if (tif.Xdown && tif.Ydown) rather.m.exprs.changes <- append(rather.m.exprs.changes, "Xdown.Ydown")
      }
    }
    user.settings$exprs.logfc <- rather.m.exprs.changes
  }
  ## Location
  # for X
  if (!is.null(sub.sel.X.Location) && length(sub.sel.X.Location) != 0) {
    user.settings[["X.Location"]] <- as.character(Tc.Cap.simple.vec(sub.sel.X.Location))
  }
  # for Y
  if (!is.null(sub.sel.Y.Location) && length(sub.sel.Y.Location) != 0) {
    user.settings[["Y.Location"]] <- as.character(Tc.Cap.simple.vec(sub.sel.Y.Location))
  }
  ## Type
  # for X
  if (!is.null(sub.sel.X.Type) && length(sub.sel.X.Type) != 0) {
    user.settings[["X.Type"]] <- as.character(Tc.Cap.simple.vec(sub.sel.X.Type))
  }
  # for Y
  if (!is.null(sub.sel.Y.Type) && length(sub.sel.Y.Type) != 0) {
    user.settings[["Y.Type"]] <- as.character(Tc.Cap.simple.vec(sub.sel.Y.Type))
  }
  ## user type
  # check
  if (!is.null(user.type.database) && 
      (class(user.type.database) == "data.frame" && nrow(user.type.database) > 0) &&
      (length(sub.sel.user.type.colname) == 1 && (sub.sel.user.type.colname %in% colnames(user.type.database) == TRUE))) {
    # authorized check
    if (sum(c("GeneID", "Gene.name") %in% colnames(user.type.database)) != 2) {
      stop("User type error: not a authorized database generated by routine process defined in this package.")
    }
    # for X
    if (!is.null(sub.sel.X.user.type) && length(sub.sel.X.user.type) != 0) {
      user.settings[["X.user.type"]] <- sub.sel.X.user.type
    }
    # for Y
    if (!is.null(sub.sel.Y.user.type) && length(sub.sel.Y.user.type) != 0) {
      user.settings[["Y.user.type"]] <- sub.sel.Y.user.type
    }
  } else {
    tmp.check.db.exist <- !is.null(user.type.database) && (class(user.type.database) == "data.frame" && nrow(user.type.database) > 0)
    tmp.check.colname.valid <- (tmp.check.db.exist == TRUE) && (length(sub.sel.user.type.colname) == 1 && (sub.sel.user.type.colname %in% colnames(user.type.database) == TRUE))
    tmp.check.selection.X <- !is.null(sub.sel.X.user.type) && (length(sub.sel.X.user.type) != 0)
    tmp.check.selection.Y <- !is.null(sub.sel.Y.user.type) && (length(sub.sel.Y.user.type) != 0)
    tmp.check.sel <- tmp.check.selection.X || tmp.check.selection.Y
    if (tmp.check.sel && !tmp.check.colname.valid) {
      if (tmp.check.db.exist) {
        stop("User type error: given undefined column name.")
      } else {
        stop("User type error: invalid user-defined database.")
      }
    }
  }
  # before running, print the selection used
  cat(paste0("\n---Strategy Used---", 
    "\nexprs.change: ", paste0(user.settings$exprs.logfc, collapse = ", "), 
    "\nLocation in X: ", ifelse(is.null(user.settings[["X.Location"]]), "-", paste0(user.settings[["X.Location"]], collapse = ", ")),  
    "\nLocation score in X: ", ifelse(is.null(user.settings[["X.Location.score.limit"]]), "-", paste0(user.settings[["X.Location.score.limit"]], collapse = ", ")), 
    "\nLocation in Y: ", ifelse(is.null(user.settings[["Y.Location"]]), "-", paste0(user.settings[["Y.Location"]], collapse = ", ")), 
    "\nLocation score in Y: ", ifelse(is.null(user.settings[["Y.Location.score.limit"]]), "-", paste0(user.settings[["Y.Location.score.limit"]], collapse = ", ")), 
    "\nType in X: ", ifelse(is.null(user.settings[["X.Type"]]), "-", paste0(user.settings[["X.Type"]], collapse = ", ")), 
    "\nType in Y: ", ifelse(is.null(user.settings[["Y.Type"]]), "-", paste0(user.settings[["Y.Type"]], collapse = ", ")), 
    "\nUse user database: ", ifelse(!is.null(user.type.database), "TRUE", "FALSE"), 
    "\nUse some genes: ", ifelse(!is.null(restricted.some.genes), "TRUE", "FALSE"), 
    "\nUse some gene pairs: ", ifelse(!is.null(restricted.gene.pairs), "TRUE", "FALSE"),
    "\n"
  ))
  ### then run the analysis
  res <- Inside.AnalyzeClustersInteracts(fgenes.remapped.all, pairs.ref, 
            anno.location.ref, anno.type.ref, user.settings, 
            user.type.database = user.type.database,
            sub.sel.user.type.colname = sub.sel.user.type.colname,
            restricted.some.genes = restricted.some.genes,
            restricted.gene.pairs = restricted.gene.pairs,
            ind.colname.end.dual = ind.colname.end.dual)
  #end# return value
  res
}





#' Get result of analysis in full view
#' 
#' @description
#' This function focuses on some subset of interaction pairs, and get rather fine grained
#' result from those.
#'
#' @inheritParams Inside.DummyInteractPairsActed 
#' @param show.clusters.in.x Vector. Clusters(use cluster names) that are chosen to show in x-axis by users.
#' @param show.clusters.in.y Vector. Clusters(use cluster names) that are chosen to show in y-axis by users.
#' @param power.max.limit Numeric. Specify the upper limit of power, whose value is highly user-defined and data-dependent.
#' @param power.min.limit Numeric. Specify the lower limit of power, like \code{power.max.limit}.
#' @param hide.power.label Logic. If TRUE, the label appended for power value will be hidden, otherwise, the label will be kept.
#' @param cnt.max.limit Numeric. Specify the upper limit of count, whose value is highly user-defined and data-dependent.
#' @param cnt.min.limit Numeric. Specify the lower limit of count, like \code{cnt.max.limit}.
#' @param hide.cnt.label Logic. If TRUE, the label appended for count value will be hidden, otherwise, the label will be kept.
#' @param nodes.colour.seq Character. Given colours will be used to generate colour gradient for plotting.
#' @param nodes.colour.value.seq Numeric. If set NULL, evenly place each colour in \code{nodes.colour.seq} vector.
#' Otherwise, numeric values with the same length of \code{nodes.colour.seq} should be given to specify the positions corresponding to each colour.
#' See parameter \code{values} in \code{ggplot2::scale_colour_gradientn} for details.
#' @param label.power.options List. Options are hjust, vjust, nudge.x, nudge.y, size. See \code{ggplot2::geom_text} for details.
#' \itemize{
#'   \item hjust:   horizontal justification. Use either a string in ("left", "center", "right"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item vjust:   vertical justification. Use either a string in ("top", "middle", "bottom"), or a number between 0 and 1. See
#'                  \code{vignette("ggplot2-specs", package = "ggplot2")} for details.
#'   \item nudge.x: horizontal adjustment to nudge labels by.
#'   \item nudge.y: vertical adjustment to nudge labels by.
#'   \item size:    label size.
#' }
#' @param label.cnt.options List. Options are like \code{label.power.options}.
#' @param plot.axis.x.name Character. X-axis name when plot is shown.
#' @param plot.axis.y.name Character. Y-axis name when plot is shown.
#' 
#' @details
#' This function uses some subset of interaction pairs, and does calculations in cluster level.
#' It calculates:
#' \itemize{
#'   \item count: count of interaction pairs.
#'   \item power: strength of interaction pairs, which is calculated by formula: [TODO-formula].
#' }
#'
#' @return List. Use \code{Tool.ShowPlot()} to see the \bold{plot}, \code{Tool.WriteTables()} to save the result \bold{table} in .csv files.
#' \itemize{
#'   \item plot: the object of \pkg{ggplot2}.
#'   \item table: a list of \code{data.frame}.
#' }
#'
#' @examples
#' \dontrun{
#' # show interaction pairs in some clusters, and set power & cnt max.limit
#' GetResult.SummaryClustersInteracts(interact.pairs.acted, 
#'   show.clusters.in.x = c("Microglia", "Astrocyte"), 
#'   show.clusters.in.y = c("Endothelial cell"),
#'   power.max.limit = 1.6,
#'   cnt.max.limit = 50
#' )
#'   
#' # change the label
#' GetResult.SummaryClustersInteracts(interact.pairs.acted, 
#'   show.clusters.in.x = c("Astrocyte"), 
#'   show.clusters.in.y = NULL,  # show all clusters in y-axis
#'   power.max.limit = 1.6,
#'   cnt.max.limit = 30,
#'   label.power.options = list(nudge.y = -0.2, size = 3),
#'   label.cnt.options = list(hjust = "left", vjust = "center", nudge.x = 0.4, nudge.y = 0, size = 2)
#' )
#' }
#'
#'
#'
#' @import ggplot2
#'
#' @export
#'
GetResult.SummaryClustersInteracts <- function(
  interact.pairs.acted,
  show.clusters.in.x = NULL,
  show.clusters.in.y = NULL,
  power.max.limit = NULL,
  power.min.limit = NULL,
  hide.power.label = FALSE,
  cnt.max.limit = NULL,
  cnt.min.limit = NULL,
  hide.cnt.label = FALSE,  # hide the labels 
  nodes.colour.seq = c("#00809D", "#EEEEEE", "#C30000"),
  nodes.colour.value.seq = c(0.0, 0.5, 1.0),
  label.power.options = list(hjust = "middle", vjust = "top", nudge.x = 0, nudge.y = -0.3, size = 2),
  label.cnt.options = list(hjust = "left", vjust = "center", nudge.x = 0.3, nudge.y = 0, size = 2),
  plot.axis.x.name = "clusters-x",
  plot.axis.y.name = "clusters-y"
) {
  # param user settings
  default.label.power.options <- list(hjust = "middle", vjust = "top", nudge.x = 0, nudge.y = -0.3, size = 2)
  default.label.cnt.options <- list(hjust = "left", vjust = "center", nudge.x = 0.3, nudge.y = 0, size = 2)
  user.label.power.options <- as.list(label.power.options)
  user.label.cnt.options <- as.list(label.cnt.options)
  # missing items in user settings
  missed.names.power.options <- setdiff(names(default.label.power.options), names(user.label.power.options))
  missed.names.cnt.options <- setdiff(names(default.label.cnt.options), names(user.label.cnt.options))
  # add missing items in user settings
  for (miss.p in missed.names.power.options) {
    user.label.power.options[[miss.p]] <- default.label.power.options[[miss.p]]
  }
  for (miss.cn in missed.names.cnt.options) {
    user.label.cnt.options[[miss.cn]] <- default.label.cnt.options[[miss.cn]]
  }
  # process
  fac.clusters <- interact.pairs.acted$list.clusters  # all clusters related
  fac.x.clusters <- fac.clusters
  fac.y.clusters <- fac.clusters
  if (!is.null(show.clusters.in.x) && length(show.clusters.in.x) != 0) {  # select part of clusters to be shown in x-axis
    ind.m.x <- match(show.clusters.in.x, fac.x.clusters)
    fac.x.clusters <- fac.x.clusters[ind.m.x[which(!is.na(ind.m.x))]]
    print(paste0("Reading clusters shown in X-axis: ", paste0(fac.x.clusters, collapse = ", "), "."))
  }
  if (!is.null(show.clusters.in.y) && length(show.clusters.in.y) != 0) {  # select part of clusters to be shown in y-axis
    ind.m.y <- match(show.clusters.in.y, fac.y.clusters)
    fac.y.clusters <- fac.y.clusters[ind.m.y[which(!is.na(ind.m.y))]]
    print(paste0("Reading clusters shown in Y-axis: ", paste0(fac.y.clusters, collapse = ", "), "."))
  }
  if (length(fac.x.clusters) < 1 || length(fac.y.clusters) < 1) {  # check
    stop("Error: X-Y plot needs at least one item in each axis!")
  }
  ### data process for selected part of data
  col.x.data <- rep(fac.x.clusters, each = length(fac.y.clusters))
  col.y.data <- rep(fac.y.clusters, times = length(fac.x.clusters))
  names.part.data <- character()
  for (i.ind in 1:length(col.x.data)) {
    names.part.data <- append(names.part.data, paste0(col.x.data[i.ind], kClustersSplit, col.y.data[i.ind]))
  }
  names.ref.data <- interact.pairs.acted$name.allpairs
  ind.names.m <- match(names.part.data, names.ref.data)
  # ------
  # for cnt & power, use limit functions
  Limit.Max.inside <- function(x, eval.max) {
    res <- if (x > eval.max) eval.max else x
    res
  }
  Limit.Min.inside <- function(x, eval.min) {
    res <- if (x < eval.min) eval.min else x
    res
  }
  # ------
  # cnt
  cnt.xy.data <- interact.pairs.acted$cnt.allpairs
  cnt.part.data <- cnt.xy.data[ind.names.m]
  # cnt limit
  cnt.limit.part.data <- cnt.part.data
  if (!is.null(cnt.max.limit)) {
    cnt.limit.part.data <- sapply(cnt.limit.part.data, eval.max = cnt.max.limit, Limit.Max.inside)
  }
  if (!is.null(cnt.min.limit)) {
    cnt.limit.part.data <- sapply(cnt.limit.part.data, eval.min = cnt.min.limit, Limit.Min.inside)
  }
  # power
  power.xy.data <- interact.pairs.acted$strength.allpairs
  power.part.data <- power.xy.data[ind.names.m]
  # power limit
  power.limit.part.data <- power.part.data
  if (!is.null(power.max.limit)) {
    power.limit.part.data <- sapply(power.limit.part.data, eval.max = power.max.limit, Limit.Max.inside)
  }
  if (!is.null(power.min.limit)) {
    power.limit.part.data <- sapply(power.limit.part.data, eval.min = power.min.limit, Limit.Min.inside)    
  }
  # plot data preparation
  pairs.plot.db <- data.frame(
    pair.name = names.part.data, 
    x = col.x.data, 
    y = col.y.data,
    cnt.orig = as.integer(cnt.part.data),
    cnt.limit = as.integer(cnt.limit.part.data),
    power.orig = power.part.data,
    power.limit = power.limit.part.data,
    stringsAsFactors = FALSE
  )
  # check if it needs x-axis text rotation
  if.need.x.axis.rotate <- FALSE
  x.axis.text.len <- sapply(pairs.plot.db$x, function(x) {
    nchar(as.character(x))
  })
  if (max(x.axis.text.len) > 4) {
    if.need.x.axis.rotate <- TRUE
  }
  # plot process
  plot.res <- ggplot(pairs.plot.db, aes(x, y))
  plot.res <- plot.res + labs(x = plot.axis.x.name, y = plot.axis.y.name)
  plot.res <- plot.res + geom_point(aes(size = cnt.limit, colour = power.limit)) + 
                         scale_x_discrete(limits = fac.x.clusters, breaks = fac.x.clusters) + 
                         scale_y_discrete(limits = fac.y.clusters, breaks = fac.y.clusters) + 
                         scale_size(name = "Count") + 
                         scale_colour_gradientn(name = "Power", colours = nodes.colour.seq, values = nodes.colour.value.seq)
  # add labels for those out of range ( > cnt.max.limit or < cnt.min.limit)
  part.on.cnt.db <- data.frame(x = pairs.plot.db$x, y = pairs.plot.db$y, cnt.orig = pairs.plot.db$cnt.orig, stringsAsFactors = FALSE)
  part.on.cnt.ex.max.db <- NULL
  if (!is.null(cnt.max.limit)) {
    part.on.cnt.ex.max.db <- part.on.cnt.db[which(part.on.cnt.db$cnt.orig > cnt.max.limit), ]
  }
  part.on.cnt.ex.min.db <- NULL
  if (!is.null(cnt.min.limit)) {
    part.on.cnt.ex.min.db <- part.on.cnt.db[which(part.on.cnt.db$cnt.orig < cnt.min.limit), ]
  }
  part.on.cnt.limit.db <- rbind(part.on.cnt.ex.max.db, part.on.cnt.ex.min.db)
  if (!is.null(part.on.cnt.limit.db) && !hide.cnt.label) {
    plot.res <- plot.res + 
                geom_text(aes(label = cnt.orig), data = part.on.cnt.limit.db,
                          hjust = user.label.cnt.options$hjust, vjust = user.label.cnt.options$vjust,
                          nudge_x = user.label.cnt.options$nudge.x, nudge_y = user.label.cnt.options$nudge.y,
                          size = user.label.cnt.options$size)
  }
  # add labels for those out of range (> power.max.limit or < power.min.limit)
  part.on.power.db <- data.frame(x = pairs.plot.db$x, y = pairs.plot.db$y, power.orig = pairs.plot.db$power.orig, stringsAsFactors = FALSE)
  part.on.power.ex.max.db <- NULL
  if (!is.null(power.max.limit)) {
    part.on.power.ex.max.db <- part.on.power.db[which(part.on.power.db$power.orig > power.max.limit), ]
  }
  part.on.power.ex.min.db <- NULL
  if (!is.null(power.min.limit)) {
    part.on.power.ex.min.db <- part.on.power.db[which(part.on.power.db$power.orig < power.min.limit), ]
  }
  part.on.power.limit.db <- rbind(part.on.power.ex.max.db, part.on.power.ex.min.db)
  if (!is.null(part.on.power.limit.db) && !hide.power.label) {
    plot.res <- plot.res + 
                geom_text(aes(label = round(power.orig, 2)), 
                          data = part.on.power.limit.db,
                          hjust = user.label.power.options$hjust, vjust = user.label.power.options$vjust,
                          nudge_x = user.label.power.options$nudge.x, nudge_y = user.label.power.options$nudge.y,
                          size = user.label.power.options$size)
  }
  plot.res <- plot.res + theme_classic()  # change theme
  # rotate labels in x-axis when needed
  if (if.need.x.axis.rotate) {
    plot.res <- plot.res + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  ## construct result tables
  # table raw.res
  res.raw.table <- data.frame(
    pair.name = pairs.plot.db$pair.name,
    x = pairs.plot.db$x,
    y = pairs.plot.db$y,
    cnt = pairs.plot.db$cnt.orig,
    power = pairs.plot.db$power.orig,
    stringsAsFactors = FALSE
  )
  # table cnt & power
  len.row <- length(fac.y.clusters)
  len.col <- length(fac.x.clusters)
  res.cnt.table <- matrix(data = c(0), nrow = len.row, ncol = len.col)
  res.power.table <- matrix(data = c(0.0), nrow = len.row, ncol = len.col)
  for (i in 1:len.col) {
    for (j in 1:len.row) {
      this.index <- len.row * (i - 1) + j
      res.cnt.table[j, i] <- pairs.plot.db$cnt.orig[this.index]
      res.power.table[j, i] <- round(pairs.plot.db$power.orig[this.index], 2)
    }
  }
  # set rownames & colnames
  rownames(res.cnt.table) <- as.character(fac.y.clusters)
  colnames(res.cnt.table) <- as.character(fac.x.clusters)
  rownames(res.power.table) <- as.character(fac.y.clusters)
  colnames(res.power.table) <- as.character(fac.x.clusters)
  #end# return 
  list(plot = plot.res, table = list(raw.res.table = res.raw.table, cnt.table = res.cnt.table, power.table = res.power.table))
}
