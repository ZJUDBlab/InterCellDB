

#' Remap gene symbols to their authorized genenames
#'
#' @description
#' Use species-specific reference database to remap genes to their authorized genenames
#'
#' @param markers.all Data.frame. Feature genes that are generated from \code{Seurat::FindAllMarkers()} or 
#' similar functions in other packages.
#' @inheritParams Inside.DummyGenesRefDB
#' @param warning.given Character. It substitutes a small part of the warning sentences. In 
#' most cases, you just leave it still.
#' @param if.used.inside Logic. If used inside, some process will not run.
#'
#'
#'
#' @export
#'
DataPrep.RemapClustersMarkers  <- function(
  markers.all,
  genes.ref.db,
  warning.given = "markers",
  if.used.inside = FALSE
) {
  if (if.used.inside == FALSE) {
    # force removing factors
    markers.all$gene <- as.character(markers.all$gene)
    markers.all$cluster <- as.character(markers.all$cluster)
  }

  # set each database
  entrez.db <- genes.ref.db$gene.ncbi.db
  map.synonyms.db <- genes.ref.db$gene.synonyms.db

  # split genes to already authorized ones and un-authorized ones
  inds.raw.match <- which(markers.all$gene %in% entrez.db$Symbol_from_nomenclature_authority)
  markers.raw.match   <- markers.all[inds.raw.match, ]
  markers.raw.unmatch <- markers.all[setdiff(seq_len(nrow(markers.all)), inds.raw.match), ]
  proc.unmatch.genes <- unique(markers.raw.unmatch$gene)
  
  ## un-authorized genes have 2 ways to go:
  # 1. cannot mapping from synonyms, keep still
  # 2. mapping from synonyms
  #
  ret.d0.cannot.match <- ret.d1.map.to.diff <- character()
  unmatch.map.results <- data.frame()
  if (nrow(markers.raw.unmatch) > 0) {
    # 1> extract unmatch-able genes
    ret.d0.cannot.match <- setdiff(proc.unmatch.genes, map.synonyms.db$Synonym.each)
    
    # short warning
    if (length(ret.d0.cannot.match) > 0) {
      print(paste0("Get ", length(ret.d0.cannot.match), 
        " genes cannot be mapped from synonyms.",
        " See return value `$unmapping.genes`."))
    }

    # 2> mapping the rest
    genes.poss.map <- setdiff(proc.unmatch.genes, ret.d0.cannot.match)
    inds.map.possible <- which(map.synonyms.db$Synonym.each %in% genes.poss.map)
    
    # check un-identical synonyms. One synonyms could map to different authorized genes
    check.non.identical <- tapply(map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.possible], map.synonyms.db$Synonym.each[inds.map.possible], function(x) {x}, simplify = FALSE)
    inds.dups <- sapply(check.non.identical, function(x) { ifelse(length(x) > 1, TRUE, FALSE) })
    # short warning
    ret.d1.warning <- names(check.non.identical)[inds.dups]
    if (length(ret.d1.warning) > 0) {
      print(paste0("Get ", length(ret.d1.warning), 
        " synonyms that cannot be mapped to unified authorized genes.",
        " See return value `$one.synonym.map.to.some.genes`."))  
    }
    # detailed hint
    ret.d1.map.to.diff <- vapply(which(inds.dups == TRUE), all.non.identical = check.non.identical, 
      function(x, all.non.identical) {
        this.gene <- names(all.non.identical)[x]
        against.authorized.genes <- paste0("(", paste0(all.non.identical[[x]], collapse = ", "), ")")
        paste(this.gene, "-->", against.authorized.genes)
      },
      FUN.VALUE = character(1), USE.NAMES = FALSE
    )

    # check if mapping result the same as other mapped genes
    check.map.dups <- tapply(map.synonyms.db$Synonym.each[inds.map.possible], 
      map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.possible], 
      function(x) {x}, simplify = FALSE)
    inds.map.dups <- sapply(check.map.dups, function(x) { ifelse(length(x) > 1, TRUE, FALSE) })
    # short warning
    ret.d2.warning <- unique(as.character(unlist(lapply(which(inds.map.dups == TRUE), all.map.dups = check.map.dups, 
      function(x, all.map.dups) {
        all.map.dups[[x]]
      }
    ))))
    if (length(ret.d2.warning) > 0) {
      print(paste0("Get ", length(ret.d2.warning), 
        " synonyms that get overlap mapping genes with at least one other.",
        " See return value `$some.synonyms.map.to.one.gene`."))  
    }
    # detailed hint
    ret.d2.diff.map.one <- vapply(which(inds.map.dups == TRUE), all.map.dups = check.map.dups,
      function(x, all.map.dups) {
        this.tg <- names(all.map.dups)[x]
        from.synonyms <- paste0("(", paste0(all.map.dups[[x]], collapse = ", "), ")")
        paste(from.synonyms, "-->", this.tg)
      },
      FUN.VALUE = character(1), USE.NAMES = FALSE
    )

    # check if mapping to exist genes
    check.map.to.exist <- intersect(markers.raw.match$gene, names(check.map.dups))
    # short warning
    ret.d3.warning <- as.character(unlist(lapply(check.map.to.exist, all.map.dups = check.map.dups, 
      function(x, all.map.dups) {
        all.map.dups[[which(names(all.map.dups) == x)]]
      }
    )))
    if (length(ret.d3.warning) > 0) {
      print(paste0("Get ", length(ret.d3.warning), 
        " synonyms that get mapping to existing authorized genes.",
        " See return value `$synonyms.map.to.exist.gene`."))
    }
    # detailed hint
    ret.d3.map.to.exist <- vapply(check.map.to.exist, all.map.dups = check.map.dups, 
      function(x, all.map.dups) {
        this.to.map.genes <- all.map.dups[[which(names(all.map.dups) == x)]]
        orig.to.map <- paste0("(", paste0(this.to.map.genes, collapse = ", "), ")")
        paste(orig.to.map, "-->", x)
      },
      FUN.VALUE = character(1), USE.NAMES = FALSE
    )
    
    ## get mapping result
    # unmatch-able ones
    unmatch.map.result.0 <- data.frame(unmatched = ret.d0.cannot.match, match.res = ret.d0.cannot.match, stringsAsFactors = FALSE)
    # match-able ones, in default: the first matched gene name will be used 
    inds.map.matches <- match(genes.poss.map, map.synonyms.db$Synonym.each)
    unmatch.map.result.1 <- data.frame(unmatched = map.synonyms.db$Synonym.each[inds.map.matches], match.res = map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.matches], stringsAsFactors = FALSE)
    # collect all
    unmatch.map.results <- rbind(unmatch.map.result.0, unmatch.map.result.1)
  }
  
  # collect markers result after mapping
  markers.raw.unmatch.dummy <- left_join(markers.raw.unmatch[, "gene", drop = FALSE],
    unmatch.map.results, by = c("gene" = "unmatched"))
  markers.raw.unmatch$gene <- markers.raw.unmatch.dummy$match.res
  markers.all <- rbind(markers.raw.unmatch, markers.raw.match)

  # after remapping, genes get to be duplicate with existing ones, re-check if mapping result has duplicate genes
  fcheck.result <- lapply(unique(markers.all$cluster), ref.markers = markers.all,
    function(x, ref.markers) {
      this.c.markers <- ref.markers[which(ref.markers$cluster == x), ]
      this.c.len <- tapply(seq_along(this.c.markers$gene), this.c.markers$gene, length)
      this.f.dup.genes <- names(this.c.len)[which(this.c.len > 1)]
      # in default: remove the latter one in given data
      this.c.markers <- DoPartUnique(this.c.markers, match(c("gene", "cluster"), colnames(this.c.markers)))
      list(dup.genes = this.f.dup.genes, markers = this.c.markers)
    })
  names(fcheck.result) <- unique(markers.all$cluster)

  # ret 
  ret.markers.all <- bind_rows(lapply(fcheck.result, function(x) { x$markers }))
  # detailed hint
  ret.dx.fcheck.dup <- lapply(seq_along(fcheck.result), all.fcheck.res = fcheck.result,
    function(x, all.fcheck.res) {
      this.cluster <- names(all.fcheck.res)[x]
      this.dup.genes <- all.fcheck.res[[x]][["dup.genes"]]
      ret.detailed <- NA
      if (length(this.dup.genes) > 0) {
        ret.detailed <- paste(this.cluster, "~", paste0(this.dup.genes, collapse = ", "))
      }
      list(detailed = ret.detailed, raw = unique(this.dup.genes))
    }
  )
  un.fin.check <- as.character(sapply(ret.dx.fcheck.dup, function(x) { x$detailed }))
  un.fin.check <- un.fin.check[which(!is.na(un.fin.check))]
  if (length(un.fin.check) > 0) {
    warning("There remain genes group by clusters to be checked manually. ",
      paste0(un.fin.check, collapse = "; "),
      ". The program automatically remove duplicate ones in cluster scale.")
  }
  ret.d4.final.dup.genes <- lapply(ret.dx.fcheck.dup, function(x) { x$raw })
  names(ret.d4.final.dup.genes) <- names(fcheck.result)

  return(list(result = ret.markers.all, 
    unmapping.genes = ret.d0.cannot.match, 
    one.synonym.map.to.some.genes = ret.d1.map.to.diff, 
    some.synonyms.map.to.one.gene = ret.d2.diff.map.one, 
    synonyms.map.to.exist.gene = ret.d3.map.to.exist,
    after.map.dup.genes = ret.d4.final.dup.genes))
}







#' Replace cluster names for feature genes
#'
#' @description
#' It replaces the old cluster names of feature genes that are in integer generated automatically by \pkg{Seurat} algorithm.
#'
#' @param markers.all Data.frame. Feature genes that are from return value of \code{Seurat::FindAllMarkers()}.
#' @param cluster.names.current Vector. It specifies the names of clusters used currently.
#' @param cluster.names.replace Vector. It gives the replaced names of clusters.
#' @param colnames.cluster Character. The name of the column that defining the cluster.
#'
#' @details
#' The replace process can be applied on part of clusters, which means the length of parameters 
#' \code{cluster.names.current} & \code{cluster.names.replace} can be arbitary one with smaller or equal to the count of all clusters.
#' The replace process is taken one-by-one way, as a result, the length of \code{cluster.names.current} and \code{cluster.names.replace} 
#' must be the same, or the replacing will not be adequately applied. 
#'
#' @examples
#' \dontrun{
#'  DataPrep.ReplaceClusterName(markers.all,
#'    cluster.names.current = c(0, 1, 2),
#'    cluster.names.replace = c("Astrocytes", "Microglia", "Endothelial cells"))
#' }
#'
#'
#'
#' @export
#'
DataPrep.ReplaceClusterName <- function(
  markers.all,
  cluster.names.current,
  cluster.names.replace,
  colnames.cluster = "cluster"
) {
  if ((colnames.cluster %in% colnames(markers.all)) == FALSE) {
    stop("Selected column name defining clusters is not in the given data!")
  }
  if (length(cluster.names.current) != length(cluster.names.replace)) {
    stop("The replaced names are of different length of the current used ones.")
  }
  reserve.oldnames.col <- reserve.oldnames.col.proto <- paste(colnames.cluster, "oldv", sep = ".")
  for (try.i in 1:100) {
    reserve.oldnames.col <- paste(reserve.oldnames.col.proto, as.character(try.i), sep = ".")
    if ((reserve.oldnames.col %in% colnames(markers.all)) == FALSE) {
      break
    }
    if (try.i == 100) {
      stop("Cannot allocate proper colnames for old cluster names, program failed! Please check given data!")
    }
  }
  markers.all[, reserve.oldnames.col] <- markers.all[, colnames.cluster]
  tmp.fac <- factor(markers.all[, colnames.cluster])
  lvl.tmp.fac <- levels(tmp.fac)
  inds.match <- match(cluster.names.current, lvl.tmp.fac)
  if (length(which(is.na(inds.match))) != 0) {
    stop(paste0("Please give right current used cluster names! ",
      "Wrong given ones are ", paste0(cluster.names.current[which(is.na(inds.match))], collapse = ", "),
      "."))    
  }
  lvl.tmp.fac[inds.match] <- cluster.names.replace
  levels(tmp.fac) <- lvl.tmp.fac
  markers.all[, colnames.cluster] <- as.character(tmp.fac)
  # return
  markers.all
}
