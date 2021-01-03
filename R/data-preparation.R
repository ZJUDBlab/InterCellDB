

#' Remap gene symbols to their authorized genenames
#'
#' @description
#' Use species-specific reference database to remap genes to their authorized genenames
#'
#' @param markers.all.from.Seurat Data.frame. Feature genes that are from return value of \code{Seurat::FindAllMarkers()}.
#' @inheritParams Inside.DummyGenesRefDB
#' @param warning.given Character. It substitutes a small part of the warning sentences. In 
#' most cases, you just leave it out.
#'
#'
#'
#' @export
#'
DataPrep.RemapClustersMarkers <- function(
  markers.all.from.Seurat,
  genes.ref.db,
  warning.given = "markers"
) {
  # set each database
  entrez.db <- genes.ref.db$gene.ncbi.db
  map.synonyms.db <- genes.ref.db$gene.synonyms.db
  dup.synonyms.ref <- genes.ref.db$gene.dup.synonyms.db$Synonym.each  # character
  # check if some genes are already authorized symbols
  inds.raw.match <- which(markers.all.from.Seurat$gene %in% entrez.db$Symbol_from_nomenclature_authority)
  markers.raw.match   <- markers.all.from.Seurat[inds.raw.match, ]
  markers.raw.unmatch <- markers.all.from.Seurat[setdiff(seq_len(nrow(markers.all.from.Seurat)), inds.raw.match), ]
  ret.final.unmatch <- character()
  ret.final.synonyms <- character()
  if (nrow(markers.raw.unmatch) > 0) {  # some unmatches exist
    inds.map.match <- match(markers.raw.unmatch$gene, map.synonyms.db$Synonym.each)
    print(paste0("In unmatched ", nrow(markers.raw.unmatch), " genes, ", length(which(!is.na(inds.map.match))), " are remapped from synonyms!"))
    if (length(which(is.na(inds.map.match))) > 0) {  # even after remap, still not matchable, so just leave out, and report in warning()
      warning("These ", warning.given, " are not matched, and total number is (", 
        length(which(is.na(inds.map.match))), ") :\n ", 
        paste0(markers.raw.unmatch$gene[which(is.na(inds.map.match))], collapse = ",  "), "."
      )
    }
    tmp.gene.name.use.old <- markers.raw.unmatch$gene[which(!is.na(inds.map.match))]
    tmp.gene.name.use.new <- map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.match][which(!is.na(inds.map.match))]
    markers.raw.unmatch$gene[which(!is.na(inds.map.match))] <- tmp.gene.name.use.new  # set those matched genes
    # check if these remapped gene name are in dup.synonyms.ref
    logic.ifinddup <- which(tmp.gene.name.use.old %in% dup.synonyms.ref)
    if (length(logic.ifinddup) > 0)
      warning("Synonyms of these (", length(logic.ifinddup), ") ", warning.given, 
        " are duplicate with others, and here recommend to do mannual check-up. ", 
        "Pairs with old~new gene names are given. \n: ",
        paste0(paste(tmp.gene.name.use.old[logic.ifinddup], tmp.gene.name.use.new[logic.ifinddup], sep = "~"), collapse = ",  "), "."
      )
    # for return values
    ret.final.unmatch <- markers.raw.unmatch$gene[which(is.na(inds.map.match))]
    ret.final.synonyms <- paste(tmp.gene.name.use.old[logic.ifinddup], tmp.gene.name.use.new[logic.ifinddup], sep = "~")
  }
  return(list(result = rbind(markers.raw.match, markers.raw.unmatch), 
    unmatched = ret.final.unmatch, synonyms.dup = ret.final.synonyms))
}





#' Replace cluster names for feature genes
#'
#' @description
#' It replaces the old cluster names of feature genes that are in integer generated automatically by \pkg{Seurat} algorithm.
#'
#' @param markers.all.from.Seurat Data.frame. Feature genes that are from return value of \code{Seurat::FindAllMarkers()}.
#' @param cluster.names.current Vector. It specifies the names of clusters used currently.
#' @param cluster.names.replace Vector. It gives the replaced names of clusters.
#' @param colnames.cluster Character. [TODO]
#'
#' @details
#' The replace process can be applied on part of clusters, which means the length of parameters 
#' \code{cluster.names.current} & \code{cluster.names.replace} can be arbitary one with smaller or equal to the count of all clusters.
#' The replace process is taken one-by-one way, as a result, the length of \code{cluster.names.current} and \code{cluster.names.replace} 
#' must be the same, or the replacing will not be adequately applied. 
#'
#' @examples
#' \dontrun{
#'  DataPrep.AddClusterName(markers.all.from.Seurat,
#'    cluster.names.current = c(0, 1, 2),
#'    cluster.names.replace = c("Astrocytes", "Microglia", "Endothelial cells"))
#' }
#'
#'
#'
#' @export
#'
DataPrep.AddClusterName <- function(  # [TODO] change this function
  markers.all.from.Seurat,
  cluster.names.current,
  cluster.names.replace,
  colnames.cluster = "cluster"
) {
  if ((colnames.cluster %in% colnames(markers.all.from.Seurat)) == FALSE) {
    stop("Selected column name defining clusters is not in the given data!")
  }
  if (length(cluster.names.current) != length(cluster.names.replace)) {
    stop("The replaced names are of different length of the current used ones.")
  }
  reserve.oldnames.col <- reserve.oldnames.col.proto <- paste(colnames.cluster, "oldv", sep = ".")
  for (try.i in 1:100) {
    reserve.oldnames.col <- paste(reserve.oldnames.col.proto, as.character(try.i), sep = ".")
    if ((reserve.oldnames.col %in% colnames(markers.all.from.Seurat)) == FALSE) {
      break
    }
    if (try.i == 100) {
      stop("Cannot allocate proper colnames for old cluster names, program failed! Please check given data!")
    }
  }
  markers.all.from.Seurat[, reserve.oldnames.col] <- markers.all.from.Seurat[, colnames.cluster]
  tmp.fac <- factor(markers.all.from.Seurat[, colnames.cluster])
  lvl.tmp.fac <- levels(tmp.fac)
  inds.match <- match(cluster.names.current, lvl.tmp.fac)
  if (length(which(is.na(inds.match))) != 0) {
    stop(paste0("Please give right current used cluster names! ",
      "Wrong given ones are ", paste0(cluster.names.current[which(is.na(inds.match))], collapse = ", "),
      "."))    
  }
  lvl.tmp.fac[inds.match] <- cluster.names.replace
  levels(tmp.fac) <- lvl.tmp.fac
  markers.all.from.Seurat[, colnames.cluster] <- as.character(tmp.fac)
  # return
  markers.all.from.Seurat
}
