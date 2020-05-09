

#' Add expression values to feature genes
#'
#' @description
#' It uses the calculated gene-counts data from Seurat to calculate the expression levels of genes.
#'
#' @param markers.all.from.Seurat Data.frame. Feature genes that are from return value of \code{Seurat::FindAllMarkers()}.
#' @param object.from.Seurat Seurat object. Please follow the vignette in \url{https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html},
#' and must run the process at least after getting result of clustering algorithm (like tSNE/UMAP).
#' Then, in \code{SeuratObject@meta.data}, you can see one column named 'seurat_clusters', 
#' or you change the name of that column, and just give it in function parameter
#' @param column.clusters Character. The default value is "seurat_clustesrs", as defined by Seurat package in
#' its default settings. If you are not sure, please use \code{colnames(object.from.Seurat@meta.data)} to figure it out.
#' @param method.calc.gene.exprs Character. Either "mean" or "median".
#' It chooses which mathematical method to be used to calculate gene expression. 
#' Its default value is "mean", which means using average expression level.
#'
#'
#'
#' @importClassesFrom Seurat Seurat
#' @importFrom stats median
#'
#' @export
#'
DataPrep.AddExprsVals <- function(
  markers.all.from.Seurat,
  object.from.Seurat,
  column.clusters = "seurat_clusters", 
  method.calc.gene.exprs = "mean"
) {
  # check clusters names are the same
  clusters.in.markers <- levels(factor(markers.all.from.Seurat$cluster))
  clusters.in.metadata <- levels(factor(object.from.Seurat@meta.data[, column.clusters]))
  if (sum(clusters.in.markers %in% clusters.in.metadata) != length(clusters.in.markers)) {
    stop("Please recheck your given data, definition of clusters is not same between SeuratObject@meta.data and markers")
  }
  # get gene exprs counts
  clusters.factors <- clusters.in.markers
  counts.all <- object.from.Seurat[[object.from.Seurat@active.assay]]@counts  # count matrix, by default, colnames are cell barcodes, rownames are genes 
  cells.all.ident <- object.from.Seurat@active.ident  # gives Cell.Barcode --- clusters pairs 
  # calculate exprs from counts
  markers.all.result <- NULL
  for (onecluster in clusters.factors) {
    cells.sub.ident <- cells.all.ident[which(cells.all.ident == onecluster)]
    cells.sub.names <- names(cells.sub.ident)  # get all cell barcodes which are of `onecluster`
    inds.counts.sub <- match(cells.sub.names, colnames(counts.all))
    counts.subdata <- counts.all[, inds.counts.sub]  # get all gene exprs of these cells
    markers.onecluster <- markers.all.from.Seurat[which(markers.all.from.Seurat$cluster == onecluster), ]  # find marker genes of `onecluster`
    genes.subdata <- markers.onecluster$gene
    counts.subdata <- counts.subdata[match(genes.subdata, rownames(counts.subdata)), ]  # get marker genes and cells of `onecluster`
    # caculate gene expression levels
    if (method.calc.gene.exprs == "mean") {
      calc.exprs.subdata <- apply(counts.subdata, MARGIN = 1, mean)
    } else {
      if (method.calc.gene.exprs == "median") {
        calc.exprs.subdata <- apply(counts.subdata, MARGIN = 1, median)
      } else {
        stop(paste0("Error: undefined method name: ", method, "!"))
      }
    } 
    inds.onecluster.matchback <- match(rownames(counts.subdata), markers.onecluster$gene)  # match back by 'gene'
    markers.onecluster[inds.onecluster.matchback, "exprs"] <- calc.exprs.subdata
    # merge exprs to form new data.frame
    markers.all.result <- rbind(markers.all.result, markers.onecluster)
  }
  # return 
  markers.all.result  # colnames() = p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene exprs
}





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
  if (length(inds.raw.match) > 0) {
    markers.raw.match   <- markers.all.from.Seurat[inds.raw.match, ]
    markers.raw.unmatch <- markers.all.from.Seurat[-inds.raw.match, ]
  } else {  # if all are not matched, which is of rare condition.
    markers.raw.match   <- NULL
    markers.raw.unmatch <- markers.all.from.Seurat
  }
  if (nrow(markers.raw.unmatch) > 0) {  # some unmatches exist
    inds.map.match <- match(markers.raw.unmatch$gene, map.synonyms.db$Synonym.each)
    print(paste0("In unmatched ", nrow(markers.raw.unmatch), " genes, ", length(which(!is.na(inds.map.match))), " are remapped from synonyms!"))
    if (length(which(is.na(inds.map.match))) > 0) {  # even after remap, still not matchable, so just leave out, and report in warning()
      warning("These ", warning.given, " are not matched, \n: ", 
        paste0(markers.raw.unmatch$gene[which(is.na(inds.map.match))], collapse = ",  "), "."
      )
    }
    tmp.gene.name.use.old <- markers.raw.unmatch$gene[which(!is.na(inds.map.match))]
    tmp.gene.name.use.new <- map.synonyms.db$Symbol_from_nomenclature_authority[inds.map.match][which(!is.na(inds.map.match))]
    markers.raw.unmatch$gene[which(!is.na(inds.map.match))] <- tmp.gene.name.use.new  # set those matched genes
    # check if these remapped gene name are in dup.synonyms.ref
    logic.ifinddup <- which(tmp.gene.name.use.old %in% dup.synonyms.ref)
    if (length(logic.ifinddup) > 0)
      warning("Synonyms of these ", warning.given, " are duplicate with others, and here recommend to do mannual check-up. ", 
        "Pairs with old~new gene names are given. \n: ",
        paste0(paste(tmp.gene.name.use.old[logic.ifinddup], tmp.gene.name.use.new[logic.ifinddup], sep = "~"), collapse = ",  "), "."
      )
  }
  rbind(markers.raw.match, markers.raw.unmatch)
}





#' Replace cluster names for feature genes
#'
#' @description
#' It replaces the old cluster names of feature genes that are in integer generated automatically by \pkg{Seurat} algorithm.
#'
#' @param markers.all.from.Seurat Data.frame. Feature genes that are from return value of \code{Seurat::FindAllMarkers()}.
#' @param cluster.names.current Vector. It specifies the names of clusters used currently.
#' @param cluster.names.replace Vector. It gives the replaced names of clusters.
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
DataPrep.AddClusterName <- function(
  markers.all.from.Seurat,
  cluster.names.current,
  cluster.names.replace
) {
  if (length(cluster.names.current) != length(cluster.names.replace)) {
    stop("The replaced names are of different length of the current used ones.")
  }
  markers.all.from.Seurat[, "cluster.old"] <- markers.all.from.Seurat[, "cluster"]
  tmp.fac <- factor(markers.all.from.Seurat[, "cluster"])
  lvl.tmp.fac <- levels(tmp.fac)
  inds.match <- match(cluster.names.current, lvl.tmp.fac)
  if (length(which(is.na(inds.match))) != 0) {
    stop(paste0("Please give right current used cluster names! ",
      "Wrong given ones are ", paste0(cluster.names.current[which(is.na(inds.match))], collapse = ", "),
      "."))    
  }
  lvl.tmp.fac[inds.match] <- cluster.names.replace
  levels(tmp.fac) <- lvl.tmp.fac
  markers.all.from.Seurat[, "cluster"] <- as.character(tmp.fac)
  # return
  markers.all.from.Seurat
}
