#' @title Detect a Cluster/Annotation Column in Seurat Metadata
#'
#' @description
#' Attempts to auto-detect a cluster or cell-type annotation column from common
#' scRNA-seq metadata conventions (e.g., \code{"seurat_clusters"}). Includes
#' \code{"merged_cluster_annotations"} as a commonly used user-defined annotation column.
#'
#' @param meta A data.frame (typically \code{obj@meta.data}).
#'
#' @details
#' The function searches for the first matching column name in the following order:
#' \itemize{
#'   \item \code{"merged_cluster_annotations"}
#'   \item \code{"seurat_clusters"}
#'   \item \code{"cluster"}, \code{"cluster_id"}
#'   \item \code{"ident"}
#' }
#' If none are found, returns \code{NA_character_}.
#'
#' @return A character scalar column name if detected; otherwise \code{NA_character_}.
#'
#' @examples
#' \dontrun{
#' meta <- combined_no_IA@meta.data
#' detect_cluster_col(meta)
#' }
#'
#' @export
detect_cluster_col <- function(meta) {
  stopifnot(is.data.frame(meta))

  candidates <- c(
    "seurat_clusters",
    "cluster",
    "cluster_id",
    "ident"
  )

  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}
