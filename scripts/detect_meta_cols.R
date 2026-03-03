#' @title Detect Common Metadata Columns for scRNA-seq Workflows
#'
#' @description
#' Convenience wrapper to auto-detect common scRNA-seq metadata columns in one call:
#' group/condition, sample/donor, and cluster/annotation. Uses
#' \code{detect_group_col()}, \code{detect_sample_col()}, and \code{detect_cluster_col()}.
#'
#' @param meta A data.frame (typically \code{obj@meta.data}).
#'
#' @return A named list:
#' \describe{
#'   \item{group_col}{Detected group/condition column name (or \code{NA_character_}).}
#'   \item{sample_col}{Detected sample/donor column name (or \code{NA_character_}).}
#'   \item{cluster_col}{Detected cluster/annotation column name (or \code{NA_character_}).}
#' }
#'
#' @examples
#' \dontrun{
#' meta <- combined_no_IA@meta.data
#' detect_meta_cols(meta)
#' }
#'
#' @export
detect_meta_cols <- function(meta) {
  stopifnot(is.data.frame(meta))

  list(
    group_col = detect_group_col(meta),
    sample_col = detect_sample_col(meta),
    cluster_col = detect_cluster_col(meta)
  )
}
