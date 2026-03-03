#' @title Detect a Sample Column in Seurat Metadata
#'
#' @description
#' Attempts to auto-detect a sample identifier column from common scRNA-seq
#' metadata conventions (e.g., \code{"orig.ident"}, \code{"sample"}, \code{"donor"}).
#' This is useful for pipelines that require per-sample grouping but want to remain
#' robust across different object conventions.
#'
#' @param meta A data.frame (typically \code{obj@meta.data}).
#'
#' @details
#' The function searches for the first matching column name in the following order:
#' \itemize{
#'   \item \code{"orig.ident"}
#'   \item \code{"sample"}, \code{"sample_id"}
#'   \item \code{"donor"}, \code{"donor_id"}
#'   \item \code{"batch"}, \code{"batch_id"}
#'   \item \code{"patient"}, \code{"patient_id"}
#' }
#' If none are found, returns \code{NA_character_}.
#'
#' @return A character scalar column name if detected; otherwise \code{NA_character_}.
#'
#' @examples
#' \dontrun{
#' meta <- combined_no_IA@meta.data
#' detect_sample_col(meta)
#' }
#'
#' @export
detect_sample_col <- function(meta) {
  stopifnot(is.data.frame(meta))

  candidates <- c(
    "orig.ident",
    "sample",
    "sample_id",
    "donor",
    "donor_id",
    "batch",
    "batch_id",
    "patient",
    "patient_id"
  )

  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}
