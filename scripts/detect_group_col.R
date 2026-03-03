#' @title Detect a Group/Condition Column in Seurat Metadata
#'
#' @description
#' Attempts to auto-detect a condition/group column from common scRNA-seq metadata
#' conventions (e.g., \code{"condition"}, \code{"condition2"}, \code{"group"}, \code{"treatment"}).
#'
#' @param meta A data.frame (typically \code{obj@meta.data}).
#'
#' @details
#' The function searches for the first matching column name in the following order:
#' \itemize{
#'   \item \code{"condition"}, \code{"condition2"}
#'   \item \code{"group"}, \code{"Group"}
#'   \item \code{"treatment"}, \code{"stim"}, \code{"status"}, \code{"timepoint"}
#' }
#' If none are found, returns \code{NA_character_}.
#'
#' @return A character scalar column name if detected; otherwise \code{NA_character_}.
#'
#' @examples
#' \dontrun{
#' meta <- combined_no_IA@meta.data
#' detect_group_col(meta)
#' }
#'
#' @export
detect_group_col <- function(meta) {
  stopifnot(is.data.frame(meta))

  candidates <- c(
    "condition",
    "condition2",
    "group",
    "Group",
    "treatment",
    "stim",
    "status",
    "timepoint"
  )

  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}
