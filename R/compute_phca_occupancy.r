#' @title Compute Archetype Occupancy by Metadata Groups
#'
#' @description
#' Computes archetype program occupancy as the mean archetype weight within each
#' metadata-defined group. This is a general-purpose aggregation utility used for
#' occupancy summaries (e.g., by condition, sample, or cell-type annotations).
#'
#' @param obj A Seurat object.
#' @param weight_cols Character vector of metadata column names containing archetype
#'   weights (e.g., \code{c("A5_A1","A5_A2",...)}). These columns must exist in
#'   \code{obj@meta.data}.
#' @param group_cols Character vector of metadata column names to group by
#'   (e.g., \code{c("condition2")} or \code{c("orig.ident","condition2")}).
#'
#' @return A data.frame with one row per unique combination of \code{group_cols},
#'   containing \code{n_cells} and mean values for each \code{weight_cols} column.
#'
#' @details
#' This function performs a pooled-cell mean within each group. For replicate-aware
#' inference, prefer aggregating by sample first (e.g., include \code{orig.ident} in
#' \code{group_cols}) and then perform statistical testing at the sample level.
#'
#' @examples
#' \dontrun{
#' # Mean occupancy by condition
#' occ_cond <- compute_pcha_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   group_cols = c("condition2")
#' )
#'
#' # Replicate-aware table (sample × condition)
#' occ_sample <- compute_pcha_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   group_cols = c("orig.ident", "condition2")
#' )
#' }
#'
#' @export
compute_pcha_occupancy <- function(obj, weight_cols, group_cols) {
  # -----------------------------
  # [CHECKS] object + metadata
  # -----------------------------
  if (!methods::is(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }
  md <- obj@meta.data

  # -----------------------------
  # [CHECKS] columns
  # -----------------------------
  if (!is.character(weight_cols) || length(weight_cols) < 1L) {
    stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(group_cols) || length(group_cols) < 1L) {
    stop("`group_cols` must be a non-empty character vector.", call. = FALSE)
  }

  weight_cols <- unique(weight_cols)
  group_cols <- unique(group_cols)

  missing_w <- setdiff(weight_cols, colnames(md))
  if (length(missing_w) > 0L) {
    stop("Missing `weight_cols` in obj@meta.data: ", paste(missing_w, collapse = ", "), call. = FALSE)
  }

  missing_g <- setdiff(group_cols, colnames(md))
  if (length(missing_g) > 0L) {
    stop("Missing `group_cols` in obj@meta.data: ", paste(missing_g, collapse = ", "), call. = FALSE)
  }

  # -----------------------------
  # [CHECKS] weights numeric-ish
  # -----------------------------
  for (w in weight_cols) {
    if (!is.numeric(md[[w]])) {
      suppressWarnings(xn <- as.numeric(md[[w]]))
      if (all(is.na(xn))) {
        stop("Weight column '", w, "' is not numeric and cannot be coerced to numeric.", call. = FALSE)
      }
      md[[w]] <- xn
    }
  }

  # -----------------------------
  # [COMPUTE] pooled mean weights within groups
  # -----------------------------
  df <- md[, c(group_cols, weight_cols), drop = FALSE]
  df$.n <- 1L

  out <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_cells = sum(.data$.n),
      dplyr::across(dplyr::all_of(weight_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  as.data.frame(out, stringsAsFactors = FALSE)
}