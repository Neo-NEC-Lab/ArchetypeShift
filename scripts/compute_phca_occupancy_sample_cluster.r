#' @title Compute Archetype Occupancy per Sample and Cluster
#'
#' @description
#' Computes archetype program occupancy as the mean archetype weight within each
#' \emph{sample × condition × cluster} group. This is a replicate-aware summary that
#' avoids pooled-cell artifacts by aggregating first within biological samples.
#'
#' @param obj A Seurat object.
#' @param weight_cols Character vector of metadata column names containing archetype
#'   weights (e.g., \code{c("A5_A1","A5_A2",...)}). These columns must exist in
#'   \code{obj@meta.data} and should represent per-cell weights in \eqn{[0,1]}.
#' @param sample_col Character scalar. Metadata column defining biological samples
#'   (e.g., \code{"orig.ident"}).
#' @param group_col Character scalar. Metadata column defining experimental groups
#'   (e.g., \code{"condition2"}).
#' @param cluster_col Character scalar. Metadata column defining cluster/cell-type
#'   labels (e.g., \code{"merged_cluster_annotations"} or \code{"seurat_clusters"}).
#'
#' @return A data.frame with one row per \code{sample_col × group_col × cluster_col}
#'   combination containing \code{n_cells} and mean values for each \code{weight_cols}
#'   column.
#'
#' @details
#' This function summarizes per-cell archetype weights into replicate-aware occupancy
#' estimates. For downstream condition contrasts, users can further aggregate or test
#' these sample-level estimates across groups.
#'
#' @examples
#' \dontrun{
#' occ_sc <- compute_pcha_occupancy_sample_cluster(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   sample_col = "orig.ident",
#'   group_col = "condition2",
#'   cluster_col = "merged_cluster_annotations"
#' )
#' head(occ_sc)
#' }
#'
#' @export
compute_pcha_occupancy_sample_cluster <- function(
  obj,
  weight_cols,
  sample_col,
  group_col,
  cluster_col
) {
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
  weight_cols <- unique(weight_cols)

  for (nm in c("sample_col", "group_col", "cluster_col")) {
    val <- get(nm, inherits = FALSE)
    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      stop("`", nm, "` must be a single non-empty character string.", call. = FALSE)
    }
  }

  missing_w <- setdiff(weight_cols, colnames(md))
  if (length(missing_w) > 0L) {
    stop("Missing `weight_cols` in obj@meta.data: ", paste(missing_w, collapse = ", "), call. = FALSE)
  }

  missing_g <- setdiff(c(sample_col, group_col, cluster_col), colnames(md))
  if (length(missing_g) > 0L) {
    stop("Missing grouping columns in obj@meta.data: ", paste(missing_g, collapse = ", "), call. = FALSE)
  }

  # -----------------------------
  # [CHECKS] weight columns numeric-ish
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
  # [COMPUTE] mean weights within sample × group × cluster
  # -----------------------------
  df <- md[, c(sample_col, group_col, cluster_col, weight_cols), drop = FALSE]
  df$.n <- 1L

  out <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(sample_col, group_col, cluster_col)))) %>%
    dplyr::summarise(
      n_cells = sum(.data$.n),
      dplyr::across(dplyr::all_of(weight_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  as.data.frame(out, stringsAsFactors = FALSE)
}