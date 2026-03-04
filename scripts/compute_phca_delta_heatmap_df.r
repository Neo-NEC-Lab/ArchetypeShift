#' @title Compute Replicate-Aware Delta Occupancy for Heatmaps
#'
#' @description
#' Computes a long-format table of replicate-aware occupancy differences
#' (\code{group2 - group1}) for each \emph{cluster × archetype} combination.
#' Input data should be sample-level occupancy (e.g., output of
#' \code{compute_pcha_occupancy_sample_cluster()}), so that condition means are
#' computed across samples rather than pooled across cells.
#'
#' @param occ_sc A data.frame of sample-level occupancy with columns for sample,
#'   group/condition, cluster, and archetype weight columns.
#' @param sample_col Character scalar. Column name in \code{occ_sc} indicating samples.
#'   (Used for validation; deltas are computed after averaging within condition.)
#' @param group_col Character scalar. Column name in \code{occ_sc} indicating group/condition.
#' @param cluster_col Character scalar. Column name in \code{occ_sc} indicating cluster/cell type.
#' @param weight_cols Character vector of column names in \code{occ_sc} containing occupancy values
#'   (mean archetype weights per sample).
#' @param group1 Character scalar. Reference group (subtracted), e.g., \code{"NEC"}.
#' @param group2 Character scalar. Comparison group, e.g., \code{"NEC_HA"}.
#' @param drop_weight_prefix_regex Optional character scalar regex used to clean archetype
#'   names for display (default drops prefixes like \code{"A5_"} or \code{"A6_"}).
#'
#' @return A data.frame with columns \code{cluster}, \code{archetype}, and \code{delta},
#'   where \code{delta = mean(group2) - mean(group1)} computed across samples within each
#'   cluster.
#'
#' @details
#' This function first filters \code{occ_sc} to \code{group1} and \code{group2}. It then
#' averages occupancy values across samples within each \code{group_col × cluster_col}
#' stratum (replicate-aware). Finally, it computes \code{group2 - group1} per archetype
#' and returns a long-format table suitable for heatmap plotting.
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
#'
#' df_delta <- compute_pcha_delta_heatmap_df(
#'   occ_sc = occ_sc,
#'   sample_col = "orig.ident",
#'   group_col = "condition2",
#'   cluster_col = "merged_cluster_annotations",
#'   weight_cols = paste0("A5_A", 1:5),
#'   group1 = "NEC",
#'   group2 = "NEC_HA",
#'   drop_weight_prefix_regex = "^A5_"
#' )
#' head(df_delta)
#' }
#'
#' @export
compute_pcha_delta_heatmap_df <- function(
  occ_sc,
  sample_col,
  group_col,
  cluster_col,
  weight_cols,
  group1 = "NEC",
  group2 = "NEC_HA",
  drop_weight_prefix_regex = "^A\\d+_"
) {
  # -----------------------------
  # [CHECKS] inputs
  # -----------------------------
  if (!is.data.frame(occ_sc)) {
    stop("`occ_sc` must be a data.frame.", call. = FALSE)
  }

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

  if (!is.character(group1) || length(group1) != 1L || is.na(group1) || !nzchar(group1)) {
    stop("`group1` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(group2) || length(group2) != 1L || is.na(group2) || !nzchar(group2)) {
    stop("`group2` must be a single non-empty character string.", call. = FALSE)
  }
  if (identical(group1, group2)) {
    stop("`group1` and `group2` must be different.", call. = FALSE)
  }

  if (!is.null(drop_weight_prefix_regex)) {
    if (!is.character(drop_weight_prefix_regex) || length(drop_weight_prefix_regex) != 1L) {
      stop("`drop_weight_prefix_regex` must be a single regex string or NULL.", call. = FALSE)
    }
  }

  needed <- c(sample_col, group_col, cluster_col, weight_cols)
  missing_cols <- setdiff(needed, colnames(occ_sc))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in `occ_sc`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Ensure numeric occupancy columns
  for (w in weight_cols) {
    if (!is.numeric(occ_sc[[w]])) {
      suppressWarnings(xn <- as.numeric(occ_sc[[w]]))
      if (all(is.na(xn))) {
        stop("Occupancy column '", w, "' is not numeric and cannot be coerced to numeric.", call. = FALSE)
      }
      occ_sc[[w]] <- xn
    }
  }

  # -----------------------------
  # [FILTER] keep only the two groups
  # -----------------------------
  occ_sc2 <- occ_sc %>%
    dplyr::filter(.data[[group_col]] %in% c(group1, group2))

  if (nrow(occ_sc2) == 0L) {
    stop("No rows remain after filtering to group1/group2.", call. = FALSE)
  }

  # -----------------------------
  # [AGGREGATE] mean across samples within each group × cluster
  # -----------------------------
  mean_by_group <- occ_sc2 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_col, cluster_col)))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(weight_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # -----------------------------
  # [WIDE] separate columns by group
  # -----------------------------
  wide <- mean_by_group %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(group_col),
      values_from = dplyr::all_of(weight_cols),
      names_sep = "___"
    )

  # -----------------------------
  # [DELTA] build long table cluster × archetype
  # -----------------------------
  deltas <- lapply(weight_cols, function(w) {
    col_g1 <- paste0(w, "___", group1)
    col_g2 <- paste0(w, "___", group2)

    if (!(col_g1 %in% colnames(wide)) || !(col_g2 %in% colnames(wide))) {
      stop("Missing group columns for '", w, "'. Check that group1/group2 match values in ", group_col, ".", call. = FALSE)
    }

    archetype_label <- if (is.null(drop_weight_prefix_regex)) w else gsub(drop_weight_prefix_regex, "", w)

    data.frame(
      cluster = wide[[cluster_col]],
      archetype = archetype_label,
      delta = as.numeric(wide[[col_g2]] - wide[[col_g1]]),
      stringsAsFactors = FALSE
    )
  })

  df_delta <- do.call(rbind, deltas)

  # Basic NA guard
  if (any(is.na(df_delta$cluster))) {
    stop("Delta table contains NA clusters; check cluster_col values in occ_sc.", call. = FALSE)
  }

  df_delta
}