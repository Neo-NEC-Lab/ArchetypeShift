#' Compute replicate-aware delta occupancy for heatmaps
#'
#' Computes a long-format table of archetype occupancy differences
#' (`group2 - group1`) for each cluster-archetype combination using sample-level
#' occupancy summaries.
#'
#' This function is intended for occupancy tables derived from
#' [compute_pcha_occupancy()] where occupancy has already been summarized at the
#' sample level within clusters. It averages occupancy across samples within each
#' group-cluster stratum and then computes `group2 - group1` per archetype.
#'
#' @param occ_sc A data frame of sample-level occupancy with columns for sample,
#'   group, cluster, and archetype occupancy values.
#' @param sample_col Character scalar giving the sample identifier column in
#'   `occ_sc`. This column is validated but not explicitly used in the delta
#'   calculation after grouping.
#' @param group_col Character scalar giving the group/condition column in
#'   `occ_sc`.
#' @param cluster_col Character scalar giving the cluster/cell-state column in
#'   `occ_sc`.
#' @param weight_cols Character vector of occupancy columns in `occ_sc`,
#'   typically mean archetype weights per sample.
#' @param group1 Character scalar giving the reference group to subtract, such as
#'   `"NEC"`.
#' @param group2 Character scalar giving the comparison group, such as
#'   `"NEC_HA"`.
#' @param drop_weight_prefix_regex Optional character scalar regex used to remove
#'   shared prefixes from archetype names for display, for example `"^A5_"`. Use
#'   `NULL` to retain original column names.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{cluster}{Cluster label.}
#'   \item{archetype}{Archetype label.}
#'   \item{delta}{Replicate-aware occupancy difference computed as
#'   `mean(group2) - mean(group1)`.}
#' }
#'
#' @details
#' The function first filters the occupancy table to the two requested groups.
#' It then averages occupancy values across samples within each
#' `group_col x cluster_col` stratum and computes `group2 - group1` separately
#' for each archetype. The output is intended for plotting with
#' [plot_pcha_delta_heatmap()].
#'
#' @examples
#' \dontrun{
#' occ_sc <- compute_pcha_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   group_cols = c("orig.ident", "condition2", "merged_cluster_annotations")
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
#'
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
  # [CHECKS] data frame and columns
  # -----------------------------
  if (!is.data.frame(occ_sc)) {
    stop("`occ_sc` must be a data.frame.", call. = FALSE)
  }

  for (nm in c("sample_col", "group_col", "cluster_col")) {
    val <- get(nm, inherits = FALSE)
    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      stop("`", nm, "` must be a single non-empty character string.", call. = FALSE)
    }
  }

  if (!is.character(weight_cols) || length(weight_cols) < 1L || any(is.na(weight_cols)) || any(!nzchar(weight_cols))) {
    stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
  }
  weight_cols <- unique(weight_cols)

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
    if (!is.character(drop_weight_prefix_regex) || length(drop_weight_prefix_regex) != 1L ||
        is.na(drop_weight_prefix_regex) || !nzchar(drop_weight_prefix_regex)) {
      stop("`drop_weight_prefix_regex` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }

  needed <- c(sample_col, group_col, cluster_col, weight_cols)
  missing_cols <- setdiff(needed, colnames(occ_sc))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in `occ_sc`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # -----------------------------
  # [COERCE] occupancy columns to numeric
  # -----------------------------
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
  # [FILTER] retain only the requested groups
  # -----------------------------
  occ_sc2 <- occ_sc %>%
    dplyr::filter(rlang::.data[[group_col]] %in% c(group1, group2))

  if (nrow(occ_sc2) == 0L) {
    stop("No rows remain after filtering to `group1` and `group2`.", call. = FALSE)
  }

  if (!all(c(group1, group2) %in% unique(as.character(occ_sc2[[group_col]])))) {
    stop("Both `group1` and `group2` must be present in `", group_col, "` after filtering.", call. = FALSE)
  }

  # -----------------------------
  # [AGGREGATE] mean across samples within group x cluster
  # -----------------------------
  mean_by_group <- occ_sc2 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_col, cluster_col)))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(weight_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # -----------------------------
  # [WIDE] split values by group
  # -----------------------------
  wide <- mean_by_group %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(group_col),
      values_from = dplyr::all_of(weight_cols),
      names_sep = "___"
    )

  # -----------------------------
  # [DELTA] build long table
  # -----------------------------
  deltas <- lapply(weight_cols, function(w) {
    col_g1 <- paste0(w, "___", group1)
    col_g2 <- paste0(w, "___", group2)

    if (!(col_g1 %in% colnames(wide)) || !(col_g2 %in% colnames(wide))) {
      stop(
        "Missing group columns for '", w,
        "'. Check that `group1` and `group2` match values in `", group_col, "`.",
        call. = FALSE
      )
    }

    archetype_label <- if (is.null(drop_weight_prefix_regex)) {
      w
    } else {
      gsub(drop_weight_prefix_regex, "", w)
    }

    data.frame(
      cluster = wide[[cluster_col]],
      archetype = archetype_label,
      delta = as.numeric(wide[[col_g2]] - wide[[col_g1]]),
      stringsAsFactors = FALSE
    )
  })

  df_delta <- do.call(rbind, deltas)

  if (nrow(df_delta) == 0L) {
    stop("No rows were produced in the delta table.", call. = FALSE)
  }
  if (any(is.na(df_delta$cluster))) {
    stop("Delta table contains NA clusters; check `cluster_col` values in `occ_sc`.", call. = FALSE)
  }

  df_delta
}
