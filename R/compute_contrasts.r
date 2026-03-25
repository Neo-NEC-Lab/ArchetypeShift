#' @title Compute condition contrasts for archetype occupancy (replicate-aware)
#'
#' @description
#' Computes NEC vs NEC_HA (or any 2-group) contrasts using per-sample occupancy values
#' produced by \code{compute_occupancy()}. Supports global (sample-level) or
#' cluster-stratified (sample×cluster) occupancy tables.
#'
#' @param occ_df A data.frame such as \code{occ$sample_occ} or \code{occ$sample_cluster_occ}.
#' @param group_col Character scalar. Column in \code{occ_df} defining groups/conditions.
#' @param sample_col Character scalar. Column in \code{occ_df} defining biological replicates.
#' @param cluster_col Optional character scalar. Column defining clusters/cell-types. If NULL,
#'   contrasts are computed without cluster stratification. If provided, contrasts are
#'   computed within each cluster.
#' @param weight_cols Optional character vector of occupancy columns (e.g., \code{c("A5_A1",...)}).
#'   If NULL, attempts to auto-detect numeric archetype columns (and excludes \code{n_cells}).
#' @param group_levels Optional length-2 character vector specifying group order:
#'   \code{c(baseline, treated)}. Effect is computed as treated - baseline.
#' @param test Statistical test: \code{"welch_t"} (default) or \code{"wilcox"}.
#' @param p_adjust Multiple testing correction method passed to \code{p.adjust()} (default \code{"BH"}).
#' @param min_n_per_group Minimum number of replicates per group required to test (default 2).
#'
#' @return A data.frame with one row per archetype (and per cluster if \code{cluster_col} is set):
#'   group means, effect (treated - baseline), p-value, adjusted p-value, and n per group.
#' 
#' #' @examples
#' \dontrun{
#' res_global <- compute_contrast(
#'   occ_df = occ$sample_occ,
#'   group_col = "condition",
#'   sample_col = "sample_id",
#'   cluster_col = NULL,
#'   weight_cols = paste0("A", 1:5),
#'   group_levels = c("control", "treated"),
#'   test = "welch_t",
#'   p_adjust = "BH",
#'   min_n_per_group = 2L
#' )
#'
#' res_by_cluster <- compute_contrast(
#'   occ_df = occ$sample_cluster_occ,
#'   group_col = "condition",
#'   sample_col = "sample_id",
#'   cluster_col = "cluster",
#'   weight_cols = paste0("A", 1:5),
#'   group_levels = c("control", "treated"),
#'   test = "welch_t",
#'   p_adjust = "BH",
#'   min_n_per_group = 2L
#' )
#' }
#'
#' @export
compute_contrast <- function(
  occ_df,
  group_col,
  sample_col,
  cluster_col = NULL,
  weight_cols = NULL,
  group_levels = NULL,
  test = c("welch_t", "wilcox")[1],
  p_adjust = "BH",
  min_n_per_group = 2L
) {
  # -----------------------------
  # [CHECKS] basic inputs
  # -----------------------------
  if (!is.data.frame(occ_df)) {
    stop("`occ_df` must be a data.frame.", call. = FALSE)
  }
  for (nm in c("group_col", "sample_col")) {
    val <- get(nm, inherits = FALSE)
    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      stop("`", nm, "` must be a single non-empty character string.", call. = FALSE)
    }
    if (!val %in% colnames(occ_df)) {
      stop("Column not found in `occ_df`: ", val, call. = FALSE)
    }
  }
  if (!is.null(cluster_col)) {
    if (!is.character(cluster_col) || length(cluster_col) != 1L || is.na(cluster_col) || !nzchar(cluster_col)) {
      stop("`cluster_col` must be a single non-empty character string or NULL.", call. = FALSE)
    }
    if (!cluster_col %in% colnames(occ_df)) {
      stop("`cluster_col` not found in `occ_df`: ", cluster_col, call. = FALSE)
    }
  }

  if (!is.numeric(min_n_per_group) || length(min_n_per_group) != 1L || is.na(min_n_per_group) || min_n_per_group < 1) {
    stop("`min_n_per_group` must be a single integer >= 1.", call. = FALSE)
  }
  min_n_per_group <- as.integer(min_n_per_group)

  test <- match.arg(test, choices = c("welch_t", "wilcox"))

  # -----------------------------
  # [CHECKS] determine groups (must be 2)
  # -----------------------------
  g_raw <- as.character(occ_df[[group_col]])
  g_raw <- trimws(g_raw)
  g_uniq <- unique(g_raw[!is.na(g_raw) & nzchar(g_raw)])

  if (is.null(group_levels)) {
    if (length(g_uniq) != 2L) {
      stop(
        "`occ_df` must contain exactly 2 non-empty groups in `", group_col, "` or you must pass `group_levels`.\n",
        "Observed groups: ", paste(g_uniq, collapse = ", "),
        call. = FALSE
      )
    }
    group_levels <- g_uniq
  } else {
    if (!is.character(group_levels) || length(group_levels) != 2L || any(is.na(group_levels)) || any(!nzchar(group_levels))) {
      stop("`group_levels` must be a length-2 character vector: c(baseline, treated).", call. = FALSE)
    }
    if (!all(group_levels %in% g_uniq)) {
      stop(
        "`group_levels` must be present in `occ_df[[group_col]]`.\n",
        "Provided: ", paste(group_levels, collapse = ", "), "\n",
        "Observed: ", paste(g_uniq, collapse = ", "),
        call. = FALSE
      )
    }
  }
  baseline <- group_levels[[1]]
  treated  <- group_levels[[2]]

  # -----------------------------
  # [CHECKS] detect weight columns
  # -----------------------------
  if (is.null(weight_cols)) {
    # choose numeric columns excluding n_cells and grouping cols
    exclude <- c(group_col, sample_col, cluster_col, "n_cells")
    candidates <- setdiff(colnames(occ_df), exclude)
    is_num <- vapply(candidates, function(x) is.numeric(occ_df[[x]]), logical(1))
    weight_cols <- candidates[is_num]
    if (length(weight_cols) == 0L) {
      stop(
        "Could not auto-detect `weight_cols`. Pass `weight_cols` explicitly.",
        call. = FALSE
      )
    }
  } else {
    if (!is.character(weight_cols) || length(weight_cols) < 1L || any(is.na(weight_cols)) || any(!nzchar(weight_cols))) {
      stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
    }
    missing_w <- setdiff(weight_cols, colnames(occ_df))
    if (length(missing_w) > 0L) {
      stop("Missing `weight_cols` in `occ_df`: ", paste(missing_w, collapse = ", "), call. = FALSE)
    }
  }

  # -----------------------------
  # [HELPER] run a 2-group test on vectors
  # -----------------------------
  run_test <- function(x0, x1, test) {
    x0 <- x0[is.finite(x0)]
    x1 <- x1[is.finite(x1)]
    if (length(x0) < min_n_per_group || length(x1) < min_n_per_group) {
      return(list(p = NA_real_, stat = NA_real_))
    }
    if (test == "welch_t") {
      tt <- stats::t.test(x1, x0, var.equal = FALSE)  # treated vs baseline
      return(list(p = as.numeric(tt$p.value), stat = as.numeric(tt$statistic)))
    }
    if (test == "wilcox") {
      wt <- stats::wilcox.test(x1, x0, exact = FALSE)
      return(list(p = as.numeric(wt$p.value), stat = as.numeric(wt$statistic)))
    }
    list(p = NA_real_, stat = NA_real_)
  }

  # -----------------------------
  # [COMPUTE] contrasts (global or within cluster)
  # -----------------------------
  split_keys <- if (is.null(cluster_col)) {
    rep(".__all__", nrow(occ_df))
  } else {
    as.character(occ_df[[cluster_col]])
  }

  keys <- unique(split_keys)
  out_list <- vector("list", length(keys))
  names(out_list) <- keys

  for (ki in seq_along(keys)) {
    key <- keys[[ki]]
    idx <- which(split_keys == key)
    d <- occ_df[idx, , drop = FALSE]

    # baseline/treated rows
    i0 <- which(as.character(d[[group_col]]) == baseline)
    i1 <- which(as.character(d[[group_col]]) == treated)

    # ensure one row per sample per group is already the case; we trust compute_occupancy output
    # (if not, user should aggregate first)

    res_k <- vector("list", length(weight_cols))

    for (j in seq_along(weight_cols)) {
      wnm <- weight_cols[[j]]
      x0 <- d[[wnm]][i0]
      x1 <- d[[wnm]][i1]

      n0 <- sum(is.finite(x0))
      n1 <- sum(is.finite(x1))

      m0 <- if (n0 > 0L) mean(x0, na.rm = TRUE) else NA_real_
      m1 <- if (n1 > 0L) mean(x1, na.rm = TRUE) else NA_real_

      eff <- m1 - m0

      tst <- run_test(x0, x1, test = test)

      res_k[[j]] <- data.frame(
        archetype = wnm,
        mean_baseline = m0,
        mean_treated = m1,
        effect_treated_minus_baseline = eff,
        n_baseline = n0,
        n_treated = n1,
        stat = tst$stat,
        p_value = tst$p,
        stringsAsFactors = FALSE
      )
    }

    res_k <- do.call(rbind, res_k)

    # adjust p-values within this stratum
    res_k$p_adj <- stats::p.adjust(res_k$p_value, method = p_adjust)

    if (!is.null(cluster_col)) {
      res_k[[cluster_col]] <- key
    }

    res_k$group_baseline <- baseline
    res_k$group_treated  <- treated
    res_k$test <- test
    res_k$p_adjust <- p_adjust

    out_list[[ki]] <- res_k
  }

  out <- do.call(rbind, out_list)

  # -----------------------------
  # [RETURN] stable column order
  # -----------------------------
  if (is.null(cluster_col)) {
    out <- out[, c(
      "archetype",
      "group_baseline", "group_treated",
      "mean_baseline", "mean_treated",
      "effect_treated_minus_baseline",
      "n_baseline", "n_treated",
      "stat", "p_value", "p_adj",
      "test", "p_adjust"
    ), drop = FALSE]
  } else {
    out <- out[, c(
      cluster_col,
      "archetype",
      "group_baseline", "group_treated",
      "mean_baseline", "mean_treated",
      "effect_treated_minus_baseline",
      "n_baseline", "n_treated",
      "stat", "p_value", "p_adj",
      "test", "p_adjust"
    ), drop = FALSE]
  }

  out
}

# =============================
# [EXAMPLE] Use with your compute_occupancy() output
# =============================
# Global (sample-level) archetype shifts
# res_global <- compute_contrast(
#   occ_df = occ_k5$sample_occ,
#   group_col = "condition2",
#   sample_col = "orig.ident",
#   cluster_col = NULL,
#   weight_cols = paste0("A5_A", 1:5),
#   group_levels = c("NEC", "NEC_HA"),
#   test = "welch_t",
#   p_adjust = "BH",
#   min_n_per_group = 2L
# )
#
# Cluster-stratified archetype shifts
# res_by_cluster <- compute_contrast(
#   occ_df = occ_k5$sample_cluster_occ,
#   group_col = "condition2",
#   sample_col = "orig.ident",
#   cluster_col = "merged_cluster_annotations",
#   weight_cols = paste0("A5_A", 1:5),
#   group_levels = c("NEC", "NEC_HA"),
#   test = "welch_t",
#   p_adjust = "BH",
#   min_n_per_group = 2L
# )