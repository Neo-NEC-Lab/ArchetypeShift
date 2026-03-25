#' @title Compute archetype occupancy (replicate-aware)
#'
#' @description
#' Computes archetype program occupancy as mean per-cell archetype weight aggregated
#' within biological samples (and optionally within sample × cluster groups).
#' This avoids pooled-cell artifacts by summarizing within samples first.
#'
#' @param obj Seurat object.
#' @param weight_cols Optional character vector of metadata column names containing per-cell
#'   archetype weights (e.g., c("A5_A1","A5_A2",...)). If NULL, attempts auto-detection.
#' @param weights_mat Optional numeric matrix/data.frame of per-cell archetype weights
#'   (cells × k). Rownames must be cell IDs in obj. If provided, this is used instead of
#'   `weight_cols`.
#' @param sample_col Metadata column defining biological samples (e.g., "orig.ident").
#'   If NULL, attempts auto-detection.
#' @param group_col Metadata column defining experimental groups (e.g., "condition2").
#'   If NULL, attempts auto-detection.
#' @param cluster_col Optional metadata column defining clusters/cell types
#'   (e.g., "merged_cluster_annotations"). If NULL, no cluster-stratified table is produced.
#' @param threshold Optional numeric in [0,1]. If provided, also computes fraction of cells
#'   per group with weight >= threshold for each archetype (suffix "_frac").
#' @param min_cells Integer >= 1. Drop sample (or sample×cluster) groups with fewer cells.
#'
#' @return A list with:
#'   - sample_occ: data.frame with one row per sample × group and mean weights (+ n_cells)
#'   - sample_cluster_occ: data.frame with one row per sample × group × cluster (if cluster_col not NULL)
#'   - spec: list describing detected columns and weight names
#' 
#' #' @examples
#' 
#' obj <- combined_no_IA
#' weight_cols_k5 <- paste0("A5_A", 1:5)
#'
#' # Replicate-aware occupancy summaries
#' occ_k5 <- compute_occupancy(
#'   obj = obj,
#'   weight_cols = weight_cols_k5,
#'   sample_col = "orig.ident",
#'   group_col = "condition2",
#'   cluster_col = "merged_cluster_annotations",
#'   threshold = NULL,
#'   min_cells = 10L
#' )
#'
#' @export
compute_occupancy <- function(
  obj,
  weight_cols = NULL,
  weights_mat = NULL,
  sample_col = NULL,
  group_col = NULL,
  cluster_col = NULL,
  threshold = NULL,
  min_cells = 1L
) {
  # -----------------------------
  # [CHECKS] object
  # -----------------------------
  if (!methods::is(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }
  md <- obj@meta.data

  # -----------------------------
  # [CHECKS] min_cells / threshold
  # -----------------------------
  if (!is.numeric(min_cells) || length(min_cells) != 1L || is.na(min_cells) || min_cells < 1) {
    stop("`min_cells` must be a single integer >= 1.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
        threshold < 0 || threshold > 1) {
      stop("`threshold` must be a single numeric value in [0, 1] or NULL.", call. = FALSE)
    }
  }

  # -----------------------------
  # [HELPER] column auto-detect
  # -----------------------------
  detect_first_present <- function(candidates, cols) {
    hit <- intersect(candidates, cols)
    if (length(hit) == 0L) return(NA_character_)
    hit[[1]]
  }

  # -----------------------------
  # [AUTO] detect sample/group cols (if NULL)
  # -----------------------------
  if (is.null(sample_col)) {
    sample_col <- detect_first_present(c("orig.ident", "sample", "Sample", "patient", "Patient"), colnames(md))
  }
  if (is.null(group_col)) {
    group_col <- detect_first_present(c("condition2", "condition", "group", "Group", "treatment", "Treatment"), colnames(md))
  }

  if (!is.character(sample_col) || length(sample_col) != 1L || is.na(sample_col) || !nzchar(sample_col)) {
    stop("Could not auto-detect `sample_col`. Please pass `sample_col` explicitly.", call. = FALSE)
  }
  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col) || !nzchar(group_col)) {
    stop("Could not auto-detect `group_col`. Please pass `group_col` explicitly.", call. = FALSE)
  }
  if (!sample_col %in% colnames(md)) {
    stop("`sample_col` not found in obj@meta.data: ", sample_col, call. = FALSE)
  }
  if (!group_col %in% colnames(md)) {
    stop("`group_col` not found in obj@meta.data: ", group_col, call. = FALSE)
  }

  if (!is.null(cluster_col)) {
    if (!is.character(cluster_col) || length(cluster_col) != 1L || is.na(cluster_col) || !nzchar(cluster_col)) {
      stop("`cluster_col` must be a single non-empty character string or NULL.", call. = FALSE)
    }
    if (!cluster_col %in% colnames(md)) {
      stop("`cluster_col` not found in obj@meta.data: ", cluster_col, call. = FALSE)
    }
  }

  # -----------------------------
  # [WEIGHTS] from weights_mat or weight_cols
  # -----------------------------
  if (!is.null(weights_mat)) {
    if (is.data.frame(weights_mat)) weights_mat <- as.matrix(weights_mat)
    if (!is.matrix(weights_mat) || !is.numeric(weights_mat)) {
      stop("`weights_mat` must be a numeric matrix/data.frame (cells × k).", call. = FALSE)
    }
    if (is.null(rownames(weights_mat)) || any(!nzchar(rownames(weights_mat)))) {
      stop("`weights_mat` must have non-empty rownames matching Seurat cell IDs.", call. = FALSE)
    }
    keep_cells <- intersect(rownames(md), rownames(weights_mat))
    if (length(keep_cells) == 0L) {
      stop("No overlapping cell IDs between obj@meta.data rownames and `weights_mat` rownames.", call. = FALSE)
    }

    w <- weights_mat[keep_cells, , drop = FALSE]
    if (is.null(colnames(w)) || any(!nzchar(colnames(w)))) {
      colnames(w) <- paste0("A", seq_len(ncol(w)))
    }
    w <- as.data.frame(w, stringsAsFactors = FALSE)

    df <- md[keep_cells, c(sample_col, group_col, if (!is.null(cluster_col)) cluster_col), drop = FALSE]
    df <- cbind(df, w, stringsAsFactors = FALSE)
    weight_names <- colnames(w)
  } else {
    if (is.null(weight_cols)) {
      # auto-detect common archetype weight column patterns
      pat1 <- "^A\\d+_A\\d+$"  # e.g., A5_A1
      pat2 <- "^A\\d+$"        # e.g., A1
      weight_cols <- colnames(md)[grepl(pat1, colnames(md))]
      if (length(weight_cols) == 0L) {
        weight_cols <- colnames(md)[grepl(pat2, colnames(md))]
      }
      if (length(weight_cols) == 0L) {
        stop("Could not auto-detect `weight_cols` in obj@meta.data. Please pass `weight_cols` explicitly.", call. = FALSE)
      }
    }

    if (!is.character(weight_cols) || length(weight_cols) < 1L) {
      stop("`weight_cols` must be a non-empty character vector or NULL (for auto-detect).", call. = FALSE)
    }
    weight_cols <- unique(weight_cols)
    missing_w <- setdiff(weight_cols, colnames(md))
    if (length(missing_w) > 0L) {
      stop("Missing `weight_cols` in obj@meta.data: ", paste(missing_w, collapse = ", "), call. = FALSE)
    }

    # ensure numeric
    for (wcol in weight_cols) {
      if (!is.numeric(md[[wcol]])) {
        suppressWarnings(xn <- as.numeric(md[[wcol]]))
        if (all(is.na(xn))) {
          stop("Weight column '", wcol, "' is not numeric and cannot be coerced.", call. = FALSE)
        }
        md[[wcol]] <- xn
      }
    }

    df <- md[, c(sample_col, group_col, if (!is.null(cluster_col)) cluster_col, weight_cols), drop = FALSE]
    weight_names <- weight_cols
  }

  # -----------------------------
  # [CHECKS] weights range (soft)
  # -----------------------------
  w_all <- unlist(df[, weight_names, drop = FALSE], use.names = FALSE)
  w_all <- w_all[is.finite(w_all)]
  if (length(w_all) > 0L) {
    if (min(w_all) < -1e-6 || max(w_all) > 1 + 1e-6) {
      warning("Some weights fall outside [0,1]. Proceeding, but verify your weight columns/matrix.")
    }
  }

  # -----------------------------
  # [COMPUTE] sample-level occupancy (mean weight)
  # -----------------------------
  df$.n <- 1L

  sample_occ <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(sample_col, group_col)))) %>%
    dplyr::summarise(
      n_cells = sum(.data$.n),
      dplyr::across(dplyr::all_of(weight_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_cells >= min_cells)

  # optional threshold fractions
  if (!is.null(threshold)) {
    frac_tbl <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(sample_col, group_col)))) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(weight_names),
          ~ mean(.x >= threshold, na.rm = TRUE),
          .names = "{.col}_frac"
        ),
        .groups = "drop"
      )
    sample_occ <- dplyr::left_join(sample_occ, frac_tbl, by = c(sample_col, group_col))
  }

  # -----------------------------
  # [COMPUTE] sample × cluster occupancy (optional)
  # -----------------------------
  sample_cluster_occ <- NULL
  if (!is.null(cluster_col)) {
    sample_cluster_occ <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(sample_col, group_col, cluster_col)))) %>%
      dplyr::summarise(
        n_cells = sum(.data$.n),
        dplyr::across(dplyr::all_of(weight_names), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$n_cells >= min_cells)

    if (!is.null(threshold)) {
      frac_tbl_sc <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(sample_col, group_col, cluster_col)))) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(weight_names),
            ~ mean(.x >= threshold, na.rm = TRUE),
            .names = "{.col}_frac"
          ),
          .groups = "drop"
        )
      sample_cluster_occ <- dplyr::left_join(
        sample_cluster_occ,
        frac_tbl_sc,
        by = c(sample_col, group_col, cluster_col)
      )
    }
  }

  # -----------------------------
  # [RETURN] stable outputs + spec
  # -----------------------------
  out <- list(
    sample_occ = as.data.frame(sample_occ, stringsAsFactors = FALSE),
    sample_cluster_occ = if (is.null(sample_cluster_occ)) NULL else as.data.frame(sample_cluster_occ, stringsAsFactors = FALSE),
    spec = list(
      sample_col = sample_col,
      group_col = group_col,
      cluster_col = cluster_col,
      weight_names = weight_names,
      threshold = threshold,
      min_cells = min_cells,
      n_cells_used = nrow(df)
    )
  )
  out
}


# =============================
# [SETUP] packages + inputs
# =============================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

obj <- combined_no_IA

# =============================
# [OPTION 1] Use archetype weights already stored in Seurat metadata
# (this is the most common workflow for you)
# =============================
# If your k=5 weights are named like: A5_A1 ... A5_A5
weight_cols_k5 <- paste0("A5_A", 1:5)

occ_k5 <- compute_occupancy(
  obj = obj,
  weight_cols = weight_cols_k5,
  sample_col = "orig.ident",
  group_col = "condition2",
  cluster_col = "merged_cluster_annotations",
  threshold = NULL,
  min_cells = 10L
)

# View outputs
print(head(occ_k5$sample_occ, 10))
print(head(occ_k5$sample_cluster_occ, 10))
str(occ_k5$spec)