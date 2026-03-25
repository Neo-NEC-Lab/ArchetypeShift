#' Plot stacked archetype occupancy bars
#'
#' Creates stacked bar plots showing archetype occupancy across a grouping
#' variable such as sample, cluster, or condition. This function is intended for
#' occupancy summaries produced by [compute_pcha_occupancy()].
#'
#' @param occ_df A data frame containing grouping columns and archetype
#'   occupancy columns.
#' @param x_col Character scalar giving the column in `occ_df` to use on the
#'   x-axis, for example `"orig.ident"` or `"condition2"`.
#' @param weight_cols Character vector of occupancy columns in `occ_df`, such as
#'   `c("A5_A1", "A5_A2", ...)`.
#' @param title Optional character scalar used as the plot title. If `NULL`, a
#'   default title is generated.
#' @param drop_weight_prefix_regex Optional character scalar regex used to remove
#'   a shared prefix from archetype labels for display, for example `"^A5_"`.
#'   Use `NULL` to disable relabeling.
#' @param x_order Optional character vector specifying the display order of
#'   `x_col` levels. Levels not present in `x_order` are dropped.
#' @param facet_col Optional character scalar. If supplied, bars are faceted by
#'   this column, for example `"condition2"`.
#' @param out_png Optional character path. If supplied, the plot is saved as a
#'   PNG file.
#' @param width Numeric scalar giving the saved PNG width in inches. Default is
#'   `9.5`.
#' @param height Numeric scalar giving the saved PNG height in inches. Default is
#'   `4.5`.
#' @param dpi Integer scalar giving the saved PNG resolution. Default is `300`.
#'
#' @return A `ggplot2` stacked bar plot object.
#'
#' @details
#' This plot is typically used for descriptive summaries of archetype occupancy.
#' For replicate-aware inference, compute occupancy at the sample level and test
#' sample-level values separately.
#'
#' @examples
#' \dontrun{
#' occ_res <- compute_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   sample_col = "orig.ident",
#'   group_col = "condition2",
#'   cluster_col = "merged_cluster_annotations",
#'   min_cells = 10L
#' )
#'
#' p <- plot_pcha_occupancy_stacked(
#'   occ_df = occ_res$sample_occ,
#'   x_col = "orig.ident",
#'   weight_cols = paste0("A5_A", 1:5),
#'   facet_col = "condition2",
#'   drop_weight_prefix_regex = "^A5_",
#'   title = "A5 program occupancy by sample"
#' )
#'
#' p
#' }
#'
#' @export
plot_pcha_occupancy_stacked <- function(
  occ_df,
  x_col,
  weight_cols,
  title = NULL,
  drop_weight_prefix_regex = "^A\\d+_",
  x_order = NULL,
  facet_col = NULL,
  out_png = NULL,
  width = 9.5,
  height = 4.5,
  dpi = 300
) {
  # -----------------------------
  # [CHECKS] inputs
  # -----------------------------
  if (!is.data.frame(occ_df)) {
    stop("`occ_df` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(x_col) || length(x_col) != 1L || is.na(x_col) || !nzchar(x_col)) {
    stop("`x_col` must be a single non-empty character string.", call. = FALSE)
  }
  if (!(x_col %in% colnames(occ_df))) {
    stop("`x_col` not found in `occ_df`: ", x_col, call. = FALSE)
  }

  if (!is.character(weight_cols) || length(weight_cols) < 1L || any(is.na(weight_cols)) || any(!nzchar(weight_cols))) {
    stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
  }
  weight_cols <- unique(weight_cols)

  missing_w <- setdiff(weight_cols, colnames(occ_df))
  if (length(missing_w) > 0L) {
    stop("Missing `weight_cols` in `occ_df`: ", paste(missing_w, collapse = ", "), call. = FALSE)
  }

  if (!is.null(title)) {
    if (!is.character(title) || length(title) != 1L || is.na(title) || !nzchar(title)) {
      stop("`title` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }

  if (!is.null(drop_weight_prefix_regex)) {
    if (!is.character(drop_weight_prefix_regex) || length(drop_weight_prefix_regex) != 1L ||
        is.na(drop_weight_prefix_regex) || !nzchar(drop_weight_prefix_regex)) {
      stop("`drop_weight_prefix_regex` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }

  if (!is.null(x_order)) {
    if (!is.character(x_order) || length(x_order) < 1L || any(is.na(x_order)) || any(!nzchar(x_order))) {
      stop("`x_order` must be NULL or a non-empty character vector.", call. = FALSE)
    }
  }

  if (!is.null(facet_col)) {
    if (!is.character(facet_col) || length(facet_col) != 1L || is.na(facet_col) || !nzchar(facet_col)) {
      stop("`facet_col` must be NULL or a single non-empty character string.", call. = FALSE)
    }
    if (!(facet_col %in% colnames(occ_df))) {
      stop("`facet_col` not found in `occ_df`: ", facet_col, call. = FALSE)
    }
  }

  if (!is.null(out_png)) {
    if (!is.character(out_png) || length(out_png) != 1L || is.na(out_png) || !nzchar(out_png)) {
      stop("`out_png` must be NULL or a single non-empty character path.", call. = FALSE)
    }
  }

  if (!is.numeric(width) || length(width) != 1L || is.na(width) || width <= 0) {
    stop("`width` must be a single numeric value > 0.", call. = FALSE)
  }
  if (!is.numeric(height) || length(height) != 1L || is.na(height) || height <= 0) {
    stop("`height` must be a single numeric value > 0.", call. = FALSE)
  }
  if (!is.numeric(dpi) || length(dpi) != 1L || is.na(dpi) || dpi < 1) {
    stop("`dpi` must be a single numeric value >= 1.", call. = FALSE)
  }
  dpi <- as.integer(dpi)

  # -----------------------------
  # [COERCE] occupancy columns to numeric
  # -----------------------------
  for (w in weight_cols) {
    if (!is.numeric(occ_df[[w]])) {
      suppressWarnings(xn <- as.numeric(occ_df[[w]]))
      if (all(is.na(xn))) {
        stop("Occupancy column '", w, "' is not numeric and cannot be coerced to numeric.", call. = FALSE)
      }
      occ_df[[w]] <- xn
    }
  }

  # -----------------------------
  # [LONG] reshape to long format
  # -----------------------------
  keep_cols <- c(x_col, weight_cols)
  if (!is.null(facet_col)) {
    keep_cols <- c(x_col, facet_col, weight_cols)
  }

  df_long <- occ_df %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(weight_cols),
      names_to = "archetype",
      values_to = "occupancy"
    )

  # -----------------------------
  # [LABELS] clean archetype labels
  # -----------------------------
  if (is.null(drop_weight_prefix_regex)) {
    archetype_levels <- weight_cols
  } else {
    archetype_levels <- gsub(drop_weight_prefix_regex, "", weight_cols)
    df_long$archetype <- gsub(drop_weight_prefix_regex, "", df_long$archetype)
  }
  df_long$archetype <- factor(df_long$archetype, levels = archetype_levels)

  # -----------------------------
  # [ORDER] x-axis
  # -----------------------------
  if (!is.null(x_order)) {
    df_long[[x_col]] <- factor(df_long[[x_col]], levels = x_order)
    df_long <- df_long[!is.na(df_long[[x_col]]), , drop = FALSE]
    if (nrow(df_long) == 0L) {
      stop("No rows remain after applying `x_order`.", call. = FALSE)
    }
  }

  if (is.null(title)) {
    title <- paste0("Archetype occupancy by ", x_col)
  }

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = rlang::.data[[x_col]],
      y = rlang::.data$occupancy,
      fill = rlang::.data$archetype
    )
  ) +
    ggplot2::geom_col(position = "stack", width = 0.85) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title = title,
      x = x_col,
      y = "Mean weight",
      fill = "Archetype"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_col)), scales = "free_x")
  }

  # -----------------------------
  # [SAVE]
  # -----------------------------
  if (!is.null(out_png)) {
    ggplot2::ggsave(out_png, p, width = width, height = height, dpi = dpi)
  }

  p
}
