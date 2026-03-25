#' Plot an archetype occupancy heatmap
#'
#' Creates a heatmap showing archetype occupancy across a grouping variable such
#' as cluster, cell type, or condition. This function is intended for occupancy
#' summaries produced by [compute_pcha_occupancy()].
#'
#' @param occ_df A data frame containing the grouping column and archetype
#'   occupancy columns.
#' @param group_col Character scalar giving the column in `occ_df` to use on the
#'   y-axis, for example `"merged_cluster_annotations"`.
#' @param weight_cols Character vector of occupancy columns in `occ_df`, such as
#'   `c("A5_A1", "A5_A2", ...)`.
#' @param title Optional character scalar used as the plot title. If `NULL`, a
#'   default title is generated.
#' @param drop_weight_prefix_regex Optional character scalar regex used to remove
#'   a shared prefix from archetype labels for display, for example `"^A5_"`.
#'   Use `NULL` to disable relabeling.
#' @param group_order Optional character vector specifying the display order of
#'   `group_col` levels on the y-axis. Rows not present in `group_order` are
#'   dropped.
#' @param out_png Optional character path. If supplied, the plot is saved as a
#'   PNG file.
#' @param width Numeric scalar giving the saved PNG width in inches. Default is
#'   `7.5`.
#' @param height Numeric scalar giving the saved PNG height in inches. Default is
#'   `6.5`.
#' @param dpi Integer scalar giving the saved PNG resolution. Default is `300`.
#'
#' @return A `ggplot2` heatmap object.
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
#' p <- plot_pcha_occupancy_heatmap(
#'   occ_df = occ_res$sample_cluster_occ,
#'   group_col = "merged_cluster_annotations",
#'   weight_cols = paste0("A5_A", 1:5),
#'   drop_weight_prefix_regex = "^A5_",
#'   group_order = cluster_order,
#'   title = "A5 program occupancy by cell-type annotation"
#' )
#'
#' p
#' }
#'
#' @export
plot_pcha_occupancy_heatmap <- function(
  occ_df,
  group_col,
  weight_cols,
  title = NULL,
  drop_weight_prefix_regex = "^A\\d+_",
  group_order = NULL,
  out_png = NULL,
  width = 7.5,
  height = 6.5,
  dpi = 300
) {
  # -----------------------------
  # [CHECKS] inputs
  # -----------------------------
  if (!is.data.frame(occ_df)) {
    stop("`occ_df` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col) || !nzchar(group_col)) {
    stop("`group_col` must be a single non-empty character string.", call. = FALSE)
  }
  if (!(group_col %in% colnames(occ_df))) {
    stop("`group_col` not found in `occ_df`: ", group_col, call. = FALSE)
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

  if (!is.null(group_order)) {
    if (!is.character(group_order) || length(group_order) < 1L || any(is.na(group_order)) || any(!nzchar(group_order))) {
      stop("`group_order` must be NULL or a non-empty character vector.", call. = FALSE)
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
  df_long <- occ_df %>%
    dplyr::select(dplyr::all_of(c(group_col, weight_cols))) %>%
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
  # [ORDER] group axis
  # -----------------------------
  if (!is.null(group_order)) {
    df_long[[group_col]] <- factor(df_long[[group_col]], levels = group_order)
    df_long <- df_long[!is.na(df_long[[group_col]]), , drop = FALSE]
    if (nrow(df_long) == 0L) {
      stop("No rows remain after applying `group_order`.", call. = FALSE)
    }
  }

  if (is.null(title)) {
    title <- paste0("Archetype occupancy by ", group_col)
  }

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = rlang::.data$archetype,
      y = rlang::.data[[group_col]],
      fill = rlang::.data$occupancy
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title = title,
      x = "Archetype program",
      y = group_col,
      fill = "Mean weight"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # -----------------------------
  # [SAVE]
  # -----------------------------
  if (!is.null(out_png)) {
    ggplot2::ggsave(out_png, p, width = width, height = height, dpi = dpi)
  }

  p
}
