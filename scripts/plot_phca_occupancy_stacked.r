#' @title Plot Stacked Archetype Occupancy Bars
#'
#' @description
#' Plots stacked bar charts of archetype program occupancy from an occupancy table
#' (e.g., output of \code{compute_pcha_occupancy()}). Each bar represents one level
#' of \code{x_col} and stacks mean archetype weights across \code{weight_cols}.
#'
#' @param occ_df A data.frame containing grouping columns and archetype occupancy columns.
#' @param x_col Character scalar. Column in \code{occ_df} to use on the x-axis (e.g., \code{"orig.ident"}
#'   for sample-level bars or \code{"condition2"} for merged condition bars).
#' @param weight_cols Character vector of occupancy columns in \code{occ_df} (e.g., \code{c("A5_A1","A5_A2",...)}).
#' @param title Optional character scalar plot title. If NULL, a default title is generated.
#' @param drop_weight_prefix_regex Optional character scalar regex used to clean archetype labels
#'   for display (e.g., \code{"^A5_"}). Use NULL to disable.
#' @param x_order Optional character vector specifying the order of x-axis levels. Levels not
#'   present in \code{x_order} are dropped.
#' @param facet_col Optional character scalar. If provided, facets bars by this column (e.g.,
#'   \code{"condition2"} to facet samples by condition).
#' @param out_png Optional character path. If provided, saves the plot as a PNG.
#' @param width Numeric scalar. Width in inches for saved PNG.
#' @param height Numeric scalar. Height in inches for saved PNG.
#' @param dpi Integer scalar. Resolution for saved PNG.
#'
#' @return A ggplot object.
#'
#' @details
#' This plot is often used for descriptive summaries. For replicate-aware inference, compute
#' sample-level occupancy (include sample in grouping) and perform statistical testing on
#' sample-level values.
#'
#' @examples
#' \dontrun{
#' occ_sample <- compute_pcha_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   group_cols = c("orig.ident", "condition2")
#' )
#'
#' p <- plot_pcha_occupancy_stacked(
#'   occ_df = occ_sample,
#'   x_col = "orig.ident",
#'   weight_cols = paste0("A5_A", 1:5),
#'   facet_col = "condition2",
#'   drop_weight_prefix_regex = "^A5_",
#'   title = "A5 program occupancy by sample"
#' )
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
    stop("`x_col` not found in occ_df: ", x_col, call. = FALSE)
  }

  if (!is.character(weight_cols) || length(weight_cols) < 1L) {
    stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
  }
  weight_cols <- unique(weight_cols)

  missing_w <- setdiff(weight_cols, colnames(occ_df))
  if (length(missing_w) > 0L) {
    stop("Missing `weight_cols` in occ_df: ", paste(missing_w, collapse = ", "), call. = FALSE)
  }

  if (!is.null(facet_col)) {
    if (!is.character(facet_col) || length(facet_col) != 1L || is.na(facet_col) || !nzchar(facet_col)) {
      stop("`facet_col` must be a single non-empty character string or NULL.", call. = FALSE)
    }
    if (!(facet_col %in% colnames(occ_df))) {
      stop("`facet_col` not found in occ_df: ", facet_col, call. = FALSE)
    }
  }

  if (!is.null(drop_weight_prefix_regex)) {
    if (!is.character(drop_weight_prefix_regex) || length(drop_weight_prefix_regex) != 1L) {
      stop("`drop_weight_prefix_regex` must be a single regex string or NULL.", call. = FALSE)
    }
  }

  # Coerce occupancy columns to numeric if needed
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
  # [LONG] reshape to archetype/value rows
  # -----------------------------
  keep_cols <- c(x_col, weight_cols)
  if (!is.null(facet_col)) keep_cols <- c(x_col, facet_col, weight_cols)

  df_long <- occ_df %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(weight_cols), names_to = "archetype", values_to = "occupancy")

  # Clean archetype labels
  if (is.null(drop_weight_prefix_regex)) {
    archetype_levels <- weight_cols
    df_long$archetype <- df_long$archetype
  } else {
    archetype_levels <- gsub(drop_weight_prefix_regex, "", weight_cols)
    df_long$archetype <- gsub(drop_weight_prefix_regex, "", df_long$archetype)
  }
  df_long$archetype <- factor(df_long$archetype, levels = archetype_levels)

  # -----------------------------
  # [ORDER] x-axis
  # -----------------------------
  if (!is.null(x_order)) {
    if (!is.character(x_order) || length(x_order) < 1L) {
      stop("`x_order` must be a non-empty character vector or NULL.", call. = FALSE)
    }
    df_long[[x_col]] <- factor(df_long[[x_col]], levels = x_order)
    df_long <- df_long[!is.na(df_long[[x_col]]), , drop = FALSE]
  }

  if (is.null(title)) title <- paste0("Archetype occupancy by ", x_col)

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[x_col]], y = .data$occupancy, fill = .data$archetype)) +
    ggplot2::geom_col(position = "stack", width = 0.85) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(title = title, x = x_col, y = "Mean weight", fill = "Archetype") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_col)), scales = "free_x")
  }

  # -----------------------------
  # [SAVE] PNG if requested
  # -----------------------------
  if (!is.null(out_png)) {
    if (!is.character(out_png) || length(out_png) != 1L || is.na(out_png) || !nzchar(out_png)) {
      stop("`out_png` must be a single non-empty character path or NULL.", call. = FALSE)
    }
    ggplot2::ggsave(out_png, p, width = width, height = height, dpi = dpi)
  }

  p
}