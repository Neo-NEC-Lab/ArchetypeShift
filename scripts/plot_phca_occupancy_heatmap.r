#' @title Plot Archetype Occupancy Heatmap
#'
#' @description
#' Plots a heatmap of archetype program occupancy across a grouping variable
#' (e.g., cluster/cell-type annotation or condition). Input should be an occupancy
#' table such as the output of \code{compute_pcha_occupancy()}.
#'
#' @param occ_df A data.frame containing the grouping column and archetype occupancy columns.
#' @param group_col Character scalar. Column in \code{occ_df} to use on the y-axis
#'   (e.g., \code{"merged_cluster_annotations"}).
#' @param weight_cols Character vector of occupancy columns in \code{occ_df}
#'   (e.g., \code{c("A5_A1","A5_A2",...)}).
#' @param title Optional character scalar plot title. If NULL, a default title is generated.
#' @param drop_weight_prefix_regex Optional character scalar regex used to clean archetype labels
#'   for display (e.g., \code{"^A5_"}). Use NULL to disable.
#' @param group_order Optional character vector specifying the order of \code{group_col} levels
#'   on the y-axis. Rows not in \code{group_order} are dropped.
#' @param out_png Optional character path. If provided, saves the plot as a PNG.
#' @param width Numeric scalar. Width in inches for saved PNG.
#' @param height Numeric scalar. Height in inches for saved PNG.
#' @param dpi Integer scalar. Resolution for saved PNG.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' occ_by_cluster <- compute_pcha_occupancy(
#'   obj = combined_no_IA,
#'   weight_cols = paste0("A5_A", 1:5),
#'   group_cols = c("merged_cluster_annotations")
#' )
#'
#' p <- plot_pcha_occupancy_heatmap(
#'   occ_df = occ_by_cluster,
#'   group_col = "merged_cluster_annotations",
#'   weight_cols = paste0("A5_A", 1:5),
#'   drop_weight_prefix_regex = "^A5_",
#'   group_order = cluster_order,
#'   title = "A5 program occupancy by cell-type annotation"
#' )
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
    stop("`group_col` not found in occ_df: ", group_col, call. = FALSE)
  }

  if (!is.character(weight_cols) || length(weight_cols) < 1L) {
    stop("`weight_cols` must be a non-empty character vector.", call. = FALSE)
  }
  weight_cols <- unique(weight_cols)

  missing_w <- setdiff(weight_cols, colnames(occ_df))
  if (length(missing_w) > 0L) {
    stop("Missing `weight_cols` in occ_df: ", paste(missing_w, collapse = ", "), call. = FALSE)
  }

  if (!is.null(drop_weight_prefix_regex)) {
    if (!is.character(drop_weight_prefix_regex) || length(drop_weight_prefix_regex) != 1L) {
      stop("`drop_weight_prefix_regex` must be a single regex string or NULL.", call. = FALSE)
    }
  }

  if (!is.null(group_order)) {
    if (!is.character(group_order) || length(group_order) < 1L) {
      stop("`group_order` must be a non-empty character vector or NULL.", call. = FALSE)
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
  df_long <- occ_df %>%
    dplyr::select(dplyr::all_of(c(group_col, weight_cols))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(weight_cols), names_to = "archetype", values_to = "occupancy")

  # Clean archetype labels + set order
  if (is.null(drop_weight_prefix_regex)) {
    archetype_levels <- weight_cols
    df_long$archetype <- df_long$archetype
  } else {
    archetype_levels <- gsub(drop_weight_prefix_regex, "", weight_cols)
    df_long$archetype <- gsub(drop_weight_prefix_regex, "", df_long$archetype)
  }
  df_long$archetype <- factor(df_long$archetype, levels = archetype_levels)

  # Order group axis if requested
  if (!is.null(group_order)) {
    df_long[[group_col]] <- factor(df_long[[group_col]], levels = group_order)
    df_long <- df_long[!is.na(df_long[[group_col]]), , drop = FALSE]
  }

  if (is.null(title)) title <- paste0("Archetype occupancy by ", group_col)

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$archetype, y = .data[[group_col]], fill = .data$occupancy)) +
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