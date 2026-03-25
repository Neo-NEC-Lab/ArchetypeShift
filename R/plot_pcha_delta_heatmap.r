#' Plot a delta occupancy heatmap
#'
#' Creates a heatmap of replicate-aware archetype occupancy differences across
#' clusters and archetypes. This is typically used to visualize `group2 - group1`
#' differences computed from sample-level occupancy summaries, for example from
#' [compute_pcha_delta_heatmap_df()].
#'
#' @param df_delta A data frame with columns `cluster`, `archetype`, and `delta`.
#'   Each row represents one `cluster x archetype` occupancy difference.
#' @param cluster_order Optional character vector specifying the display order of
#'   clusters on the y-axis. Clusters not present in `cluster_order` are dropped.
#' @param archetype_order Optional character vector specifying the display order
#'   of archetypes on the x-axis.
#' @param title Character scalar used as the plot title. Defaults to
#'   `"Δ occupancy (group2 − group1) by cluster and archetype"`.
#' @param out_png Optional character path. If supplied, the plot is saved as a
#'   PNG file.
#' @param width Numeric scalar giving the saved PNG width in inches. Default is
#'   `7.8`.
#' @param height Numeric scalar giving the saved PNG height in inches. Default is
#'   `6.5`.
#' @param dpi Integer scalar giving the saved PNG resolution. Default is `300`.
#'
#' @return A `ggplot2` heatmap object.
#'
#' @examples
#' \dontrun{
#' p <- plot_pcha_delta_heatmap(
#'   df_delta = df_delta,
#'   cluster_order = cluster_order,
#'   archetype_order = paste0("A", 1:5),
#'   title = "Δ occupancy (NEC_HA − NEC) by cluster and archetype",
#'   out_png = "delta_heatmap.png"
#' )
#'
#' p
#' }
#'
#' @export
plot_pcha_delta_heatmap <- function(
  df_delta,
  cluster_order = NULL,
  archetype_order = NULL,
  title = "\u0394 occupancy (group2 \u2212 group1) by cluster and archetype",
  out_png = NULL,
  width = 7.8,
  height = 6.5,
  dpi = 300
) {
  # -----------------------------
  # [CHECKS] data frame structure
  # -----------------------------
  if (!is.data.frame(df_delta)) {
    stop("`df_delta` must be a data.frame.", call. = FALSE)
  }

  req_cols <- c("cluster", "archetype", "delta")
  missing_cols <- setdiff(req_cols, colnames(df_delta))
  if (length(missing_cols) > 0L) {
    stop("`df_delta` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!is.null(cluster_order)) {
    if (!is.character(cluster_order) || length(cluster_order) < 1L || any(is.na(cluster_order)) || any(!nzchar(cluster_order))) {
      stop("`cluster_order` must be NULL or a non-empty character vector.", call. = FALSE)
    }
  }

  if (!is.null(archetype_order)) {
    if (!is.character(archetype_order) || length(archetype_order) < 1L || any(is.na(archetype_order)) || any(!nzchar(archetype_order))) {
      stop("`archetype_order` must be NULL or a non-empty character vector.", call. = FALSE)
    }
  }

  if (!is.character(title) || length(title) != 1L || is.na(title) || !nzchar(title)) {
    stop("`title` must be a single non-empty character string.", call. = FALSE)
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
  # [COERCE] delta column to numeric
  # -----------------------------
  if (!is.numeric(df_delta$delta)) {
    suppressWarnings(xn <- as.numeric(df_delta$delta))
    if (all(is.na(xn))) {
      stop("`df_delta$delta` is not numeric and cannot be coerced.", call. = FALSE)
    }
    df_delta$delta <- xn
  }

  dfp <- df_delta

  # -----------------------------
  # [ORDER] cluster axis
  # -----------------------------
  if (!is.null(cluster_order)) {
    dfp$cluster <- factor(dfp$cluster, levels = cluster_order)
    dfp <- dfp[!is.na(dfp$cluster), , drop = FALSE]
    if (nrow(dfp) == 0L) {
      stop("No rows remain after applying `cluster_order`.", call. = FALSE)
    }
  }

  # -----------------------------
  # [ORDER] archetype axis
  # -----------------------------
  if (!is.null(archetype_order)) {
    dfp$archetype <- factor(dfp$archetype, levels = archetype_order)
  }

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(
    dfp,
    ggplot2::aes(
      x = rlang::.data$archetype,
      y = rlang::.data$cluster,
      fill = rlang::.data$delta
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title = title,
      x = "Archetype program",
      y = "Cell-type annotation",
      fill = "\u0394 mean weight"
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
