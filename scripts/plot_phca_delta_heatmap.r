#' @title Plot Delta Occupancy Heatmap
#'
#' @description
#' Plots a heatmap of replicate-aware archetype occupancy differences (Δ occupancy)
#' across \emph{clusters × archetypes}. This is typically used to visualize
#' \code{group2 - group1} differences computed from sample-level occupancy
#' (e.g., output of \code{compute_pcha_delta_heatmap_df()}).
#'
#' @param df_delta A data.frame with columns \code{cluster}, \code{archetype}, and
#'   \code{delta}. Each row represents one \code{cluster × archetype} value.
#' @param cluster_order Optional character vector specifying the y-axis order for clusters.
#'   If provided, clusters not in \code{cluster_order} are dropped.
#' @param archetype_order Optional character vector specifying the x-axis order for archetypes.
#' @param title Character scalar plot title.
#' @param out_png Optional character path. If provided, saves the plot as a PNG.
#' @param width Numeric scalar. Width in inches for saved PNG.
#' @param height Numeric scalar. Height in inches for saved PNG.
#' @param dpi Integer scalar. Resolution for saved PNG.
#'
#' @return A ggplot object.
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
#' p
#' }
#'
#' @export
plot_pcha_delta_heatmap <- function(
  df_delta,
  cluster_order = NULL,
  archetype_order = NULL,
  title = "Δ occupancy (group2 − group1) by cluster and archetype",
  out_png = NULL,
  width = 7.8,
  height = 6.5,
  dpi = 300
) {
  # -----------------------------
  # [CHECKS] df structure
  # -----------------------------
  if (!is.data.frame(df_delta)) {
    stop("`df_delta` must be a data.frame.", call. = FALSE)
  }
  req_cols <- c("cluster", "archetype", "delta")
  missing_cols <- setdiff(req_cols, colnames(df_delta))
  if (length(missing_cols) > 0L) {
    stop("`df_delta` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Coerce delta to numeric if needed
  if (!is.numeric(df_delta$delta)) {
    suppressWarnings(xn <- as.numeric(df_delta$delta))
    if (all(is.na(xn))) stop("`df_delta$delta` is not numeric and cannot be coerced.", call. = FALSE)
    df_delta$delta <- xn
  }

  dfp <- df_delta

  # -----------------------------
  # [ORDER] cluster axis
  # -----------------------------
  if (!is.null(cluster_order)) {
    if (!is.character(cluster_order) || length(cluster_order) < 1L) {
      stop("`cluster_order` must be a non-empty character vector or NULL.", call. = FALSE)
    }
    dfp$cluster <- factor(dfp$cluster, levels = cluster_order)
    dfp <- dfp[!is.na(dfp$cluster), , drop = FALSE]
  }

  # -----------------------------
  # [ORDER] archetype axis
  # -----------------------------
  if (!is.null(archetype_order)) {
    if (!is.character(archetype_order) || length(archetype_order) < 1L) {
      stop("`archetype_order` must be a non-empty character vector or NULL.", call. = FALSE)
    }
    dfp$archetype <- factor(dfp$archetype, levels = archetype_order)
  }

  # -----------------------------
  # [PLOT]
  # -----------------------------
  p <- ggplot2::ggplot(dfp, ggplot2::aes(x = .data$archetype, y = .data$cluster, fill = .data$delta)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title = title,
      x = "Archetype program",
      y = "Cell-type annotation",
      fill = "Δ mean weight"
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