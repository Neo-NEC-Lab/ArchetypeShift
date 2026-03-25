#' Plot archetype weight QC histograms
#'
#' Creates a histogram for a selected quality-control metric derived from
#' per-cell archetype weights, such as maximum archetype weight or normalized
#' entropy.
#'
#' This function is intended for use with QC summaries produced from archetype
#' weight matrices, for example by `compute_pcha_weight_qc_from_W()`.
#'
#' @param qc_df A data frame containing per-cell archetype QC metrics.
#' @param value_col Character scalar giving the column in `qc_df` to plot.
#' @param title Optional character scalar used as the plot title. If `NULL`,
#'   the title defaults to `"QC: <value_col>"`.
#' @param bins Integer giving the number of histogram bins. Default is `50`.
#'
#' @return A `ggplot2` histogram object.
#'
#' @examples
#' \dontrun{
#' qcA5 <- compute_pcha_weight_qc_from_W(W5)
#'
#' p_qc_max <- plot_pcha_qc_hist(
#'   qc_df = qcA5,
#'   value_col = "max_weight",
#'   title = "Archetype Max Weight Per Cell"
#' )
#'
#' p_qc_ent <- plot_pcha_qc_hist(
#'   qc_df = qcA5,
#'   value_col = "entropy_norm",
#'   title = "Archetype Entropy"
#' )
#' }
#'
#' @export
plot_pcha_qc_hist <- function(qc_df, value_col, title = NULL, bins = 50L) {
  # -----------------------------
  # [CHECKS] inputs
  # -----------------------------
  if (!is.data.frame(qc_df)) {
    stop("`qc_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(value_col) || length(value_col) != 1L || is.na(value_col) || !nzchar(value_col)) {
    stop("`value_col` must be a single non-empty character string.", call. = FALSE)
  }
  if (!(value_col %in% colnames(qc_df))) {
    stop("`value_col` not found in `qc_df`: ", value_col, call. = FALSE)
  }
  if (!is.null(title)) {
    if (!is.character(title) || length(title) != 1L || is.na(title) || !nzchar(title)) {
      stop("`title` must be NULL or a single non-empty character string.", call. = FALSE)
    }
  }
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 1) {
    stop("`bins` must be a single integer >= 1.", call. = FALSE)
  }
  bins <- as.integer(bins)

  if (is.null(title)) {
    title <- paste0("QC: ", value_col)
  }

  ggplot2::ggplot(qc_df, ggplot2::aes(x = .data[[value_col]])) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(title = title, x = value_col, y = "Cells") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}
