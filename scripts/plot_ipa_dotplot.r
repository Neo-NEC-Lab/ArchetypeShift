# =============================
# Plot IPA-style dotplots from annotate_ipa() results
# - x-axis: IPA z-score (direction/magnitude)
# - color:  -log10(p-value) from enrichment (IPA-like significance)
# - size:   overlap gene count (optional but useful)
# - works for both pathways and upstreams outputs from annotate_ipa()
# =============================

#' @title Plot IPA enrichment dotplot (z-score vs significance)
#'
#' @description
#' Creates an IPA-style dotplot for either canonical pathways or upstream regulators
#' from \code{annotate_ipa()} output. Encodes:
#'   - x: IPA z-score
#'   - color: -log10(p-value) from enrichment
#'   - size: overlap gene count
#'
#' @param enrich_df A data.frame from \code{annotate_ipa()} (e.g., \code{ipa_res$pathways} or \code{ipa_res$upstreams}).
#' @param archetype Character scalar archetype ID to plot (e.g., \code{"A5_A3"}).
#' @param cluster Optional character scalar. If not NULL and \code{enrich_df} contains a \code{cluster} column,
#'   filter to this cluster.
#' @param top_n Integer >= 1. Number of top terms to show after filtering.
#' @param min_overlap Integer >= 0. Minimum overlap to include.
#' @param p_max Optional numeric. If not NULL, keep only rows with \code{p_value <= p_max}.
#' @param use_padj Logical. If TRUE, color uses \code{p_adj} instead of \code{p_value}.
#' @param title Optional character scalar. If NULL, a default title is generated.
#' @param out_png Optional file path. If provided, saves plot as PNG.
#' @param width,height Numeric inches for PNG.
#' @param dpi Integer DPI for PNG.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' p_path_a2 <- plot_ipa_dotplot(
#'   enrich_df = ipa_res$pathways,
#'   archetype = "A5_A2",
#'   cluster = NULL,
#'   top_n = 20L,
#'   min_overlap = 3L,
#'   p_max = 0.05,
#'   use_padj = FALSE,
#'   order_by = "zscore",
#'   out_png = "A2_pathways_dotplot_ordered.png",
#'   width = 8,
#'   height = 5,
#'   dpi = 300,
#'   title = "A2 Canonical Pathways"
#' )
#' p_path_a2
#' }
#' @export
 plot_ipa_dotplot <- function(
    enrich_df,
    archetype,
    cluster = NULL,
    top_n = 20L,
    min_overlap = 2L,
    p_max = NULL,
    use_padj = FALSE,
    order_by = c("zscore", "neglog10p", "overlap")[1],
    title = NULL,
    out_png = NULL,
    width = 8,
    height = 5,
    dpi = 300
 ) {
   # -----------------------------
   # [CHECKS] inputs
   # -----------------------------
   if (!is.data.frame(enrich_df)) stop("`enrich_df` must be a data.frame.", call. = FALSE)
   if (!is.character(archetype) || length(archetype) != 1L || is.na(archetype) || !nzchar(archetype)) {
     stop("`archetype` must be a single non-empty character string.", call. = FALSE)
   }
   if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) stop("`top_n` must be >= 1.", call. = FALSE)
   if (!is.numeric(min_overlap) || length(min_overlap) != 1L || is.na(min_overlap) || min_overlap < 0) stop("`min_overlap` must be >= 0.", call. = FALSE)
   top_n <- as.integer(top_n)
   min_overlap <- as.integer(min_overlap)
   
   if (!is.null(cluster)) {
     if (!is.character(cluster) || length(cluster) != 1L || is.na(cluster) || !nzchar(cluster)) {
       stop("`cluster` must be NULL or a single non-empty character string.", call. = FALSE)
     }
     if (!("cluster" %in% colnames(enrich_df))) {
       stop("`cluster` was provided but enrich_df has no 'cluster' column.", call. = FALSE)
     }
   }
   
   if (!is.logical(use_padj) || length(use_padj) != 1L) stop("`use_padj` must be TRUE/FALSE.", call. = FALSE)
   order_by <- match.arg(order_by, choices = c("zscore", "neglog10p", "overlap"))
   
   req <- c("term", "archetype", "overlap", "p_value", "zscore")
   miss <- setdiff(req, colnames(enrich_df))
   if (length(miss) > 0L) stop("enrich_df missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
   if (use_padj && !("p_adj" %in% colnames(enrich_df))) stop("use_padj=TRUE requires column 'p_adj'.", call. = FALSE)
   
   # -----------------------------
   # [FILTER] archetype (+ cluster)
   # -----------------------------
   d <- enrich_df[enrich_df$archetype == archetype, , drop = FALSE]
   if (!is.null(cluster)) d <- d[d$cluster == cluster, , drop = FALSE]
   
   d$overlap <- suppressWarnings(as.integer(d$overlap))
   d$p_value <- suppressWarnings(as.numeric(d$p_value))
   d$zscore <- suppressWarnings(as.numeric(d$zscore))
   if (use_padj) d$p_adj <- suppressWarnings(as.numeric(d$p_adj))
   
   d <- d[is.finite(d$overlap) & d$overlap >= min_overlap, , drop = FALSE]
   d <- d[is.finite(d$p_value) & d$p_value > 0, , drop = FALSE]
   d <- d[is.finite(d$zscore), , drop = FALSE]
   
   if (!is.null(p_max)) {
     if (!is.numeric(p_max) || length(p_max) != 1L || is.na(p_max) || p_max <= 0) {
       stop("`p_max` must be NULL or a single numeric > 0.", call. = FALSE)
     }
     d <- d[d$p_value <= p_max, , drop = FALSE]
   }
   
   if (nrow(d) == 0L) stop("No rows remain after filtering.", call. = FALSE)
   
   # -----------------------------
   # [SELECT] top terms by significance (IPA-like)
   # -----------------------------
   d <- d[order(d$p_value, -abs(d$zscore), -d$overlap), , drop = FALSE]
   d <- d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
   
   # Color value
   pcol <- if (use_padj) d$p_adj else d$p_value
   pcol[pcol <= 0] <- NA_real_
   d$neglog10p <- -log10(pcol)
   
   # -----------------------------
   # [ORDER] terms by requested key (to get the "continuous line" look)
   # -----------------------------
   if (order_by == "zscore") {
     d <- d[order(d$zscore, decreasing = FALSE), , drop = FALSE]
   } else if (order_by == "neglog10p") {
     d <- d[order(d$neglog10p, decreasing = TRUE), , drop = FALSE]
   } else if (order_by == "overlap") {
     d <- d[order(d$overlap, decreasing = TRUE), , drop = FALSE]
   }
   
   d$term <- factor(as.character(d$term), levels = d$term)
   
   # -----------------------------
   # [TITLE]
   # -----------------------------
   if (is.null(title)) {
     base <- if ("ipa_class" %in% colnames(d)) unique(as.character(d$ipa_class))[1] else "IPA terms"
     if (!is.null(cluster)) {
       title <- paste0("Archetype ", archetype, ": ", base, " (", cluster, ")")
     } else {
       title <- paste0("Archetype ", archetype, ": ", base)
     }
   }
   
   # -----------------------------
   # [PLOT]
   # -----------------------------
   p <- ggplot2::ggplot(d, ggplot2::aes(x = zscore, y = term)) +
     ggplot2::geom_point(ggplot2::aes(size = overlap, color = neglog10p), alpha = 0.9) +
     ggplot2::labs(
       title = title,
       x = "IPA z-score",
       y = NULL,
       color = if (use_padj) "-log10(FDR)" else "-log10(p-value)",
       size = "Gene count"
     ) +
     ggplot2::scale_color_gradient2(
       low = "#4575B4",
       mid = "#7B3294",
       high = "#D73027",
       midpoint = median(d$neglog10p, na.rm = TRUE)
     ) +
     ggplot2::theme_classic() +
     ggplot2::theme(
       plot.title = ggplot2::element_text(face = "bold"),
       axis.text.y = ggplot2::element_text(size = 9),
       panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
       panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.18),
       panel.grid.minor = ggplot2::element_blank()
     )
   
   if (!is.null(out_png)) {
     ggplot2::ggsave(filename = out_png, plot = p, width = width, height = height, dpi = dpi, units = "in")
   }
   
   p
 }

# =============================
# Example usage (upstream regulators, program-level)
# - Uses color = -log10(p_value) and x = zscore
# =============================
# p_up_a3 <- plot_ipa_dotplot(
#   enrich_df = ipa_res_k5_global_paper$upstreams,
#   archetype = "A5_A3",
#   cluster = NULL,
#   top_n = 20L,
#   min_overlap = 3L,
#   p_max = 0.05,
#   use_padj = FALSE,
#   out_png = "A5_A3_upstream_dotplot.png",
#   width = 8,
#   height = 5,
#   dpi = 300
# )
# p_up_a3