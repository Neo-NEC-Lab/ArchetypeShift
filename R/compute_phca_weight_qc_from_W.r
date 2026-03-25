#' Compute per-cell QC metrics from an archetype weight matrix
#'
#' Computes simple quality-control summaries from a per-cell archetype weight
#' matrix, including row sums, maximum archetype weight, normalized entropy,
#' and near-zero sparsity.
#'
#' This helper is intended for inspecting the sharpness and stability of
#' archetype assignments after PCHA fitting. The returned data frame can be used
#' directly with [plot_pcha_qc_hist()] to visualize per-cell QC distributions.
#'
#' @param W Numeric matrix of archetype weights with dimensions `cells x k`.
#'   Rows should correspond to cells and columns to archetypes.
#' @param near_zero Numeric threshold used to define near-zero weights when
#'   computing `sparsity_near_zero`. Default is `1e-3`.
#'
#' @return A `data.frame` with one row per cell and the following columns:
#' \describe{
#'   \item{cell_id}{Cell identifier taken from `rownames(W)`.}
#'   \item{row_sum}{Sum of archetype weights across each cell.}
#'   \item{max_weight}{Maximum archetype weight for each cell.}
#'   \item{entropy_norm}{Normalized entropy of the archetype weight vector for
#'   each cell, scaled to the range `[0, 1]`. Lower values indicate more
#'   concentrated archetype assignments.}
#'   \item{sparsity_near_zero}{Fraction of archetype weights below
#'   `near_zero` for each cell.}
#' }
#'
#' @details
#' Negative weights, if present, are truncated to zero before entropy is
#' computed. Entropy is calculated from row-normalized non-negative weights
#' using a small numerical offset to avoid taking the logarithm of zero.
#'
#' @examples
#' \dontrun{
#' qcA5 <- compute_pcha_weight_qc_from_W(W5)
#'
#' head(qcA5)
#'
#' plot_pcha_qc_hist(
#'   qc_df = qcA5,
#'   value_col = "max_weight",
#'   title = "Archetype Max Weight Per Cell"
#' )
#' }
#'
#' @export
compute_pcha_weight_qc_from_W <- function(W, near_zero = 1e-3) {
  # -----------------------------
  # [CHECKS] inputs
  # -----------------------------
  if (!is.matrix(W)) {
    stop("`W` must be a matrix with dimensions cells x k.", call. = FALSE)
  }
  if (!is.numeric(W)) {
    stop("`W` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(W) < 1L || ncol(W) < 2L) {
    stop("`W` must have at least 1 row and 2 columns.", call. = FALSE)
  }
  if (is.null(rownames(W)) || any(!nzchar(rownames(W)))) {
    stop("`W` must have non-empty rownames corresponding to cell IDs.", call. = FALSE)
  }
  if (!is.numeric(near_zero) || length(near_zero) != 1L || is.na(near_zero) || near_zero < 0) {
    stop("`near_zero` must be a single numeric value >= 0.", call. = FALSE)
  }
  if (anyNA(W)) {
    stop("`W` must not contain NA values.", call. = FALSE)
  }

  # -----------------------------
  # [COMPUTE] QC summaries
  # -----------------------------
  rs <- rowSums(W)
  mx <- apply(W, 1L, max)

  eps <- 1e-12
  W_pos <- pmax(W, 0)
  W_norm <- W_pos / pmax(rowSums(W_pos), eps)
  ent <- -rowSums(W_norm * log(pmax(W_norm, eps)))
  ent_norm <- ent / log(ncol(W))
  spars <- rowMeans(W < near_zero)

  data.frame(
    cell_id = rownames(W),
    row_sum = rs,
    max_weight = mx,
    entropy_norm = ent_norm,
    sparsity_near_zero = spars,
    stringsAsFactors = FALSE
  )
}
 

