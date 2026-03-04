#' @title Extract per-cell archetype weights from a ParetoTI fit
#'
#' @description
#' Extracts the archetype weight matrix from a ParetoTI \code{pch_fit} object and
#' returns a standardized \code{cells x k} matrix (\code{W}) with consistent
#' row/column names. ParetoTI stores weights as \code{fit$S} in \code{k x cells}
#' orientation, so this function transposes to \code{cells x k}.
#'
#' @param fit A ParetoTI \code{pch_fit} object (or list-like object) containing \code{fit$S}.
#' @param cell_ids Optional character vector of length \code{ncol(fit$S)} providing cell IDs
#'   to use as rownames of the returned matrix. If \code{NULL}, \code{colnames(fit$S)} are used.
#' @param prefix Character scalar prefix for archetype column names (default \code{"A"}).
#'
#' @return A numeric matrix of archetype weights with dimensions \code{cells x k}.
#'   Row names are cell IDs and column names are \code{paste0(prefix, 1:k)}.
#'
#' @examples
#' \dontrun{
#' W <- pcha_extract_weights(fit = fit_k5)
#' dim(W)  # cells x k
#' head(colnames(W))
#' }
#'
#' @export
pcha_extract_weights <- function(fit, cell_ids = NULL, prefix = "A") {
  # -----------------------------
  # Validate inputs
  # -----------------------------
  if (is.null(fit) || !is.list(fit)) {
    stop("fit must be a non-null list-like object (e.g., a ParetoTI pch_fit).", call. = FALSE)
  }
  if (!("S" %in% names(fit))) {
    stop("fit must contain element $S (k x cells).", call. = FALSE)
  }
  S <- fit$S
  if (!is.matrix(S)) stop("fit$S must be a matrix.", call. = FALSE)
  if (nrow(S) < 2L || ncol(S) < 1L) stop("fit$S must have dimensions (k >= 2) x (cells >= 1).", call. = FALSE)
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) || !nzchar(prefix)) {
    stop("prefix must be a non-empty character scalar.", call. = FALSE)
  }

  # -----------------------------
  # Resolve cell IDs
  # -----------------------------
  if (is.null(cell_ids)) {
    if (is.null(colnames(S))) {
      stop(
        "cell_ids is NULL and colnames(fit$S) is NULL. Provide cell_ids to name/align cells.",
        call. = FALSE
      )
    }
    cell_ids_use <- colnames(S)
  } else {
    if (!is.character(cell_ids)) stop("cell_ids must be a character vector.", call. = FALSE)
    if (length(cell_ids) != ncol(S)) stop("cell_ids must have length ncol(fit$S).", call. = FALSE)
    cell_ids_use <- cell_ids
  }

  # -----------------------------
  # Standardize orientation: cells x k
  # -----------------------------
  W <- t(S)
  rownames(W) <- cell_ids_use
  colnames(W) <- paste0(prefix, seq_len(ncol(W)))

  W
}