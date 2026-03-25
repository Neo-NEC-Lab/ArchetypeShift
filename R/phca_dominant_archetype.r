#' @title Assign dominant archetype per cell
#'
#' @description
#' Given a \code{cells x k} archetype weight matrix, assigns each cell to its
#' dominant archetype (argmax) and computes simple confidence metrics:
#' max weight, second-highest weight, and the margin between them.
#' Optionally marks low-confidence cells as \code{"Unassigned"} using thresholds.
#'
#' @param W Numeric matrix of archetype weights with dimensions \code{cells x k}.
#'   Rows should correspond to cells and columns to archetypes.
#' @param prefix Character scalar prefix used to generate archetype labels
#'   (default \code{"A"}). If \code{colnames(W)} already exist, they are used directly.
#' @param min_max_weight Optional numeric scalar. If provided, cells with max weight
#'   < \code{min_max_weight} are labeled \code{"Unassigned"}.
#' @param min_margin Optional numeric scalar. If provided, cells with
#'   (max - second) < \code{min_margin} are labeled \code{"Unassigned"}.
#' @param unassigned_label Character scalar label for low-confidence assignments
#'   (default \code{"Unassigned"}).
#'
#' @return A data.frame with one row per cell and columns:
#'   \code{cell_id}, \code{dominant}, \code{max_weight}, \code{second_weight}, \code{margin}.
#'
#' @examples
#' \dontrun{
#' W <- pcha_extract_weights(fit_k5)
#' dom <- pcha_dominant_archetype(W, min_max_weight = 0.4, min_margin = 0.1)
#' table(dom$dominant)
#' }
#'
#' @export
pcha_dominant_archetype <- function(
  W,
  prefix = "A",
  min_max_weight = NULL,
  min_margin = NULL,
  unassigned_label = "Unassigned"
) {
  # -----------------------------
  # Validate inputs
  # -----------------------------
  if (!is.matrix(W)) stop("W must be a matrix (cells x k).", call. = FALSE)
  if (nrow(W) < 1L || ncol(W) < 2L) stop("W must have dimensions (cells >= 1) x (k >= 2).", call. = FALSE)
  if (anyNA(W)) stop("W contains NA values.", call. = FALSE)
  if (!is.numeric(W)) stop("W must be numeric.", call. = FALSE)

  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) || !nzchar(prefix)) {
    stop("prefix must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(unassigned_label) || length(unassigned_label) != 1L || is.na(unassigned_label) || !nzchar(unassigned_label)) {
    stop("unassigned_label must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.null(min_max_weight) && (!is.numeric(min_max_weight) || length(min_max_weight) != 1L || is.na(min_max_weight))) {
    stop("min_max_weight must be a numeric scalar or NULL.", call. = FALSE)
  }
  if (!is.null(min_margin) && (!is.numeric(min_margin) || length(min_margin) != 1L || is.na(min_margin))) {
    stop("min_margin must be a numeric scalar or NULL.", call. = FALSE)
  }

  # -----------------------------
  # Resolve archetype labels
  # -----------------------------
  arch_names <- colnames(W)
  if (is.null(arch_names)) arch_names <- paste0(prefix, seq_len(ncol(W)))

  # -----------------------------
  # Dominant and second-best weights
  # -----------------------------
  max_weight <- apply(W, 1, max)
  dom_i <- apply(W, 1, which.max)

  second_weight <- vapply(seq_len(nrow(W)), function(i) {
    w <- W[i, ]
    w[dom_i[i]] <- -Inf
    max(w)
  }, numeric(1))

  margin <- max_weight - second_weight

  dominant <- arch_names[dom_i]

  # -----------------------------
  # Optional low-confidence filtering
  # -----------------------------
  low_conf <- rep(FALSE, length(dominant))
  if (!is.null(min_max_weight)) low_conf <- low_conf | (max_weight < min_max_weight)
  if (!is.null(min_margin)) low_conf <- low_conf | (margin < min_margin)
  dominant[low_conf] <- unassigned_label

  # -----------------------------
  # Build output table
  # -----------------------------
  cell_id <- rownames(W)
  if (is.null(cell_id)) cell_id <- paste0("cell_", seq_len(nrow(W)))

  data.frame(
    cell_id = cell_id,
    dominant = dominant,
    max_weight = as.numeric(max_weight),
    second_weight = as.numeric(second_weight),
    margin = as.numeric(margin),
    stringsAsFactors = FALSE
  )
}
