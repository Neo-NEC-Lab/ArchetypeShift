#' @title Extract robust k-scan model-selection metrics from ParetoTI output
#'
#' @description
#' Extracts a standardized model-selection metrics table from a ParetoTI
#' \code{k_pch_fit} object (as returned by \code{ParetoTI::k_fit_pch()}).
#' This function is robust to different internal representations:
#' (1) \code{k_scan$pch_fits} may be a list of per-k \code{pch_fit} objects, or
#' (2) \code{k_scan$pch_fits} may contain metric vectors. Metrics from
#' \code{k_scan$summary} (if present) are merged in to fill missing values.
#'
#' @param k_scan A ParetoTI \code{k_pch_fit} object (list-like) containing at least
#'   \code{$pch_fits} and optionally \code{$summary}.
#'
#' @return A data.frame with one row per \code{k} and columns (when available):
#'   \code{k}, \code{SSE}, \code{varexpl}, \code{t_ratio}, \code{arc_vol}, \code{hull_vol},
#'   \code{total_var}, and \code{var_explained} (canonical variance explained).
#'
#' @examples
#' \dontrun{
#' metrics_tbl <- pcha_kscan_metrics(k_scan)
#' print(metrics_tbl)
#' }
#'
#' @export
pcha_kscan_metrics <- function(k_scan) {
  # -----------------------------
  # Inputs
  # -----------------------------
  if (is.null(k_scan) || !is.list(k_scan)) stop("k_scan must be list-like.", call. = FALSE)
  if (is.null(k_scan$pch_fits) || !is.list(k_scan$pch_fits)) stop("k_scan$pch_fits must be a list.", call. = FALSE)

  pf <- k_scan$pch_fits
  pf_names <- names(pf)
  if (is.null(pf_names)) pf_names <- rep("", length(pf))

  # -----------------------------
  # Helper: identify fit-like entries
  # -----------------------------
  is_fit_like <- function(x) {
    is.list(x) && any(c("SSE", "varexpl", "t_ratio", "arc_vol", "hull_vol", "S", "XC") %in% names(x))
  }

  fit_idx <- which(vapply(pf, is_fit_like, logical(1)))

  # -----------------------------
  # Case A: list of pch_fit objects
  # -----------------------------
  if (length(fit_idx) > 0L) {
    rows <- lapply(fit_idx, function(i) {
      fit <- pf[[i]]

      # k: prefer nrow(S), else parse from element name
      k_val <- NA_integer_
      if (!is.null(fit$S) && is.matrix(fit$S)) {
        k_val <- as.integer(nrow(fit$S))
      } else if (nzchar(pf_names[i])) {
        k_guess <- suppressWarnings(as.integer(gsub("[^0-9]+", "", pf_names[i])))
        if (is.finite(k_guess)) k_val <- as.integer(k_guess)
      }

      data.frame(
        k = as.integer(k_val),
        SSE = if (!is.null(fit$SSE)) as.numeric(fit$SSE) else NA_real_,
        varexpl = if (!is.null(fit$varexpl)) as.numeric(fit$varexpl) else NA_real_,
        t_ratio = if (!is.null(fit$t_ratio)) as.numeric(fit$t_ratio) else NA_real_,
        arc_vol = if (!is.null(fit$arc_vol)) as.numeric(fit$arc_vol) else NA_real_,
        hull_vol = if (!is.null(fit$hull_vol)) as.numeric(fit$hull_vol) else NA_real_,
        total_var = if (!is.null(fit$total_var)) as.numeric(fit$total_var) else NA_real_,
        stringsAsFactors = FALSE
      )
    })

    metrics <- do.call(rbind, rows)

  } else {
    # -----------------------------
    # Case B: metric vectors stored directly under pch_fits
    # -----------------------------
    # infer k
    k_vec <- NULL
    if (!is.null(pf$SSE) && !is.null(names(pf$SSE)) && all(nzchar(names(pf$SSE)))) {
      k_vec <- suppressWarnings(as.integer(names(pf$SSE)))
    } else if (!is.null(pf$varexpl)) {
      k_vec <- seq_along(pf$varexpl)
    } else if (!is.null(k_scan$summary) && ("k" %in% colnames(k_scan$summary))) {
      k_vec <- as.integer(k_scan$summary$k)
    } else {
      stop("Cannot infer k values from k_scan$pch_fits.", call. = FALSE)
    }

    metrics <- data.frame(k = as.integer(k_vec), stringsAsFactors = FALSE)

    add_vec <- function(df, colname, vec) {
      if (is.null(vec)) return(df)
      v <- suppressWarnings(as.numeric(vec))

      if (!is.null(names(vec)) && all(nzchar(names(vec)))) {
        kk <- suppressWarnings(as.integer(names(vec)))
        m <- match(df$k, kk)
        df[[colname]] <- v[m]
      } else {
        df[[colname]] <- v
      }
      df
    }

    metrics <- add_vec(metrics, "SSE", pf$SSE)
    metrics <- add_vec(metrics, "arc_vol", pf$arc_vol)
    metrics <- add_vec(metrics, "hull_vol", pf$hull_vol)
    metrics <- add_vec(metrics, "t_ratio", pf$t_ratio)
    metrics <- add_vec(metrics, "varexpl", pf$varexpl)
    metrics <- add_vec(metrics, "total_var", pf$total_var)
  }

  # -----------------------------
  # Merge in k_scan$summary where present (fills missing varexpl/t_ratio/total_var)
  # -----------------------------
  if (!is.null(k_scan$summary)) {
    summ <- as.data.frame(k_scan$summary, stringsAsFactors = FALSE)
    if ("k" %in% colnames(summ)) {
      summ$k <- as.integer(summ$k)

      metrics <- merge(metrics, summ, by = "k", all.x = TRUE, suffixes = c("", "_summary"))

      fill_from <- function(col) {
        col_s <- paste0(col, "_summary")
        if (!(col_s %in% colnames(metrics))) return()
        if (!(col %in% colnames(metrics))) metrics[[col]] <- NA_real_
        cur <- suppressWarnings(as.numeric(metrics[[col]]))
        add <- suppressWarnings(as.numeric(metrics[[col_s]]))
        to_fill <- !is.finite(cur) & is.finite(add)
        cur[to_fill] <- add[to_fill]
        metrics[[col]] <- cur
      }

      fill_from("varexpl")
      fill_from("t_ratio")
      fill_from("total_var")

      drop_cols <- grep("_summary$", colnames(metrics), value = TRUE)
      if (length(drop_cols) > 0) metrics <- metrics[, setdiff(colnames(metrics), drop_cols), drop = FALSE]
    }
  }

  # -----------------------------
  # Canonical var explained
  # -----------------------------
  tv <- if ("total_var" %in% colnames(metrics)) suppressWarnings(as.numeric(metrics$total_var)) else rep(NA_real_, nrow(metrics))
  vx <- if ("varexpl" %in% colnames(metrics)) suppressWarnings(as.numeric(metrics$varexpl)) else rep(NA_real_, nrow(metrics))
  metrics$var_explained <- ifelse(is.finite(tv), tv, vx)

  metrics <- metrics[order(metrics$k), , drop = FALSE]
  rownames(metrics) <- NULL
  metrics
}
