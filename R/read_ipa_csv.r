#' Read an IPA export CSV with simple encoding fallback
#'
#' Reads a CSV file exported from Ingenuity Pathway Analysis (IPA), with a
#' fallback encoding strategy for files that may contain non-UTF-8 characters.
#'
#' This helper is intended for canonical pathway and upstream regulator export
#' tables that are later passed to [annotate_ipa()]. The function first attempts
#' a standard CSV import and, if that fails, retries using a fallback file
#' encoding.
#'
#' @param path Character scalar giving the path to the IPA CSV file.
#' @param strings_as_factors Logical; passed to [utils::read.csv()]. Default is
#'   `FALSE`.
#' @param check_names Logical; passed to [utils::read.csv()] as `check.names`.
#'   Default is `FALSE`.
#' @param fallback_encoding Character scalar giving the fallback encoding used if
#'   the initial read fails. Default is `"latin1"`.
#' @param verbose Logical; if `TRUE`, prints a message when the fallback encoding
#'   is used. Default is `TRUE`.
#'
#' @return A `data.frame` containing the imported IPA table.
#'
#' @details
#' IPA exports are sometimes saved with encodings that can cause a standard
#' `read.csv()` call to fail or import text incorrectly. This helper provides a
#' lightweight fallback mechanism without imposing additional package dependencies.
#'
#' @examples
#' \dontrun{
#' pathways_df <- read_ipa_csv("IPA_pathways.csv")
#' upstream_df <- read_ipa_csv("IPA_upstream.csv")
#' }
#'
#' @export
read_ipa_csv <- function(
  path,
  strings_as_factors = FALSE,
  check_names = FALSE,
  fallback_encoding = "latin1",
  verbose = TRUE
) {
  # -----------------------------
  # [CHECKS] input path and options
  # -----------------------------
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  if (!is.logical(strings_as_factors) || length(strings_as_factors) != 1L) {
    stop("`strings_as_factors` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(check_names) || length(check_names) != 1L) {
    stop("`check_names` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.character(fallback_encoding) || length(fallback_encoding) != 1L ||
      is.na(fallback_encoding) || !nzchar(fallback_encoding)) {
    stop("`fallback_encoding` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be TRUE/FALSE.", call. = FALSE)
  }

  # -----------------------------
  # [READ] standard attempt, then fallback encoding
  # -----------------------------
  out <- tryCatch(
    utils::read.csv(
      file = path,
      stringsAsFactors = strings_as_factors,
      check.names = check_names
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    if (verbose) {
      message("Standard CSV read failed; retrying with fileEncoding='", fallback_encoding, "'.")
    }

    out <- utils::read.csv(
      file = path,
      stringsAsFactors = strings_as_factors,
      check.names = check_names,
      fileEncoding = fallback_encoding
    )
  }

  out
}
