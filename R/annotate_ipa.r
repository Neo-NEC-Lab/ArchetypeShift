#' Annotate archetype signatures with IPA pathways and upstream regulators
#'
#' Annotates archetype-associated gene signatures by testing overlap or enrichment
#' against IPA canonical pathways and IPA upstream regulator target sets.
#'
#' The function accepts a named list of archetype gene signatures and two IPA-style
#' tables: one for canonical pathways and one for upstream regulators. Each IPA
#' term is converted into a gene set, optionally stratified by cluster, and then
#' compared against each archetype signature using either Fisher's exact test or
#' overlap-only summaries.
#'
#' @param signatures A non-empty named list of character vectors. Each element
#'   represents the gene signature for one archetype, and list names are treated
#'   as archetype identifiers.
#' @param pathways_df A data frame containing IPA canonical pathway results. Must
#'   include columns `pathway` and `molecules`. If `by_cluster = TRUE`, it must
#'   also contain a `cluster` column.
#' @param upstream_df A data frame containing IPA upstream regulator results. Must
#'   include columns `regulator` and `targets`. If `by_cluster = TRUE`, it must
#'   also contain a `cluster` column.
#' @param universe Optional character vector defining the background gene universe
#'   for enrichment testing. If `NULL`, the universe is constructed from the union
#'   of genes present in the signatures and IPA gene sets.
#' @param by_cluster Logical; if `TRUE`, pathway and regulator gene sets are built
#'   separately within each cluster. Default is `FALSE`.
#' @param method Enrichment method. `"fisher"` performs one-sided Fisher's exact
#'   tests, whereas `"overlap"` returns overlap summaries without significance
#'   testing. Default is `"fisher"`.
#' @param min_overlap Integer >= 0. Minimum number of overlapping genes required
#'   between an archetype signature and an IPA gene set to retain the result.
#'   Default is `1L`.
#' @param p_adjust Multiple-testing correction method passed to
#'   [stats::p.adjust()]. Default is `"BH"`.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `TRUE`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{pathways}{A `data.frame` of archetype enrichments against IPA canonical
#'   pathways.}
#'   \item{upstreams}{A `data.frame` of archetype enrichments against IPA upstream
#'   regulator target sets.}
#' }
#'
#' Returned tables include overlap counts, set sizes, expected overlap, ratios,
#' odds ratios, p-values, adjusted p-values, and retained IPA metadata such as
#' z-score, IPA p-value, and ratio where available.
#'
#' @details
#' The IPA input tables are converted into gene-set representations by parsing
#' comma-separated gene lists from `molecules` or `targets`. When `method = "fisher"`,
#' enrichment is tested using one-sided Fisher's exact tests against the supplied
#' or inferred gene universe. P-values are adjusted separately within IPA class
#' (and within cluster when `by_cluster = TRUE`).
#'
#' @examples
#' \dontrun{
#' ipa_res_k5 <- annotate_ipa(
#'   signatures = sig_k5$signatures,
#'   pathways_df = IPA_pathways,
#'   upstream_df = IPA_upstream,
#'   universe = NULL,
#'   by_cluster = TRUE,
#'   method = "fisher",
#'   min_overlap = 2L,
#'   p_adjust = "BH",
#'   verbose = TRUE
#' )
#'
#' head(ipa_res_k5$pathways)
#' head(ipa_res_k5$upstreams)
#' }
#'
#' @export
annotate_ipa <- function(
  signatures,
  pathways_df,
  upstream_df,
  universe = NULL,
  by_cluster = FALSE,
  method = c("fisher", "overlap")[1],
  min_overlap = 1L,
  p_adjust = "BH",
  verbose = TRUE
) {
  # -----------------------------
  # [CHECKS] signatures
  # -----------------------------
  if (is.null(signatures) || !is.list(signatures) || length(signatures) < 1L) {
    stop("`signatures` must be a non-empty named list of character vectors.", call. = FALSE)
  }
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must be a named list (names are archetype IDs).", call. = FALSE)
  }
  for (nm in names(signatures)) {
    if (!is.character(signatures[[nm]])) {
      stop("Each signatures[[archetype]] must be a character vector of genes.", call. = FALSE)
    }
  }

  # -----------------------------
  # [CHECKS] IPA inputs
  # -----------------------------
  if (!is.data.frame(pathways_df)) {
    stop("`pathways_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(upstream_df)) {
    stop("`upstream_df` must be a data.frame.", call. = FALSE)
  }

  req_p <- c("pathway", "molecules")
  req_u <- c("regulator", "targets")

  miss_p <- setdiff(req_p, colnames(pathways_df))
  miss_u <- setdiff(req_u, colnames(upstream_df))

  if (length(miss_p) > 0L) {
    stop("pathways_df missing columns: ", paste(miss_p, collapse = ", "), call. = FALSE)
  }
  if (length(miss_u) > 0L) {
    stop("upstream_df missing columns: ", paste(miss_u, collapse = ", "), call. = FALSE)
  }

  method <- match.arg(method, choices = c("fisher", "overlap"))

  if (!is.numeric(min_overlap) || length(min_overlap) != 1L || is.na(min_overlap) || min_overlap < 0) {
    stop("`min_overlap` must be a single integer >= 0.", call. = FALSE)
  }
  min_overlap <- as.integer(min_overlap)

  if (!is.logical(by_cluster) || length(by_cluster) != 1L) {
    stop("`by_cluster` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be TRUE/FALSE.", call. = FALSE)
  }

  # -----------------------------
  # [HELPER] parse comma-separated genes
  # -----------------------------
  parse_genes <- function(x) {
    if (length(x) == 0L) {
      return(character(0))
    }
    x <- x[!is.na(x)]
    if (length(x) == 0L) {
      return(character(0))
    }
    parts <- unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    unique(parts)
  }

  # -----------------------------
  # [HELPER] build IPA gene sets table
  # -----------------------------
  build_ipa_sets <- function(df, term_col, genes_col, ipa_class, by_cluster) {
    has_cluster <- "cluster" %in% colnames(df)

    if (by_cluster && !has_cluster) {
      stop("by_cluster=TRUE but IPA table has no 'cluster' column.", call. = FALSE)
    }

    key_cols <- character(0)
    if (by_cluster) {
      key_cols <- c(key_cols, "cluster")
    }
    key_cols <- c(key_cols, term_col)

    key_df <- df[, key_cols, drop = FALSE]
    for (cc in key_cols) {
      key_df[[cc]] <- as.character(key_df[[cc]])
    }

    key_string <- apply(key_df, 1, function(r) paste(r, collapse = "||"))
    groups <- split(seq_len(nrow(df)), key_string)

    out <- vector("list", length(groups))
    out_i <- 0L

    for (gk in names(groups)) {
      ii <- groups[[gk]]

      term_val <- unique(as.character(df[[term_col]][ii]))
      term_val <- term_val[!is.na(term_val) & nzchar(term_val)]
      if (length(term_val) == 0L) {
        next
      }
      term_val <- term_val[[1]]

      cluster_val <- NA_character_
      if (by_cluster) {
        cluster_val <- unique(as.character(df[["cluster"]][ii]))
        cluster_val <- cluster_val[!is.na(cluster_val) & nzchar(cluster_val)]
        cluster_val <- if (length(cluster_val) == 0L) NA_character_ else cluster_val[[1]]
      }

      geneset <- parse_genes(df[[genes_col]][ii])
      if (length(geneset) == 0L) {
        next
      }

      z_val <- NA_real_
      if ("zscore" %in% colnames(df)) {
        zz <- suppressWarnings(as.numeric(df[["zscore"]][ii]))
        zz <- zz[is.finite(zz)]
        if (length(zz) > 0L) {
          z_val <- zz[which.max(abs(zz))][[1]]
        }
      }

      p_val <- NA_real_
      if (ipa_class == "pathway" && "mlog10p" %in% colnames(df)) {
        mp <- suppressWarnings(as.numeric(df[["mlog10p"]][ii]))
        mp <- mp[is.finite(mp)]
        if (length(mp) > 0L) {
          p_val <- 10^(-max(mp))
        }
      }
      if (ipa_class == "upstream" && "pvalue" %in% colnames(df)) {
        pv <- suppressWarnings(as.numeric(df[["pvalue"]][ii]))
        pv <- pv[is.finite(pv)]
        if (length(pv) > 0L) {
          p_val <- min(pv)
        }
      }

      ratio_val <- NA_real_
      if (ipa_class == "pathway" && "ratio" %in% colnames(df)) {
        rr <- suppressWarnings(as.numeric(df[["ratio"]][ii]))
        rr <- rr[is.finite(rr)]
        if (length(rr) > 0L) {
          ratio_val <- max(rr)
        }
      }

      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        ipa_class = ipa_class,
        term = term_val,
        cluster = if (by_cluster) cluster_val else NA_character_,
        set_size = length(geneset),
        zscore = z_val,
        ipa_p = p_val,
        ratio = ratio_val,
        stringsAsFactors = FALSE
      )
      out[[out_i]]$geneset <- list(geneset)
    }

    out <- out[seq_len(out_i)]
    if (length(out) == 0L) {
      return(data.frame(
        ipa_class = character(0),
        term = character(0),
        cluster = character(0),
        set_size = integer(0),
        zscore = numeric(0),
        ipa_p = numeric(0),
        ratio = numeric(0),
        geneset = I(list()),
        stringsAsFactors = FALSE
      ))
    }

    do.call(rbind, out)
  }

  # -----------------------------
  # [BUILD] IPA gene sets
  # -----------------------------
  if (verbose) {
    message("Building IPA gene sets ...")
  }

  ipa_path_sets <- build_ipa_sets(
    pathways_df,
    term_col = "pathway",
    genes_col = "molecules",
    ipa_class = "pathway",
    by_cluster = by_cluster
  )
  ipa_up_sets <- build_ipa_sets(
    upstream_df,
    term_col = "regulator",
    genes_col = "targets",
    ipa_class = "upstream",
    by_cluster = by_cluster
  )

  # -----------------------------
  # [UNIVERSE] define background gene universe
  # -----------------------------
  if (is.null(universe)) {
    u1 <- unique(unlist(signatures, use.names = FALSE))
    u2 <- unique(unlist(ipa_path_sets$geneset, use.names = FALSE))
    u3 <- unique(unlist(ipa_up_sets$geneset, use.names = FALSE))
    universe <- unique(c(u1, u2, u3))
  }

  universe <- unique(trimws(as.character(universe)))
  universe <- universe[nzchar(universe)]
  if (length(universe) < 10L) {
    stop("`universe` is too small after cleaning.", call. = FALSE)
  }

  universe_size <- length(universe)

  # -----------------------------
  # [HELPER] enrichment for one IPA table
  # -----------------------------
  run_enrich <- function(ipa_tbl) {
    if (nrow(ipa_tbl) == 0L) {
      return(data.frame())
    }

    res <- list()
    res_i <- 0L

    for (a in names(signatures)) {
      a_genes <- unique(trimws(signatures[[a]]))
      a_genes <- a_genes[nzchar(a_genes)]
      a_genes <- intersect(a_genes, universe)
      a_size <- length(a_genes)
      if (a_size == 0L) {
        next
      }

      for (i in seq_len(nrow(ipa_tbl))) {
        term <- ipa_tbl$term[[i]]
        geneset <- ipa_tbl$geneset[[i]]
        geneset <- intersect(unique(trimws(geneset)), universe)
        set_size <- length(geneset)
        if (set_size == 0L) {
          next
        }

        overlap_genes <- intersect(a_genes, geneset)
        overlap <- length(overlap_genes)
        if (overlap < min_overlap) {
          next
        }

        expected <- (a_size * set_size) / universe_size

        a_only <- a_size - overlap
        set_only <- set_size - overlap
        neither <- universe_size - (overlap + a_only + set_only)

        if (min(c(overlap, a_only, set_only, neither)) < 0) {
          next
        }

        pval <- NA_real_
        stat <- NA_real_
        odds_ratio <- NA_real_

        if (method == "fisher") {
          mat <- matrix(c(overlap, a_only, set_only, neither), nrow = 2, byrow = TRUE)
          ft <- suppressWarnings(stats::fisher.test(mat, alternative = "greater"))
          pval <- as.numeric(ft$p.value)

          aa <- overlap + 0.5
          bb <- a_only + 0.5
          cc <- set_only + 0.5
          dd <- neither + 0.5
          odds_ratio <- (aa * dd) / (bb * cc)
          stat <- odds_ratio
        }

        gene_ratio <- overlap / a_size
        set_ratio <- overlap / set_size

        res_i <- res_i + 1L
        res[[res_i]] <- data.frame(
          ipa_class = ipa_tbl$ipa_class[[i]],
          term = term,
          cluster = ipa_tbl$cluster[[i]],
          archetype = a,
          overlap = overlap,
          a_size = a_size,
          set_size = set_size,
          universe_size = universe_size,
          expected = expected,
          gene_ratio = gene_ratio,
          set_ratio = set_ratio,
          odds_ratio = odds_ratio,
          log2_or = ifelse(is.finite(odds_ratio) && odds_ratio > 0, log2(odds_ratio), NA_real_),
          p_value = pval,
          zscore = ipa_tbl$zscore[[i]],
          ipa_p = ipa_tbl$ipa_p[[i]],
          ratio = ipa_tbl$ratio[[i]],
          stringsAsFactors = FALSE
        )
      }
    }

    if (res_i == 0L) {
      return(data.frame())
    }
    out <- do.call(rbind, res)

    if (method == "fisher" && nrow(out) > 0L) {
      if (by_cluster) {
        key <- paste(out$ipa_class, out$cluster, sep = "||")
      } else {
        key <- out$ipa_class
      }
      out$p_adj <- NA_real_
      for (kk in unique(key)) {
        ii <- which(key == kk)
        out$p_adj[ii] <- stats::p.adjust(out$p_value[ii], method = p_adjust)
      }
    } else {
      out$p_adj <- NA_real_
    }

    out$method <- method
    out$p_adjust <- p_adjust
    out
  }

  # -----------------------------
  # [RUN] enrichment
  # -----------------------------
  if (verbose) {
    message("Running enrichment (method=", method, ") ...")
  }
  pathways_res <- run_enrich(ipa_path_sets)
  upstream_res <- run_enrich(ipa_up_sets)

  # -----------------------------
  # [RETURN]
  # -----------------------------
  list(
    pathways = pathways_res,
    upstreams = upstream_res
  )
}
