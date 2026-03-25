#' Compute correlation-based archetype gene signatures
#'
#' Computes archetype-associated gene signatures by correlating per-cell gene
#' expression with per-cell archetype weights. For each archetype, genes are
#' ranked by correlation coefficient and the top positively correlated genes are
#' returned as the signature.
#'
#' This function supports archetype weights supplied either as metadata columns
#' in a Seurat object or as an external `cells x k` matrix. Genes can be filtered
#' by minimum detection frequency and by default exclusion patterns
#' (for example mitochondrial, ribosomal, and highly abundant non-coding genes),
#' with optional user-specified regex-based exclusions.
#'
#' @param obj A Seurat object.
#' @param assay Character scalar giving the assay name. Default is `"RNA"`.
#' @param slot Character scalar giving the assay slot to use for expression
#'   values. Default is `"data"`.
#' @param weight_cols Optional character vector of metadata columns containing
#'   per-cell archetype weights, such as `c("A5_A1", "A5_A2", ...)`. If `NULL`,
#'   the function attempts to auto-detect archetype weight columns in
#'   `obj@meta.data`.
#' @param weights_mat Optional numeric matrix or data frame of archetype weights
#'   with dimensions `cells x k`. Row names must match Seurat cell IDs. If
#'   supplied, this is used instead of `weight_cols`.
#' @param method Correlation method, either `"spearman"` (default) or `"pearson"`.
#' @param top_n Integer giving the number of top positively correlated genes to
#'   retain per archetype. Default is `100L`.
#' @param min_detect_frac Numeric in `(0, 1]` giving the minimum fraction of cells
#'   in which a gene must be detected (`> 0`) to be retained. Default is `0.05`.
#' @param exclude_gene_default Logical; if `TRUE`, excludes genes matching a
#'   default set of patterns (`^MT-`, `^RPL`, `^RPS`, `^MALAT1$`, `^NEAT1$`)
#'   after detection filtering. Default is `TRUE`.
#' @param exclude_gene_regex Optional character vector of additional regex
#'   patterns used to exclude genes after detection filtering. Default is `NULL`.
#' @param max_genes Optional integer. If supplied, randomly subsamples the
#'   filtered gene set to this size for faster development runs. Default is `NULL`.
#' @param seed Integer random seed used only when `max_genes` is not `NULL`.
#'   Default is `1L`.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `TRUE`.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{ranked}{A `data.frame` with columns `gene`, `archetype`, `cor`,
#'   and `n_cells_used`, containing correlation-ranked genes for each archetype.}
#'   \item{signatures}{A named list of character vectors containing the top
#'   `top_n` positively correlated genes for each archetype.}
#'   \item{spec}{A list describing the inputs and filtering settings used,
#'   including assay, slot, method, filtering thresholds, archetype names,
#'   and matrix dimensions.}
#' }
#'
#' @details
#' Expression values are extracted from the requested Seurat assay and slot,
#' aligned to the cells for which archetype weights are available, filtered by
#' detection frequency, optionally filtered by gene exclusion patterns, and then
#' converted to a dense matrix for correlation calculation. Correlations are
#' computed independently for each archetype across all retained genes.
#'
#' @examples
#' \dontrun{
#' sig_k5 <- compute_signature_cor(
#'   obj = combined_no_IA,
#'   assay = "RNA",
#'   slot = "data",
#'   weight_cols = paste0("A5_A", 1:5),
#'   method = "spearman",
#'   top_n = 100L,
#'   min_detect_frac = 0.05,
#'   exclude_gene_default = TRUE,
#'   exclude_gene_regex = c("^HBA", "^HBB"),
#'   verbose = TRUE
#' )
#'
#' head(sig_k5$signatures[["A5_A5"]], 25)
#'
#' ranked_A5 <- sig_k5$ranked[sig_k5$ranked$archetype == "A5_A5", ]
#' ranked_A5 <- ranked_A5[order(ranked_A5$cor, decreasing = TRUE), ]
#' head(ranked_A5, 20)
#' }
#'
#' @export
compute_signature_cor <- function(
  obj,
  assay = "RNA",
  slot = "data",
  weight_cols = NULL,
  weights_mat = NULL,
  method = c("spearman", "pearson")[1],
  top_n = 100L,
  min_detect_frac = 0.05,
  exclude_gene_default = TRUE,
  exclude_gene_regex = NULL,
  max_genes = NULL,
  seed = 1L,
  verbose = TRUE
) {
  # -----------------------------
  # [CHECKS] object + assay/slot
  # -----------------------------
  if (!methods::is(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
    stop("`assay` must be a single non-empty character string.", call. = FALSE)
  }
  if (!assay %in% names(obj@assays)) {
    stop("Assay not found in Seurat object: ", assay, call. = FALSE)
  }
  if (!is.character(slot) || length(slot) != 1L || is.na(slot) || !nzchar(slot)) {
    stop("`slot` must be a single non-empty character string.", call. = FALSE)
  }
  method <- match.arg(method, choices = c("spearman", "pearson"))

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    stop("`top_n` must be a single integer >= 1.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(min_detect_frac) || length(min_detect_frac) != 1L || is.na(min_detect_frac) ||
      min_detect_frac <= 0 || min_detect_frac > 1) {
    stop("`min_detect_frac` must be a single numeric value in (0,1].", call. = FALSE)
  }

  if (!is.logical(exclude_gene_default) || length(exclude_gene_default) != 1L) {
    stop("`exclude_gene_default` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.null(exclude_gene_regex)) {
    if (!is.character(exclude_gene_regex) || length(exclude_gene_regex) < 1L ||
        any(is.na(exclude_gene_regex)) || any(!nzchar(exclude_gene_regex))) {
      stop("`exclude_gene_regex` must be NULL or a non-empty character vector of regex patterns.", call. = FALSE)
    }
  }

  if (!is.null(max_genes)) {
    if (!is.numeric(max_genes) || length(max_genes) != 1L || is.na(max_genes) || max_genes < 1) {
      stop("`max_genes` must be a single integer >= 1 or NULL.", call. = FALSE)
    }
    max_genes <- as.integer(max_genes)
  }
  seed <- as.integer(seed)

  # -----------------------------
  # [WEIGHTS] from weights_mat or weight_cols
  # -----------------------------
  md <- obj@meta.data
  cell_ids <- colnames(obj)

  if (!is.null(weights_mat)) {
    if (is.data.frame(weights_mat)) {
      weights_mat <- as.matrix(weights_mat)
    }
    if (!is.matrix(weights_mat) || !is.numeric(weights_mat)) {
      stop("`weights_mat` must be a numeric matrix/data.frame (cells x k).", call. = FALSE)
    }
    if (is.null(rownames(weights_mat)) || any(!nzchar(rownames(weights_mat)))) {
      stop("`weights_mat` must have non-empty rownames matching Seurat cell IDs.", call. = FALSE)
    }
    keep_cells <- intersect(cell_ids, rownames(weights_mat))
    if (length(keep_cells) == 0L) {
      stop("No overlapping cell IDs between Seurat object and `weights_mat` rownames.", call. = FALSE)
    }
    W <- weights_mat[keep_cells, , drop = FALSE]
    if (is.null(colnames(W)) || any(!nzchar(colnames(W)))) {
      colnames(W) <- paste0("A", seq_len(ncol(W)))
    }
    W <- as.matrix(W)
    use_cells <- rownames(W)
    arch_names <- colnames(W)
  } else {
    if (is.null(weight_cols)) {
      pat <- "^A\\d+_A\\d+$"
      weight_cols <- colnames(md)[grepl(pat, colnames(md))]
      if (length(weight_cols) == 0L) {
        stop("Could not auto-detect `weight_cols` in obj@meta.data. Pass `weight_cols` explicitly.", call. = FALSE)
      }
    }
    if (!is.character(weight_cols) || length(weight_cols) < 1L ||
        any(is.na(weight_cols)) || any(!nzchar(weight_cols))) {
      stop("`weight_cols` must be a non-empty character vector or NULL (auto-detect).", call. = FALSE)
    }
    missing_w <- setdiff(weight_cols, colnames(md))
    if (length(missing_w) > 0L) {
      stop("Missing `weight_cols` in obj@meta.data: ", paste(missing_w, collapse = ", "), call. = FALSE)
    }

    W <- md[, weight_cols, drop = FALSE]
    for (w in weight_cols) {
      if (!is.numeric(W[[w]])) {
        suppressWarnings(xn <- as.numeric(W[[w]]))
        if (all(is.na(xn))) {
          stop("Weight column '", w, "' cannot be coerced to numeric.", call. = FALSE)
        }
        W[[w]] <- xn
      }
    }
    W <- as.matrix(W)
    use_cells <- rownames(W)
    arch_names <- colnames(W)
  }

  # -----------------------------
  # [EXPR] get expression matrix aligned to use_cells
  # -----------------------------
  if (verbose) {
    message("Using assay=", assay, " slot=", slot, " method=", method, " top_n=", top_n)
    message("Cells used: ", length(use_cells), " | Archetypes: ", length(arch_names))
  }

  X <- Seurat::GetAssayData(obj, assay = assay, slot = slot)

  keep_cells2 <- intersect(use_cells, colnames(X))
  if (length(keep_cells2) == 0L) {
    stop("No overlapping cell IDs between weights and expression matrix.", call. = FALSE)
  }

  W <- W[keep_cells2, , drop = FALSE]
  X <- X[, keep_cells2, drop = FALSE]

  # -----------------------------
  # [FILTER] genes by detection fraction
  # -----------------------------
  det_frac <- Matrix::rowMeans(X > 0)
  keep_genes <- names(det_frac)[det_frac >= min_detect_frac]
  if (length(keep_genes) == 0L) {
    stop("No genes pass min_detect_frac = ", min_detect_frac, ".", call. = FALSE)
  }
  X <- X[keep_genes, , drop = FALSE]

  # -----------------------------
  # [FILTER] exclude unwanted gene patterns
  # -----------------------------
  exclude_default <- c(
    "^MT-",
    "^RPL",
    "^RPS",
    "^MALAT1$",
    "^NEAT1$"
  )

  patterns <- character(0)
  if (isTRUE(exclude_gene_default)) {
    patterns <- c(patterns, exclude_default)
  }
  if (!is.null(exclude_gene_regex)) {
    patterns <- c(patterns, exclude_gene_regex)
  }

  if (length(patterns) > 0L) {
    pat <- paste(patterns, collapse = "|")
    drop_genes <- grepl(pat, rownames(X), perl = TRUE)
    X <- X[!drop_genes, , drop = FALSE]
  }

  if (nrow(X) == 0L) {
    stop("All genes were filtered out after exclude_gene_default/exclude_gene_regex filtering.", call. = FALSE)
  }

  # -----------------------------
  # [OPTIONAL] subsample genes for speed
  # -----------------------------
  if (!is.null(max_genes) && nrow(X) > max_genes) {
    set.seed(seed)
    keep_sub <- sample(rownames(X), size = max_genes, replace = FALSE)
    X <- X[keep_sub, , drop = FALSE]
  }

  # -----------------------------
  # [COMPUTE] correlation per gene per archetype
  # -----------------------------
  ranked_list <- vector("list", length(arch_names))
  names(ranked_list) <- arch_names

  X_dense <- as.matrix(X)

  for (a in arch_names) {
    w <- as.numeric(W[, a])

    if (verbose) {
      message("Computing correlations for ", a, " ...")
    }

    cors <- apply(
      X_dense,
      1,
      function(g) {
        stats::cor(g, w, method = method, use = "pairwise.complete.obs")
      }
    )

    ranked_list[[a]] <- data.frame(
      gene = names(cors),
      archetype = a,
      cor = as.numeric(cors),
      n_cells_used = length(w),
      stringsAsFactors = FALSE
    )
  }

  ranked <- do.call(rbind, ranked_list)
  ranked <- ranked[is.finite(ranked$cor), , drop = FALSE]

  # -----------------------------
  # [SIGNATURES] top-N positive genes per archetype
  # -----------------------------
  signatures <- vector("list", length(arch_names))
  names(signatures) <- arch_names

  for (a in arch_names) {
    d <- ranked[ranked$archetype == a, , drop = FALSE]
    d <- d[order(d$cor, decreasing = TRUE), , drop = FALSE]
    d <- d[d$cor > 0, , drop = FALSE]
    signatures[[a]] <- head(d$gene, top_n)
  }

  # -----------------------------
  # [RETURN]
  # -----------------------------
  list(
    ranked = ranked,
    signatures = signatures,
    spec = list(
      assay = assay,
      slot = slot,
      method = method,
      top_n = top_n,
      min_detect_frac = min_detect_frac,
      exclude_gene_default = exclude_gene_default,
      exclude_gene_regex = exclude_gene_regex,
      max_genes = max_genes,
      seed = seed,
      archetypes = arch_names,
      n_cells = ncol(X_dense),
      n_genes = nrow(X_dense)
    )
  )
}
