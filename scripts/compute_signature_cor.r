#' @title Compute correlation-based archetype gene signatures
#'
#' @description
#' Computes archetype gene signatures by correlating per-cell gene expression with
#' per-cell archetype weights. For each archetype, genes are ranked by correlation
#' (default Spearman) and the top-N positively correlated genes define the signature.
#'
#' @param obj Seurat object.
#' @param assay Character scalar assay name (default "RNA").
#' @param slot Character scalar expression slot (default "data").
#' @param weight_cols Optional character vector of metadata columns containing per-cell
#'   archetype weights (e.g., c("A5_A1","A5_A2",...)). If NULL, tries auto-detection.
#' @param weights_mat Optional numeric matrix/data.frame of per-cell weights (cells × k).
#'   Rownames must be Seurat cell IDs. If provided, this is used instead of `weight_cols`.
#' @param method Correlation method: "spearman" (default) or "pearson".
#' @param top_n Integer >= 1. Number of top positively correlated genes per archetype
#'   to return as the signature.
#' @param min_detect_frac Numeric in (0,1]. Filter genes detected (>0) in at least this
#'   fraction of cells (default 0.05).
#' @param max_genes Optional integer. If provided, randomly subsample genes to this number
#'   after filtering (useful for quick dev runs).
#' @param seed Integer seed used only when max_genes is not NULL.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with:
#'   - ranked: data.frame with columns gene, archetype, cor, n_cells_used
#'   - signatures: named list of character vectors (top_n genes per archetype)
#'   - spec: list describing inputs used (assay/slot/method/top_n/etc.)
#' 
#' #' @examples
#' \dontrun{
#' weight_cols <- paste0("A5_A", 1:5)
#'
#' sig <- compute_signature_cor(
#'   obj = obj,
#'   assay = "RNA",
#'   slot = "data",
#'   weight_cols = weight_cols,
#'   weights_mat = NULL,
#'   method = "spearman",
#'   top_n = 100L,
#'   min_detect_frac = 0.05,
#'   max_genes = NULL,
#'   seed = 1L,
#'   verbose = TRUE
#' )
#'
#' # Top signature genes for archetype 1
#' head(sig$signatures[[weight_cols[1]]], 20)
#'
#' # Ranked genes for archetype 1 (highest correlations first)
#' ranked_a1 <- sig$ranked[sig$ranked$archetype == weight_cols[1], ]
#' ranked_a1 <- ranked_a1[order(ranked_a1$cor, decreasing = TRUE), ]
#' head(ranked_a1, 20)
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
    if (is.data.frame(weights_mat)) weights_mat <- as.matrix(weights_mat)
    if (!is.matrix(weights_mat) || !is.numeric(weights_mat)) {
      stop("`weights_mat` must be a numeric matrix/data.frame (cells × k).", call. = FALSE)
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
      pat <- "^A\\d+_A\\d+$" # e.g. A5_A1
      weight_cols <- colnames(md)[grepl(pat, colnames(md))]
      if (length(weight_cols) == 0L) {
        stop("Could not auto-detect `weight_cols` in obj@meta.data. Pass `weight_cols` explicitly.", call. = FALSE)
      }
    }
    if (!is.character(weight_cols) || length(weight_cols) < 1L || any(is.na(weight_cols)) || any(!nzchar(weight_cols))) {
      stop("`weight_cols` must be a non-empty character vector or NULL (auto-detect).", call. = FALSE)
    }
    missing_w <- setdiff(weight_cols, colnames(md))
    if (length(missing_w) > 0L) {
      stop("Missing `weight_cols` in obj@meta.data: ", paste(missing_w, collapse = ", "), call. = FALSE)
    }

    # coerce numeric
    W <- md[, weight_cols, drop = FALSE]
    for (w in weight_cols) {
      if (!is.numeric(W[[w]])) {
        suppressWarnings(xn <- as.numeric(W[[w]]))
        if (all(is.na(xn))) stop("Weight column '", w, "' cannot be coerced to numeric.", call. = FALSE)
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

  # Align cells
  keep_cells2 <- intersect(use_cells, colnames(X))
  if (length(keep_cells2) == 0L) {
    stop("No overlapping cell IDs between weights and expression matrix.", call. = FALSE)
  }

  # subset and align
  W <- W[keep_cells2, , drop = FALSE]
  X <- X[, keep_cells2, drop = FALSE]

  # -----------------------------
  # [FILTER] genes by detection fraction
  # -----------------------------
  # Detection defined as >0 in the chosen slot
  det_frac <- Matrix::rowMeans(X > 0)
  keep_genes <- names(det_frac)[det_frac >= min_detect_frac]
  if (length(keep_genes) == 0L) {
    stop("No genes pass min_detect_frac = ", min_detect_frac, ".", call. = FALSE)
  }
  X <- X[keep_genes, , drop = FALSE]

  if (!is.null(max_genes) && length(keep_genes) > max_genes) {
    set.seed(seed)
    keep_sub <- sample(rownames(X), size = max_genes, replace = FALSE)
    X <- X[keep_sub, , drop = FALSE]
  }

  # -----------------------------
  # [COMPUTE] correlation per gene per archetype
  # - Efficient approach:
  #   compute cor(gene_vector, archetype_weight_vector) for each archetype separately.
  #   Use apply over genes; OK for moderate gene counts.
  # -----------------------------
  ranked_list <- vector("list", length(arch_names))
  names(ranked_list) <- arch_names

  # Ensure X is a base matrix for apply speed (can be memory heavy; keep filtered)
  X_dense <- as.matrix(X)

  for (a in arch_names) {
    w <- W[, a]
    w <- as.numeric(w)

    if (verbose) message("Computing correlations for ", a, " ...")

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
      max_genes = max_genes,
      seed = seed,
      archetypes = arch_names,
      n_cells = ncol(X_dense),
      n_genes = nrow(X_dense)
    )
  )
}