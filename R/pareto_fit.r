#' Fit ParetoTI archetype models from a Seurat object
#'
#' Runs a reproducible ParetoTI/PCHA workflow on a Seurat object by:
#' defining a frozen gene universe, selecting highly variable genes (HVGs),
#' computing PCA on normalized expression, running a ParetoTI k-scan, and
#' optionally fitting selected target-k archetype models.
#'
#' This function is designed as the core fitting step in the ArchetypeShift
#' workflow. It writes the main model artifacts to a run-specific output
#' directory and returns a compact summary bundle describing the analysis.
#'
#' @param obj A Seurat object.
#' @param cfg Optional configuration list returned by [check_input()]. If `NULL`,
#'   configuration is resolved internally using `assay`, `group_col`,
#'   `sample_col`, and `cluster_col`.
#' @param assay Character scalar giving the assay name. Default is `"RNA"`.
#' @param group_col Optional metadata column giving the condition/group label,
#'   used only when `cfg = NULL`.
#' @param sample_col Optional metadata column giving the sample identifier, used
#'   only when `cfg = NULL`.
#' @param cluster_col Optional metadata column giving the cluster annotation,
#'   used only when `cfg = NULL`.
#' @param out_dir Character scalar giving the output directory root. A run
#'   folder is created within this directory.
#' @param run_id Optional character scalar run identifier. If `NULL`, a
#'   timestamp-based identifier is used.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is
#'   `TRUE`.
#' @param seed Integer random seed for reproducibility. Default is `1`.
#' @param gene_slot Character scalar giving the Seurat slot used to define the
#'   gene universe, typically `"counts"`.
#' @param min_cells_expressed Integer >= 1 giving the minimum number of cells in
#'   which a gene must be expressed to enter the gene universe.
#' @param exclude_mt Logical; if `TRUE`, excludes mitochondrial genes from the
#'   gene universe. Default is `TRUE`.
#' @param exclude_ribo Logical; if `TRUE`, excludes ribosomal genes from the
#'   gene universe. Default is `TRUE`.
#' @param exclude_malat1 Logical; if `TRUE`, excludes `MALAT1` from the gene
#'   universe. Default is `TRUE`.
#' @param exclude_gene_regex Optional character vector of additional regex
#'   patterns used to exclude genes from the gene universe.
#' @param hvgs_n Integer >= 1 giving the number of HVGs used for PCA. Default is
#'   `3000L`.
#' @param hvgs_per_group Logical; if `TRUE`, computes HVGs per group and then
#'   combines them according to `hvgs_method`. Default is `TRUE`.
#' @param hvgs_method Method used to combine per-group HVGs, either
#'   `"hvg_union"` or `"hvg_intersect"`.
#' @param hvgs_exclude_mt Logical; if `TRUE`, excludes mitochondrial genes from
#'   HVG selection. Default is `TRUE`.
#' @param hvgs_exclude_regex Optional character vector of regex patterns used to
#'   exclude genes during HVG selection.
#' @param pcs_use Integer >= 2 giving the number of principal components used
#'   for ParetoTI fitting. Default is `30L`.
#' @param ks Integer vector of k values used for the ParetoTI k-scan. Default is
#'   `3:10`.
#' @param k_targets Optional integer vector of target k values to fit with
#'   [ParetoTI::fit_pch()]. Default is `c(4L, 5L)`.
#' @param do_bootstrap Logical; if `TRUE`, enables bootstrap fitting in
#'   [ParetoTI::fit_pch()] for each `k_targets` value. Default is `FALSE`.
#' @param bootstrap_n Integer >= 1 giving the number of bootstrap replicates
#'   used when `do_bootstrap = TRUE`. Default is `20L`.
#' @param bootstrap_prop Numeric in `(0, 1]` giving the fraction of cells used
#'   per bootstrap replicate when `do_bootstrap = TRUE`. Default is `0.70`.
#' @param bootstrap_type Bootstrap mode passed to [ParetoTI::fit_pch()], either
#'   `"m"` or `"s"`.
#'
#' @return A named list with analysis metadata and key output paths:
#' \describe{
#'   \item{run_dir}{Run directory path.}
#'   \item{cfg}{Resolved configuration list.}
#'   \item{ks}{Integer vector used for k-scan.}
#'   \item{k_targets}{Integer vector of fitted target k values.}
#'   \item{pcs_use}{Number of PCs used for ParetoTI fitting.}
#'   \item{hvgs}{Character vector of HVGs used for PCA.}
#'   \item{gene_universe}{Character vector of genes in the frozen gene universe.}
#'   \item{pca_rds}{Path to the saved PCA object.}
#'   \item{kscan_rds}{Path to the saved ParetoTI k-scan object.}
#'   \item{archetype_fit_paths}{Character vector of saved target-k fit paths.}
#'   \item{seed}{Random seed used.}
#'   \item{do_bootstrap}{Logical indicating whether bootstrap fitting was used.}
#'   \item{bootstrap_n}{Bootstrap replicate count.}
#'   \item{bootstrap_prop}{Bootstrap sampling proportion.}
#'   \item{bootstrap_type}{Bootstrap mode used.}
#' }
#'
#' @details
#' The gene universe is constructed from `gene_slot` and written to disk so that
#' downstream enrichment analyses can use a stable background. HVGs are selected
#' independently from normalized expression, intersected with the frozen gene
#' universe, and then used for PCA via [irlba::prcomp_irlba()]. The PCA matrix is
#' transposed to PCs x cells before ParetoTI fitting.
#'
#' The function checks for required R package dependencies and for Python modules
#' required by ParetoTI through `reticulate`, including `py_pcha` and
#' `geosketch`.
#'
#' @examples
#' \dontrun{
#' cfg <- check_input(
#'   obj = combined_no_IA,
#'   assay = "RNA",
#'   group_col = "condition2",
#'   sample_col = "orig.ident",
#'   cluster_col = "merged_cluster_annotations"
#' )
#'
#' fit_res <- pareto_fit(
#'   obj = combined_no_IA,
#'   cfg = cfg,
#'   out_dir = "outputs",
#'   ks = 3:10,
#'   k_targets = c(4L, 5L),
#'   do_bootstrap = FALSE
#' )
#'
#' fit_res$run_dir
#' fit_res$archetype_fit_paths
#' }
#'
#' @export
pareto_fit <- function(
  obj,
  cfg = NULL,
  assay = "RNA",
  group_col = NULL,
  sample_col = NULL,
  cluster_col = NULL,
  out_dir = "outputs",
  run_id = NULL,
  verbose = TRUE,
  seed = 1,
  gene_slot = "counts",
  min_cells_expressed = 10L,
  exclude_mt = TRUE,
  exclude_ribo = TRUE,
  exclude_malat1 = TRUE,
  exclude_gene_regex = NULL,
  hvgs_n = 3000L,
  hvgs_per_group = TRUE,
  hvgs_method = c("hvg_union", "hvg_intersect"),
  hvgs_exclude_mt = TRUE,
  hvgs_exclude_regex = NULL,
  pcs_use = 30L,
  ks = 3:10,
  k_targets = c(4L, 5L),
  do_bootstrap = FALSE,
  bootstrap_n = 20L,
  bootstrap_prop = 0.70,
  bootstrap_type = c("m", "s")
) {
  # ----------------------------
  # [CHECKS] core inputs
  # ----------------------------
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }
  if (!is.null(cfg) && !is.list(cfg)) {
    stop("`cfg` must be NULL or a list returned by `check_input()`.", call. = FALSE)
  }
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
    stop("`assay` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)) {
    stop("`out_dir` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.null(run_id) && (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) || !nzchar(run_id))) {
    stop("`run_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single numeric value.", call. = FALSE)
  }
  seed <- as.integer(seed)
  if (!is.character(gene_slot) || length(gene_slot) != 1L || is.na(gene_slot) || !nzchar(gene_slot)) {
    stop("`gene_slot` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.numeric(min_cells_expressed) || length(min_cells_expressed) != 1L ||
      is.na(min_cells_expressed) || min_cells_expressed < 1) {
    stop("`min_cells_expressed` must be a single integer >= 1.", call. = FALSE)
  }
  min_cells_expressed <- as.integer(min_cells_expressed)

  logical_args <- list(
    exclude_mt = exclude_mt,
    exclude_ribo = exclude_ribo,
    exclude_malat1 = exclude_malat1,
    hvgs_per_group = hvgs_per_group,
    hvgs_exclude_mt = hvgs_exclude_mt,
    do_bootstrap = do_bootstrap
  )
  for (nm in names(logical_args)) {
    if (!is.logical(logical_args[[nm]]) || length(logical_args[[nm]]) != 1L) {
      stop("`", nm, "` must be TRUE/FALSE.", call. = FALSE)
    }
  }

  check_regex_arg <- function(x, arg_name) {
    if (!is.null(x)) {
      if (!is.character(x) || length(x) < 1L || any(is.na(x)) || any(!nzchar(x))) {
        stop("`", arg_name, "` must be NULL or a non-empty character vector.", call. = FALSE)
      }
    }
  }
  check_regex_arg(exclude_gene_regex, "exclude_gene_regex")
  check_regex_arg(hvgs_exclude_regex, "hvgs_exclude_regex")

  if (!is.numeric(hvgs_n) || length(hvgs_n) != 1L || is.na(hvgs_n) || hvgs_n < 1) {
    stop("`hvgs_n` must be a single integer >= 1.", call. = FALSE)
  }
  hvgs_n <- as.integer(hvgs_n)

  hvgs_method <- match.arg(hvgs_method)
  bootstrap_type <- match.arg(bootstrap_type)

  if (!is.numeric(pcs_use) || length(pcs_use) != 1L || is.na(pcs_use) || pcs_use < 2) {
    stop("`pcs_use` must be a single integer >= 2.", call. = FALSE)
  }
  pcs_use <- as.integer(pcs_use)

  if (!is.numeric(ks) || length(ks) < 1L || any(is.na(ks))) {
    stop("`ks` must be a non-empty numeric vector.", call. = FALSE)
  }
  ks <- sort(unique(as.integer(ks)))
  if (any(ks < 2L)) {
    stop("All values in `ks` must be >= 2.", call. = FALSE)
  }

  if (!is.null(k_targets)) {
    if (!is.numeric(k_targets) || any(is.na(k_targets))) {
      stop("`k_targets` must be NULL or a numeric vector.", call. = FALSE)
    }
    k_targets <- sort(unique(as.integer(k_targets)))
    if (any(k_targets < 2L)) {
      stop("All values in `k_targets` must be >= 2.", call. = FALSE)
    }
  }

  if (!is.numeric(bootstrap_n) || length(bootstrap_n) != 1L || is.na(bootstrap_n) || bootstrap_n < 1) {
    stop("`bootstrap_n` must be a single integer >= 1.", call. = FALSE)
  }
  bootstrap_n <- as.integer(bootstrap_n)

  if (!is.numeric(bootstrap_prop) || length(bootstrap_prop) != 1L ||
      is.na(bootstrap_prop) || bootstrap_prop <= 0 || bootstrap_prop > 1) {
    stop("`bootstrap_prop` must be a single numeric value in (0, 1].", call. = FALSE)
  }
  bootstrap_prop <- as.numeric(bootstrap_prop)

  set.seed(seed)

  if (is.null(run_id)) {
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }

  # ----------------------------
  # [CONFIG]
  # ----------------------------
  if (is.null(cfg)) {
    cfg <- check_input(
      obj = obj,
      assay = assay,
      group_col = group_col,
      sample_col = sample_col,
      cluster_col = cluster_col,
      verbose = verbose
    )
  } else {
    if (is.null(cfg$assay) || !is.character(cfg$assay) || length(cfg$assay) != 1L || !nzchar(cfg$assay)) {
      stop("When `cfg` is supplied, it must contain a valid `assay` entry.", call. = FALSE)
    }
  }

  # ----------------------------
  # [DEPENDENCIES]
  # ----------------------------
  if (!requireNamespace("ParetoTI", quietly = TRUE)) {
    stop("Package 'ParetoTI' is required.", call. = FALSE)
  }
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for PCA.", call. = FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required because ParetoTI depends on Python modules.", call. = FALSE)
  }
  if (!reticulate::py_module_available("py_pcha")) {
    stop("Python module 'py_pcha' is not available to reticulate.", call. = FALSE)
  }
  if (!reticulate::py_module_available("geosketch")) {
    stop("Python module 'geosketch' is not available to reticulate.", call. = FALSE)
  }

  # ----------------------------
  # [OUTPUT]
  # ----------------------------
  out_root <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  run_dir <- file.path(out_root, paste0("run_", run_id))

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "bundle"), recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message("run_dir: ", run_dir)
  }

  # ----------------------------
  # [GENE UNIVERSE]
  # ----------------------------
  gene_universe <- make_gene_universe(
    obj = obj,
    assay = cfg$assay,
    slot = gene_slot,
    exclude_mt = exclude_mt,
    exclude_ribo = exclude_ribo,
    exclude_malat1 = exclude_malat1,
    exclude_regex = exclude_gene_regex,
    min_cells_expressed = min_cells_expressed,
    verbose = verbose,
    write_path = file.path(run_dir, "bundle", "gene_universe.txt")
  )

  # ----------------------------
  # [HVG]
  # ----------------------------
  hvgs <- select_hvgs(
    obj = obj,
    assay = cfg$assay,
    method = hvgs_method,
    n_features = hvgs_n,
    per_group = hvgs_per_group,
    group_col = cfg$group_col,
    exclude_mt = hvgs_exclude_mt,
    exclude_regex = hvgs_exclude_regex,
    verbose = verbose
  )

  hvgs <- intersect(hvgs, gene_universe)
  if (length(hvgs) < 100L) {
    stop(
      "Too few HVGs remain after intersecting with the gene universe (n=",
      length(hvgs),
      ").",
      call. = FALSE
    )
  }

  writeLines(hvgs, con = file.path(run_dir, "bundle", paste0("hvgs_", length(hvgs), ".txt")))

  # ----------------------------
  # [EXTRACT] normalized expression for PCA
  # ----------------------------
  old_assay <- Seurat::DefaultAssay(obj)
  on.exit(Seurat::DefaultAssay(obj) <- old_assay, add = TRUE)
  Seurat::DefaultAssay(obj) <- cfg$assay

  expr_mat <- Seurat::GetAssayData(obj, assay = cfg$assay, slot = "data")
  if (nrow(expr_mat) == 0L || ncol(expr_mat) == 0L) {
    stop(
      "Assay '", cfg$assay, "' slot 'data' is empty. Run NormalizeData() first.",
      call. = FALSE
    )
  }

  hvgs_in_data <- intersect(hvgs, rownames(expr_mat))
  if (length(hvgs_in_data) < 100L) {
    stop(
      "Too few HVGs are present in slot='data' after intersection (n=",
      length(hvgs_in_data),
      ").",
      call. = FALSE
    )
  }

  if (verbose && length(hvgs_in_data) < length(hvgs)) {
    message(
      "Dropping ",
      length(hvgs) - length(hvgs_in_data),
      " HVGs not found in slot='data'. Using n=",
      length(hvgs_in_data),
      "."
    )
  }
  hvgs <- hvgs_in_data

  expr_hvg <- expr_mat[hvgs, , drop = FALSE]

  # ----------------------------
  # [PCA]
  # ----------------------------
  pca_res <- irlba::prcomp_irlba(
    t(as.matrix(expr_hvg)),
    n = pcs_use,
    center = TRUE,
    scale. = FALSE
  )

  pcs_mat <- as.matrix(pca_res$x)
  rownames(pcs_mat) <- colnames(expr_hvg)
  colnames(pcs_mat) <- paste0("PC", seq_len(ncol(pcs_mat)))

  pca_for_pareto <- t(pcs_mat)
  pca_rds <- file.path(run_dir, "models", "pca_irlba.rds")
  saveRDS(pca_res, pca_rds)

  # ----------------------------
  # [K-SCAN]
  # ----------------------------
  if (verbose) {
    message("Running k_fit_pch for ks = ", paste(ks, collapse = ","))
  }

  k_scan <- ParetoTI::k_fit_pch(
    data = pca_for_pareto,
    ks = ks,
    bootstrap = FALSE,
    simplex = FALSE,
    var_in_dims = FALSE,
    normalise_var = TRUE
  )

  kscan_rds <- file.path(
    run_dir,
    "models",
    paste0("pch_k_scan_", min(ks), "_", max(ks), ".rds")
  )
  saveRDS(k_scan, kscan_rds)

  # ----------------------------
  # [TARGET K FITS]
  # ----------------------------
  archetype_fit_paths <- character(0)

  if (!is.null(k_targets) && length(k_targets) > 0L) {
    for (k in k_targets) {
      if (verbose) {
        message("Fitting archetypes for k = ", k)
      }

      fit_k <- ParetoTI::fit_pch(
        data = pca_for_pareto,
        k = k,
        bootstrap = isTRUE(do_bootstrap),
        bootstrap_N = if (isTRUE(do_bootstrap)) bootstrap_n else 0L,
        sample_prop = if (isTRUE(do_bootstrap)) bootstrap_prop else 1,
        bootstrap_type = bootstrap_type,
        seed = seed,
        simplex = FALSE,
        var_in_dims = FALSE,
        normalise_var = TRUE
      )

      fit_path <- file.path(run_dir, "models", paste0("pch_fit_k", k, ".rds"))
      saveRDS(fit_k, fit_path)
      archetype_fit_paths <- c(archetype_fit_paths, fit_path)
    }
  }

  # ----------------------------
  # [RETURN]
  # ----------------------------
  res <- list(
    run_dir = run_dir,
    cfg = cfg,
    ks = ks,
    k_targets = k_targets,
    pcs_use = pcs_use,
    hvgs = hvgs,
    gene_universe = gene_universe,
    pca_rds = pca_rds,
    kscan_rds = kscan_rds,
    archetype_fit_paths = archetype_fit_paths,
    seed = seed,
    do_bootstrap = isTRUE(do_bootstrap),
    bootstrap_n = bootstrap_n,
    bootstrap_prop = bootstrap_prop,
    bootstrap_type = bootstrap_type
  )

  saveRDS(res, file.path(run_dir, "models", "run_bundle.rds"))
  res
}
