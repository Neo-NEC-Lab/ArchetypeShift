#' @title Define a Consistent Gene Universe for Archetype Analysis
#'
#' @description
#' Constructs a reproducible background gene universe from a Seurat object for use in
#' archetype fitting and enrichment analyses. Applies configurable gene filters and a
#' minimum detection threshold to remove rarely expressed genes.
#'
#' @param obj Seurat object.
#' @param assay Character assay name (e.g., \code{"RNA"}).
#' @param slot Seurat assay slot to use (typically \code{"counts"} or \code{"data"}).
#' @param exclude_mt Logical; if TRUE, drops genes matching \code{^MT-} or \code{^mt-}.
#' @param exclude_ribo Logical; if TRUE, drops genes matching \code{^(RPL|RPS)}.
#' @param exclude_malat1 Logical; if TRUE, drops the gene \code{"MALAT1"} (exact match).
#' @param exclude_regex Optional character vector of additional regex patterns to exclude.
#' @param min_cells_expressed Integer; keep genes detected (nonzero) in at least this many cells.
#' @param verbose Logical; print summary counts.
#' @param write_path Optional file path. If provided, writes the resulting gene universe
#'   (one gene per line) to this file.
#'
#' @details
#' The gene universe is used as the background for enrichment calculations (e.g., Fisher's exact test),
#' and should match the feature space used when computing archetype gene scores/loadings.
#' To avoid mismatch issues, freeze this universe once per run and reuse it downstream.
#'
#' @return Character vector of gene symbols.
#'
#' @examples
#' \dontrun{
#' gene_universe <- make_gene_universe(
#'   obj = sobj,
#'   assay = "RNA",
#'   slot = "counts",
#'   exclude_mt = TRUE,
#'   exclude_ribo = TRUE,
#'   exclude_malat1 = TRUE,
#'   min_cells_expressed = 10
#' )
#' }
#'
#' @export
make_gene_universe <- function(
  obj,
  assay = "RNA",
  slot = "counts",
  exclude_mt = TRUE,
  exclude_ribo = TRUE,
  exclude_malat1 = TRUE,
  exclude_regex = NULL,
  min_cells_expressed = 10L,
  verbose = TRUE,
  write_path = NULL
) {
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!assay %in% names(obj@assays)) stop("Assay '", assay, "' not found in object.")
  if (!slot %in% c("counts", "data", "scale.data")) stop("slot must be one of: counts, data, scale.data")

  old_assay <- Seurat::DefaultAssay(obj)
  on.exit(Seurat::DefaultAssay(obj) <- old_assay, add = TRUE)
  Seurat::DefaultAssay(obj) <- assay

  mat <- Seurat::GetAssayData(obj, assay = assay, slot = slot)
  if (is.null(rownames(mat))) stop("Expression matrix has no rownames (gene identifiers).")

  genes_all <- rownames(mat)
  keep <- rep(TRUE, length(genes_all))

  if (exclude_mt) {
    keep <- keep & !grepl("^MT-|^mt-", genes_all)
  }
  if (exclude_ribo) {
    keep <- keep & !grepl("^(RPL|RPS)", genes_all)
  }
  if (exclude_malat1) {
    keep <- keep & genes_all != "MALAT1"
  }
  if (!is.null(exclude_regex)) {
    if (!is.character(exclude_regex)) stop("exclude_regex must be a character vector of regex patterns or NULL.")
    for (rgx in exclude_regex) {
      keep <- keep & !grepl(rgx, genes_all)
    }
  }

  gene_universe <- genes_all[keep]

  if (!is.integer(min_cells_expressed)) min_cells_expressed <- as.integer(min_cells_expressed)
  if (min_cells_expressed < 1L) stop("min_cells_expressed must be >= 1")

  detected_cells <- Matrix::rowSums(mat[gene_universe, , drop = FALSE] > 0)
  gene_universe <- gene_universe[detected_cells >= min_cells_expressed]

  if (verbose) {
    message("Genes total: ", length(genes_all))
    message("Genes after filters: ", sum(keep))
    message("Genes after min_cells_expressed (>= ", min_cells_expressed, "): ", length(gene_universe))
  }

  if (!is.null(write_path)) {
    write_path <- normalizePath(write_path, winslash = "/", mustWork = FALSE)
    dir.create(dirname(write_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(gene_universe, con = write_path)
    if (verbose) message("Wrote gene universe: ", write_path)
  }

  gene_universe
}
