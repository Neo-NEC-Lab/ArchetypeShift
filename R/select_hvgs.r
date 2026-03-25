#' @title Select Highly Variable Genes for PCA
#'
#' @description
#' Selects highly variable genes (HVGs) from a Seurat object for dimensionality reduction.
#' Optionally computes HVGs per group and combines them using intersection or union.
#'
#' @param obj Seurat object.
#' @param assay Assay to use (e.g., \code{"RNA"}).
#' @param method How to combine per-group HVGs when \code{per_group=TRUE}. One of
#'   \code{"hvg_intersect"} or \code{"hvg_union"}.
#' @param n_features Integer; target number of HVGs to return.
#' @param per_group Logical; if TRUE, compute HVGs separately per group and combine.
#' @param group_col Metadata column to define groups when \code{per_group=TRUE}.
#' @param exclude_mt Logical; if TRUE, drops genes matching \code{^MT-} or \code{^mt-}.
#' @param exclude_regex Optional character vector of additional regex patterns to exclude
#'   (applied after HVG selection). For example, \code{c("^(RPL|RPS)", "^HBA", "^HBB")}.
#' @param verbose Logical; print brief status messages.
#'
#' @details
#' - Union retains group-specific signals (often preferred when biology differs by group).
#' - Intersection emphasizes shared variability across groups (more conservative).
#' - If \code{exclude_mt=FALSE}, mitochondrial genes are allowed.
#' - \code{exclude_regex} provides additional customizable filtering.
#'
#' @return Character vector of HVG gene names (length \code{<= n_features}).
#'
#' @examples
#' \dontrun{
#' # Keep MT genes (not typical, but allowed)
#' hvgs_keep_mt <- select_hvgs(
#'   obj = sobj,
#'   assay = "RNA",
#'   per_group = TRUE,
#'   group_col = "condition2",
#'   method = "hvg_union",
#'   n_features = 3000,
#'   exclude_mt = FALSE
#' )
#'
#' # Exclude MT + ribosomal
#' hvgs_no_ribo <- select_hvgs(
#'   obj = sobj,
#'   assay = "RNA",
#'   per_group = FALSE,
#'   group_col = "condition2",
#'   method = "hvg_union",
#'   n_features = 3000,
#'   exclude_mt = TRUE,
#'   exclude_regex = c("^(RPL|RPS)")
#' )
#' }
#'
#' @export
select_hvgs <- function(
  obj,
  assay = "RNA",
  method = c("hvg_intersect", "hvg_union"),
  n_features = 3000,
  per_group = TRUE,
  group_col = "condition",
  exclude_mt = TRUE,
  exclude_regex = NULL,
  verbose = TRUE
) {
  method <- match.arg(method)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!assay %in% names(obj@assays)) stop("Assay '", assay, "' not found in object.")

  old_assay <- DefaultAssay(obj)
  on.exit(DefaultAssay(obj) <- old_assay, add = TRUE)
  DefaultAssay(obj) <- assay

  meta <- obj@meta.data

  if (per_group) {
    if (!group_col %in% colnames(meta)) stop("group_col '", group_col, "' not found in meta.data")
    groups <- as.character(meta[[group_col]])
    if (all(is.na(groups))) stop("group_col '", group_col, "' is all NA")

    group_cells <- split(colnames(obj), groups)
    if (verbose) message("Computing HVGs per group (", length(group_cells), " groups) using method=", method)

    hvg_sets <- lapply(names(group_cells), function(g) {
      cells <- group_cells[[g]]
      x <- subset(obj, cells = cells)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
      Seurat::VariableFeatures(x)
    })
    genes <- if (method == "hvg_intersect") Reduce(intersect, hvg_sets) else Reduce(union, hvg_sets)
    if (length(genes) > n_features) genes <- genes[seq_len(n_features)]
  } else {
    if (verbose) message("Computing HVGs on all cells (n_features=", n_features, ")")
    x <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
    genes <- Seurat::VariableFeatures(x)
  }

  if (exclude_mt) {
    genes <- genes[!grepl("^MT-|^mt-", genes)]
  }

  if (!is.null(exclude_regex)) {
    if (!is.character(exclude_regex)) stop("exclude_regex must be a character vector of regex patterns or NULL.")
    for (rgx in exclude_regex) {
      genes <- genes[!grepl(rgx, genes)]
    }
  }

  unique(genes)
}