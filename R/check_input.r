#' @title Check Input Seurat Object and Required Metadata
#'
#' @description
#' Validates that a Seurat object contains a requested assay and the metadata columns
#' required for downstream archetype fitting and annotation workflows. The function can
#' auto-detect commonly used group, sample, and cluster columns from \code{obj@meta.data}
#' via \code{detect_meta_cols()}.
#'
#' @param obj Seurat object.
#' @param assay Character assay name (e.g., \code{"RNA"} or \code{"SCT"}).
#' @param group_col Optional metadata column name for condition/group (e.g., \code{"condition2"}).
#'   If \code{NULL}, attempts to auto-detect using \code{detect_meta_cols()}.
#' @param sample_col Optional sample column name (e.g., \code{"orig.ident"}). If \code{NULL},
#'   attempts to auto-detect using \code{detect_meta_cols()}.
#' @param cluster_col Optional cluster/cell-type column name. If \code{NULL}, attempts to
#'   auto-detect using \code{detect_meta_cols()}.
#' @param min_cells_per_sample Minimum cells required per sample to issue a warning
#'   (does not stop execution).
#' @param verbose Logical; print resolved columns and basic information.
#'
#' @return A named list with resolved configuration values:
#' \describe{
#'   \item{assay}{Assay name used.}
#'   \item{group_col}{Group/condition column name.}
#'   \item{sample_col}{Sample column name.}
#'   \item{cluster_col}{Cluster column name.}
#'   \item{min_cells_per_sample}{Minimum cells per sample used for warnings.}
#' }
#'
#' @examples
#' \dontrun{
#' cfg <- check_input(
#'   obj = sobj,
#'   assay = "RNA",
#'   group_col = "condition2",
#'   sample_col = "orig.ident",
#'   cluster_col = "merged_cluster_annotations"
#' )
#'
#' # Or rely on auto-detection:
#' cfg_auto <- check_input(
#'   obj = sobj,
#'   assay = "RNA",
#'   group_col = NULL,
#'   sample_col = NULL,
#'   cluster_col = NULL
#' )
#' }
#'
#' @export
check_input <- function(
  obj,
  assay = "RNA",
  group_col = NULL,
  sample_col = NULL,
  cluster_col = NULL,
  min_cells_per_sample = 100,
  verbose = TRUE
) {
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!assay %in% names(obj@assays)) {
    stop(
      "Assay '", assay, "' not found in object. Available: ",
      paste(names(obj@assays), collapse = ", ")
    )
  }

  meta <- obj@meta.data

  # Resolve group/sample/cluster via a single wrapper call
  hits <- detect_meta_cols(meta)
  if (is.null(group_col)) group_col <- hits$group_col
  if (is.null(sample_col)) sample_col <- hits$sample_col
  if (is.null(cluster_col)) cluster_col <- hits$cluster_col

  # Resolve group_col
  if (is.null(group_col) || length(group_col) != 1 || !nzchar(group_col) || is.na(group_col)) {
    stop("Could not auto-detect a group/condition column; please supply `group_col`.")
  }
  if (!group_col %in% colnames(meta)) stop("group_col '", group_col, "' not found in meta.data")
  if (all(is.na(meta[[group_col]]))) stop("group_col '", group_col, "' is all NA.")

  # Resolve sample_col
  if (is.null(sample_col) || length(sample_col) != 1 || !nzchar(sample_col) || is.na(sample_col)) {
    stop("Could not auto-detect a sample column; please supply `sample_col`.")
  }
  if (!sample_col %in% colnames(meta)) stop("sample_col '", sample_col, "' not found in meta.data")
  if (all(is.na(meta[[sample_col]]))) stop("sample_col '", sample_col, "' is all NA.")

  # Resolve cluster_col
  if (is.null(cluster_col) || length(cluster_col) != 1 || !nzchar(cluster_col) || is.na(cluster_col)) {
    stop("Could not auto-detect a cluster column; please supply `cluster_col`.")
  }
  if (!cluster_col %in% colnames(meta)) stop("cluster_col '", cluster_col, "' not found in meta.data")
  if (all(is.na(meta[[cluster_col]]))) stop("cluster_col '", cluster_col, "' is all NA.")

  if (verbose) {
    message(
      "Using assay=", assay,
      "; group_col=", group_col,
      "; sample_col=", sample_col,
      "; cluster_col=", cluster_col
    )
  }

  tbl <- table(meta[[sample_col]])
  if (any(tbl < min_cells_per_sample)) {
    small <- paste(names(tbl)[tbl < min_cells_per_sample], collapse = ", ")
    warning("Some samples have < ", min_cells_per_sample, " cells: ", small)
  }

  list(
    assay = assay,
    group_col = group_col,
    sample_col = sample_col,
    cluster_col = cluster_col,
    min_cells_per_sample = min_cells_per_sample
  )
}