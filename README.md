# ArchetypeShift

ArchetypeShift is an R package for fitting, annotating, and comparing transcriptional archetype programs in single-cell RNA-seq data, with particular emphasis on biological interpretation through correlation-based gene signatures, IPA canonical pathways, upstream regulators, and cluster/state localization. It is designed for Seurat-based workflows and supports reproducible preprocessing, ParetoTI/PCHA model fitting, archetype annotation, and downstream condition- and cluster-level comparisons. The package is especially useful when the goal is not only to identify archetypes, but also to interpret them as biologically meaningful programs across cell states and experimental conditions.

## Installation

You can install ArchetypeShift with:

```r
devtools::install_github("Neo-NEC-Lab/ArchetypeShift")
library(ArchetypeShift)
```

## What This Package Does

This workflow is designed to help you:

- validate Seurat inputs and metadata columns
- define a reproducible gene universe
- select highly variable genes for PCA
- fit ParetoTI/PCHA archetype models across a scan of `k`
- extract per-cell archetype weights
- compute correlation-based archetype gene signatures
- summarize archetype occupancy at the sample and sample-by-cluster levels
- test replicate-aware condition contrasts
- generate publication-ready QC and summary plots

The package was developed for epithelial archetype analysis in scRNA-seq, but the core functions are general enough to be adapted to other cell systems.

## Core Workflow

At a high level, the repository supports the following analysis path:

1. Validate the Seurat object and resolve metadata columns with `check_input()`
2. Define a stable gene universe with `make_gene_universe()`
3. Select HVGs with `select_hvgs()`
4. Fit ParetoTI archetype models with `pareto_fit()`
5. Extract weights with `pcha_extract_weights()`
6. Add dominant archetype assignments with `phca_dominant_archetype()`
7. Compute gene signatures with `compute_signature_cor()`
8. Summarize archetype occupancy with `compute_occupancy()`
9. Test replicate-aware contrasts with `compute_contrast()`
10. Generate downstream occupancy, delta, QC, and enrichment plots

## Repository Contents

Some of the main scripts in this repository are:

- `check_input.r`: validates Seurat objects and required metadata
- `detect_group_col.r`, `detect_meta_cols.r`, `detect_sample_col.r`, `detect_cluster_col.R`: helper functions for metadata auto-detection
- `make_gene_universe.r`: defines a frozen background gene universe
- `select_hvgs.r`: selects HVGs for PCA
- `pareto_fit2.r`: main ParetoTI/PCHA fitting workflow with k-scan and bootstrap support
- `phca_extract_weights.r`: extracts a standardized cells x k archetype weight matrix
- `phca_dominant_archetype.r`: assigns dominant archetypes per cell
- `phca_scan_metrics.r`: summarizes k-scan metrics
- `compute_signature_cor.r`: computes correlation-based archetype signatures
- `compute_occupancy.r`: replicate-aware sample and sample-by-cluster occupancy summaries
- `compute_contrasts.r`: replicate-aware condition contrasts for archetype occupancy
- `compute_phca_delta_heatmap_df.r`, `compute_phca_occupancy.r`, `compute_phca_occupancy_sample_cluster.r`: helper summaries for downstream plotting
- `plot_phca_delta_heatmap.r`, `plot_phca_occupancy_heatmap.r`, `plot_phca_occupancy_stacked.r`, `plot_phca_qc_hist.r`: plotting utilities
- `plot_ipa_dotplot.r`, `annotate_ipa().R`: downstream pathway visualization/annotation helpers

## Requirements

The workflow expects:

- R
- a Seurat object with normalized expression
- common metadata columns for:
  - group/condition
  - sample/replicate
  - cluster/cell type

Likely R package dependencies include:

- `Seurat`
- `Matrix`
- `dplyr`
- `ggplot2`
- `tibble`
- `ParetoTI`
- `reticulate`

`pareto_fit()` also checks for Python modules required by ParetoTI:

- `py_pcha`
- `geosketch`

## Quick Start

Because this repo is currently script-based, source the functions you need before running the workflow.

```r
source("check_input.r")
source("make_gene_universe.r")
source("select_hvgs.r")
source("pareto_fit2.r")
source("phca_extract_weights.r")
source("phca_dominant_archetype.r")
source("compute_signature_cor.r")
source("compute_occupancy.r")
source("compute_contrasts.r")
```

## Example Workflow

```r
# Validate inputs
cfg <- check_input(
  obj = sobj,
  assay = "RNA",
  group_col = "condition",
  sample_col = "orig.ident",
  cluster_col = "clusters"
)

# Fit archetype models
fit_res <- pareto_fit(
  obj = sobj,
  cfg = cfg,
  out_dir = "outputs",
  ks = 3:10,
  k_targets = c(4, 5),
  do_bootstrap = TRUE,
  bootstrap_n = 20,
  bootstrap_prop = 0.70
)

# Extract weights from a chosen fit
W5 <- pcha_extract_weights(
  fit = fit_res$fits[["k5"]],
  cell_ids = colnames(combined_no_IA),
  prefix = "A"
)

# Add weights to Seurat metadata
sobj@meta.data[colnames(W5)] <- W5[rownames(sobj@meta.data), , drop = FALSE]

# Add dominant archetype call
sobj$dominant_archetype <- phca_dominant_archetype(
  weights_mat = W5
)

# Compute correlation-based signatures
sig_k5 <- compute_signature_cor(
  obj = sobj,
  assay = "RNA",
  slot = "data",
  weights_mat = W5,
  method = "spearman",
  top_n = 100L
)

# Replicate-aware occupancy summaries
occ_k5 <- compute_occupancy(
  obj = sobj,
  weights_mat = W5,
  sample_col = "orig.ident",
  group_col = "condition",
  cluster_col = "clusters",
  min_cells = 10L
)

# Condition contrasts
contrast_k5 <- compute_contrast(
  occ_df = occ_k5$sample_cluster_occ,
  group_col = "condition",
  sample_col = "orig.ident",
  cluster_col = "clusters",
  group_levels = c("NEC", "NEC_HA"),
  test = "welch_t"
)
```

## Outputs

The fitting workflow writes organized outputs to a run directory, including:

- figures
- summary tables
- saved model objects
- frozen gene universe and HVG bundles

This makes runs easier to reproduce and compare across different values of `k`.

## Typical Use Cases

This repository is especially useful for analyses such as:

- identifying transcriptional archetypes in scRNA-seq
- comparing archetype occupancy between conditions
- studying cluster-specific redistribution of archetype programs
- deriving and ranking archetype-associated genes
- connecting archetype programs to pathway analyses and figure generation

## Status

This repository is under active development. The current codebase is functional for analysis workflows, but the repository is still organized primarily as sourced R scripts rather than a finalized installed package.

## Citation

If you use ArchetypeShift in your manuscript, please cite:
Future Publication Name

**Primary Contact**: Adam Wilson
**Lab**: Neo NEC Lab, University of Oklahoma Health Campus
**Email**: adam-wilson@ou.edu

## Data Availability

Due to storage constraints, raw and processed single-cell data are stored on the OU HPC system.
Metadata and example outputs can be shared upon request.

-
