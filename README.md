# ArchetypeShift

ArchetypeShift is an R package for fitting, annotating, and comparing transcriptional archetype programs in single-cell RNA-seq data, with particular emphasis on biological interpretation through correlation-based gene signatures, IPA canonical pathways, upstream regulators, and cluster/state localization. It is designed for Seurat-based workflows and supports reproducible preprocessing, ParetoTI/PCHA model fitting, archetype annotation, and downstream condition- and cluster-level comparisons. The package is especially useful when the goal is not only to identify archetypes, but also to interpret them as biologically meaningful programs across cell states and experimental conditions. The package was developed for epithelial archetype analysis in scRNA-seq, but the core functions are general enough to be adapted to other cell systems.

## Installation

You can install ArchetypeShift with:

```r
devtools::install_github("Neo-NEC-Lab/ArchetypeShift")
library(ArchetypeShift)
```

## Core Workflow

The package supports the following analysis path:

1. Validate the Seurat object and resolve metadata columns with `check_input()`
2. Define a stable gene universe with `make_gene_universe()`
3. Select HVGs with `select_hvgs()`
4. Fit ParetoTI archetype models with `pareto_fit()`
5. Extract weights with `pcha_extract_weights()`
6. Add dominant archetype assignments with `phca_dominant_archetype()`
7. Compute gene signatures with `compute_signature_cor()`
8. Annotate acrchetypes with IPA using `annotate_ipa()`
9. Summarize archetype occupancy with `compute_occupancy()`
10. Test replicate-aware contrasts with `compute_contrast()`
11. Generate downstream occupancy, delta, QC, and enrichment plots

## Repository Contents

Some of the main functions in this package are:

- `check_input.r`: validates Seurat objects and required metadata
- `detect_group_col.r`, `detect_meta_cols.r`, `detect_sample_col.r`, `detect_cluster_col.r`: helper functions for metadata auto-detection
- `make_gene_universe.r`: defines a frozen background gene universe
- `select_hvgs.r`: selects HVGs for PCA
- `pareto_fit.r`: main ParetoTI/PCHA fitting workflow with k-scan support
- `phca_extract_weights.r`: extracts a standardized cells x k archetype weight matrix
- `phca_dominant_archetype.r`: assigns dominant archetypes per cell
- `phca_scan_metrics.r`: summarizes k-scan metrics
- `compute_signature_cor.r`: computes correlation-based archetype signatures
- `read_ipa_csv.r`, `annotate_ipa.r`: import and annotate archetypes using IPA pathways and upstream regulators
- `compute_occupancy.r`: replicate-aware sample and sample-by-cluster occupancy summaries
- `compute_contrasts.r`: replicate-aware condition contrasts for archetype occupancy
- `compute_pcha_delta_heatmap_df.r`: helper for replicate-aware delta occupancy heatmaps
- `compute_pcha_weight_qc_from_W.r`: computes per-cell QC summaries from archetype weights
- `plot_pcha_delta_heatmap.r`, `plot_pcha_occupancy_heatmap.r`, `plot_pcha_occupancy_stacked.r`, `plot_pcha_qc_hist.r`, `plot_ipa_dotplot.r`: plotting utilities

## Requirements

The workflow expects:

- R
- a Seurat object with normalized expression
- common metadata columns for:
  - group/condition
  - sample/replicate
  - cluster/cell type

Some R package dependencies include:
```r
library(Seurat)
libary(Matrix)
libary(dplyr)
library(ggplot2)
library(tibble)
library(reticulate)
libary(ParetoTI)
```

pareto_fit() also checks for Python modules required by ParetoTI:

- `py_pcha`
- `geosketch`

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
  cell_ids = colnames(sobj),
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

# Read IPA exports
IPA_pathways <- read_ipa_csv("IPA_pathways.csv")
IPA_upstream <- read_ipa_csv("IPA_upstream.csv")

# Annotate archetypes with IPA pathways and upstream regulators
ipa_res_k5 <- annotate_ipa(
  signatures = sig_k5$signatures,
  pathways_df = IPA_pathways,
  upstream_df = IPA_upstream,
  universe = NULL,
  by_cluster = TRUE,
  method = "fisher",
  min_overlap = 2L,
  p_adjust = "BH"
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

## Example Uses

This package is useful for analyses such as:

- Identifying transcriptional archetypes in scRNA-seq
- Comparing archetype occupancy between conditions
- Studying cluster-specific redistribution of archetype programs
- Deriving and ranking archetype-associated genes
- Connecting archetype programs to pathway analyses and figure generation

## Citation

If you use ArchetypeShift in your manuscript, please cite:  
Future Publication Name

**Primary Contact**: Adam Wilson   
**Lab**: Neo NEC Lab, University of Oklahoma Health Campus  
**Email**: adam-wilson@ou.edu  

## Data Availability

De-identified raw FASTQ files and processed Cell Ranger outputs generated in this study will be publicly available through GEO. The source code and example workflow are available in this repository.
