# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a single-file R Markdown pipeline for analyzing snRNA-seq (single-nuclei RNA sequencing) data from fetal and pediatric lung tissue. The analysis investigates developmental trajectories of lung alveolar cells.

**Primary file:** `snRNA_analysis_reorganized.Rmd` (~2,274 lines)

## Running the Pipeline

```r
# Render full HTML report
rmarkdown::render("snRNA_analysis_reorganized.Rmd")
```

Or run interactively chunk-by-chunk in RStudio. The pipeline saves RDS checkpoints at key stages so expensive steps can be skipped on re-runs.

**Memory requirement:** ≥40GB RAM (`options(future.globals.maxSize = 40 * 1024^3)` is set globally).

## Input Data

Pre-processed Seurat objects expected at:
- `.../Sequencing data/20250528/TF_fetal_nuc.rds`
- `.../Sequencing data/20250528/TF_ped_nuc.rds`

Key metadata columns on input objects: `batch`, `sex`, `age`, `age_days`, `group` (fetal/pediatric), `TF_annotation` (and `_lvl1/2/3`), `percent.mt`.

## Pipeline Architecture

The Rmd is organized into sequential sections with RDS checkpoints between major stages:

| Stage | Lines | Output RDS |
|-------|-------|------------|
| Setup & libraries | 1–105 | — |
| Data input | 109–141 | — |
| Individual integration (fetal + ped) | 144–256 | `TF_fetal_nuc_cluster.rds`, `TF_pediatric_nuc_subset_cluster.rds` |
| Lineage annotation | 259–343 | — |
| Joint fetal–pediatric integration | 346–545 | `TF_fetal_pediatric_cluster.rds` |
| Epithelial subsetting & re-integration | 546–803 | — |
| Alveolar lineage analysis | 798–1235 | — |
| Slingshot trajectory inference | 1271–1740 | — |
| CellChat cell-cell communication | 1815–2088 | — |
| Cell composition analysis | 2089–2254 | — |
| Save all objects | 2258–2274 | Multiple |

## Key Methods

- **Integration:** SCTransform normalization + RPCA batch correction; 30 PCA dims, resolution 0.5
- **Trajectory:** Slingshot from Bud Tip progenitors; weighted pseudotime across 4 lineages; tradeSeq GAMs for pseudotime DE
- **Communication:** CellChat with human L-R database; pathway-level (BMP, TGFb, FGF) analysis
- **Annotation:** Manual cluster-to-celltype mapping; 4 broad lineages (Epithelium, Mesenchyme, Immune, Endothelium); fine types include AT1, AT2, LAP, Transitional, Bud Tip

## Epithelial Subsetting — Two Parallel Approaches

The pipeline maintains **two variants** throughout sections 6–9:
- **pre-merge** (`TF_fetal_pediatric_epi_pre`, `pre_alve`): subset from individual objects, then merge and integrate
- **post-merge** (`TF_fetal_pediatric_epi_post`, `post_alve`): subset from joint integrated object

Both are carried through trajectory and CellChat analyses for comparison.

## Known Code Issues (FLAGGED comments)

- **~Line 641–644:** Contaminant gene filters in post-merge approach were written as bare expressions without assignment — filters were not applied as written in original code
- **~Line 1621:** Original code referenced undefined variable `curves`; should be `getCurves(sds_BT)`
- **~Line 1656:** Original code referenced `sce2`; should be `sce_post`

## Core R Packages

```r
# Single-cell
Seurat, SeuratWrappers, sctransform, SingleCellExperiment

# Trajectory
slingshot, tradeSeq, BiocParallel

# Communication
CellChat, NMF, ggalluvial

# Visualization
ggplot2, patchwork, pheatmap, ComplexHeatmap, viridis, ggrepel, cowplot, circlize

# DE / markers
presto, writexl

# Utilities
dplyr, tidyverse, tibble, R.utils
```
