# snRNA-seq Fetal–Pediatric Lung Pipeline

A modular, reproducible [`{targets}`](https://docs.ropensci.org/targets/) pipeline for single-nucleus RNA sequencing analysis of fetal and pediatric human lung tissue. The analysis characterizes developmental trajectories of alveolar epithelial cells using Seurat, Slingshot, tradeSeq, CellChat, and cell composition analyses.

> **Two workflows are available** — choose the one that fits your needs:
> - **Interactive Rmd** ([`snRNA_interactive_analysis.Rmd`](snRNA_interactive_analysis.Rmd)) — single file, run chunk-by-chunk in RStudio, plots appear inline. Best for exploration and parameter tuning.
> - **`{targets}` pipeline** (`_targets.R` + `R/`) — fully automated, cached, reproducible. Best for final runs and HPC.
>
> The original monolithic R Markdown workflow is preserved in [`archive/`](archive/).

---

## Scientific Background

This pipeline analyzes pre-processed Seurat objects from fetal and pediatric lung snRNA-seq data to investigate:

- Broad lineage composition (Epithelium, Mesenchyme, Immune, Endothelium)
- Alveolar epithelial maturation (Bud Tip → AT2 → Transitional → AT1)
- Pseudotime trajectories using Slingshot, with tradeSeq for pseudotime DE
- Cell–cell communication (BMP, TGFb, FGF pathways) using CellChat
- Cell type composition changes across fetal and pediatric age

**Two parallel epithelial/alveolar branches** are maintained throughout:

| Branch | Description |
|--------|-------------|
| `pre` (pre-merge) | Epithelial cells subset from individual fetal and pediatric objects before joint integration |
| `post` (post-merge) | Epithelial cells subset from the jointly integrated fetal–pediatric object |

Results from both branches are tracked and compared throughout the pipeline.

---

## Interactive R Markdown

[`snRNA_interactive_analysis.Rmd`](snRNA_interactive_analysis.Rmd) is a self-contained alternative to the `{targets}` pipeline. Open it in RStudio and run chunks one at a time — plots render inline.

### How it works

1. **Edit the CONFIG chunk** (Section 0) — every parameter lives there, nowhere else:
   - `PATH_fetal_raw` / `PATH_pediatric_raw` — input RDS paths
   - `CKPT_*` — checkpoint paths (set `NULL` to compute, or a file path to load a saved object and skip that step)
   - `SAVE_*` — where to write objects after computing (`NULL` = don't save)
   - Integration params: `DIMS`, `RESOLUTION`, `VARS_TO_REGRESS`, `PED_DOWNSAMPLE_N`
   - Annotation maps: `ANNOT_FETAL_LINEAGE`, `ANNOT_JOINT_CELLTYPE`, etc.
   - Trajectory, tradeSeq, CellChat, and plot settings

2. **Run the CONFIG chunk**, then step through sections in order.

3. **Load pre-computed objects** at any stage by setting the matching `CKPT_*` variable:

   | Variable | Object type |
   |----------|-------------|
   | `CKPT_fetal_integrated` | Seurat (fetal, integrated) |
   | `CKPT_joint_integrated` | Seurat (joint fetal–pediatric) |
   | `CKPT_epi_pre` / `CKPT_epi_post` | Seurat (epithelial subset) |
   | `CKPT_alve_pre` / `CKPT_alve_post` | Seurat (alveolar subset) |
   | `CKPT_sds_pre` / `CKPT_sds_post` | `SlingshotDataSet` |
   | `CKPT_sce_pre` / `CKPT_sce_post` | `SingleCellExperiment` (tradeSeq) |
   | `CKPT_cellchat_pre` / `CKPT_cellchat_post` | `CellChat` object |

---

## Repository Structure

```
.
├── snRNA_interactive_analysis.Rmd  # Interactive single-file workflow
├── _targets.R                  # Pipeline dependency graph (~50 targets)
├── config/
│   ├── paths.yaml              # Input/output paths — edit before running
│   ├── parameters.yaml         # Analysis parameters (dims, resolution, etc.)
│   ├── markers.yaml            # Marker gene sets by cell type
│   ├── annotations.yaml        # Cluster-to-celltype mapping tables
│   └── plotting.yaml           # Figure dimensions and color palettes
├── R/
│   ├── 00_utils.R              # Shared helpers: scale01, audit writers, validators
│   ├── 01_io.R                 # Load/save Seurat objects
│   ├── 02_qc.R                 # QC filtering
│   ├── 03_integration.R        # SCTransform + RPCA integration
│   ├── 04_annotation.R         # Cluster annotation from config maps
│   ├── 05_subsetting.R         # Epithelial and alveolar subsetting (pre + post)
│   ├── 06_alveolar_analysis.R  # UMAPs, dotplots, DEG
│   ├── 07_trajectory_slingshot.R # Slingshot trajectory inference
│   ├── 08_tradeseq.R           # GAM fitting and pseudotime DE tests
│   ├── 09_cellchat.R           # CellChat cell–cell communication
│   ├── 10_composition.R        # Cell composition tables and plots
│   ├── 11_plotting.R           # UMAP panels, DEG heatmaps
│   └── 12_audit.R              # Per-module audit file writers
├── reports/
│   ├── final_report.qmd        # Final HTML report (reads pipeline outputs)
│   └── module_audit_report.qmd # Per-module audit summary report
├── scripts/
│   ├── run_pipeline.R          # Run full pipeline
│   ├── run_module.R            # Run a single named module
│   ├── inspect_targets.R       # Print target status and outdated targets
│   └── clean_outputs.R         # Safely clean output directories
├── tests/
│   ├── testthat.R
│   └── testthat/               # 60 tests for config and utility functions
├── archive/
│   └── snRNA_analysis_reorganized_original.Rmd  # Original workflow — preserved
└── outputs/                    # Generated files (gitignored, recreated by pipeline)
    ├── rds/  plots/  tables/  audit/  logs/  reports/
```

---

## Quickstart

### Option A — Interactive Rmd (exploration / parameter tuning)

1. Open `snRNA_interactive_analysis.Rmd` in RStudio.
2. Edit the **CONFIG chunk** (Section 0): set `PATH_fetal_raw`, `PATH_pediatric_raw`, and any `CKPT_*` paths for pre-computed objects you already have.
3. Run the CONFIG chunk, then step through sections with **Ctrl/Cmd + Shift + Enter**. Plots appear inline.

### Option B — `{targets}` pipeline (automated / reproducible)

**1. Set your data paths** — edit `config/paths.yaml`:

```yaml
input:
  fetal_rds: "/path/to/TF_fetal_nuc.rds"
  pediatric_rds: "/path/to/TF_ped_nuc.rds"
```

**2. Run the full pipeline:**

```r
library(targets)
tar_make()
```

or from the terminal:

```bash
Rscript scripts/run_pipeline.R
```

**3. Run a single module:**

```bash
Rscript scripts/run_module.R trajectory
Rscript scripts/run_module.R cellchat
Rscript scripts/run_module.R final_report
```

Available modules: `config`, `input`, `individual_integration`, `annotation`, `joint_integration`, `epithelial`, `alveolar`, `trajectory`, `tradeseq`, `cellchat`, `composition`, `final_plots`, `final_tables`, `audit_report`, `final_report`

**4. Check pipeline status:**

```r
library(targets)
tar_outdated()     # which targets need rerunning
tar_visnetwork()   # interactive dependency graph
```

---

## Input Data Requirements

Two pre-processed Seurat objects (Seurat v5) with these metadata columns:

| Column | Description |
|--------|-------------|
| `batch` | Sample/batch identifier |
| `sex` | Donor sex |
| `age` | Age label (e.g., `"103 days"`) |
| `age_days` | Age in days (numeric) |
| `group` | `"fetal"` or `"pediatric"` |
| `TF_annotation` | Cell type annotation |
| `TF_annotation_lvl1/2/3` | Hierarchical annotations |
| `percent.mt` | Mitochondrial read fraction |

---

## Configuration

All analysis parameters live in `config/` — no need to edit R source files.

| File | What to edit |
|------|-------------|
| `config/paths.yaml` | Data paths, output directories |
| `config/parameters.yaml` | PCA dims, resolution, trajectory start, tradeSeq knots/cores |
| `config/markers.yaml` | Marker gene sets for cell type identification |
| `config/annotations.yaml` | Cluster-to-celltype maps (update after re-clustering) |
| `config/plotting.yaml` | Figure sizes, DPI, color palettes |

---

## Outputs

```
outputs/
├── plots/
│   ├── integration/   annotation/   epithelial/   alveolar/
│   ├── trajectory/    tradeseq/     cellchat/     composition/
├── tables/
│   ├── markers/       trajectory/   tradeseq/
│   ├── cellchat/      composition/
├── audit/             # Per-module YAML + CSV + Markdown audit files
└── reports/           # Rendered HTML reports
```

Filenames use stable branch-labelled conventions:

```
pre_alveolar_umap_panel.png
post_alveolar_umap_panel.png
pre_tradeseq_startvsend_results.csv
post_tradeseq_association_results.csv
```

---

## Audit System

Every major pipeline stage writes audit records to `outputs/audit/`:

- `<module>_audit.yaml` — machine-readable
- `<module>_audit.csv` — spreadsheet-importable
- `<module>_audit.md` — human-readable table

Each audit records: timestamp, cell counts before/after, metadata columns, assays, reductions, parameters used, filtering decisions, and warnings.

`outputs/audit/code_fix_log.yaml` documents three bugs fixed from the original Rmd during refactoring.

Generate the consolidated audit HTML report:

```bash
Rscript scripts/run_module.R audit_report
```

---

## Known Bugs Fixed from Original Rmd

| Issue | Location | Fix |
|-------|----------|-----|
| Contaminant gene filters in post-merge epithelial subsetting were bare expressions — never assigned or applied | `R/05_subsetting.R` | Wrapped in `filter_contaminants()` which explicitly returns the filtered object |
| Undefined variable `curves` in tradeSeq `fitGAM` call | `R/08_tradeseq.R` | `run_fitgam()` takes `sds_obj` as an explicit argument |
| Undefined variable `sce2` in `startVsEndTest` | `R/08_tradeseq.R` | `run_startvsend_test()` takes `sce` as an explicit argument |

---

## Memory Requirements

**≥ 40 GB RAM** recommended. The most memory-intensive steps:

- Joint fetal–pediatric integration (`joint_integrated`)
- tradeSeq GAM fitting (`pre_tradeseq`, `post_tradeseq`)
- CellChat pipeline (`pre_cellchat`, `post_cellchat`)

Set globally:

```r
options(future.globals.maxSize = 40 * 1024^3)
```

This is already set in `scripts/run_pipeline.R`.

---

## Package Setup

```r
# Bioconductor packages
BiocManager::install(c(
  "Seurat", "SeuratWrappers", "sctransform",
  "slingshot", "SingleCellExperiment", "tradeSeq", "BiocParallel",
  "ComplexHeatmap"
))

# CRAN packages
install.packages(c(
  "targets", "tarchetypes",
  "ggplot2", "patchwork", "viridis", "ggridges", "ggrepel",
  "ggnewscale", "cowplot", "circlize", "pheatmap", "RColorBrewer",
  "dplyr", "tidyr", "tibble", "tidytext",
  "presto", "writexl",
  "CellChat", "NMF", "ggalluvial",
  "yaml", "R.utils"
))
```

If using `renv`, restore from the lockfile (if present):

```r
renv::restore()
```

---

## Running Tests

```r
library(testthat)
test_dir("tests/testthat")
```

Or from the terminal:

```bash
Rscript tests/testthat.R
```

60 tests cover config file loading, required key validation, utility functions, and audit file round-trips.

---

## Troubleshooting

**Outdated targets after editing R code:** Run `tar_outdated()` — targets will rebuild automatically on the next `tar_make()`.

**Cluster annotation mismatch after re-clustering:** Update the relevant map in `config/annotations.yaml` (e.g., `joint_celltype`, `alveolar_pre_celltype2`).

**`integrated.dr` reduction not found:** IntegrateLayers failed silently — check the integration audit file. The pipeline falls back to `pca`.

**CellChat fails with too few cells:** Lower `cellchat.min_cells` in `config/parameters.yaml`.

**tradeSeq is slow:** Reduce `tradeseq.cores` or `tradeseq.nknots` in `config/parameters.yaml`. Consider running on HPC.

---

## Original R Markdown

The original single-file workflow is archived at:

```
archive/snRNA_analysis_reorganized_original.Rmd
```

It is never deleted or modified. Use it as the scientific reference for all pipeline logic.

---

## For AI Coding Agents

See [`AGENTS.md`](AGENTS.md) for detailed instructions on how to safely modify targets, R functions, parameters, and audit outputs in this repository.
