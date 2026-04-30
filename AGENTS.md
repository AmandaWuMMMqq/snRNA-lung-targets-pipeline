# AGENTS.md

This file gives durable instructions for AI coding agents working on this repository.

## Repository purpose

This repository contains a modular, reproducible snRNA-seq analysis pipeline for fetal and pediatric lung tissue. The original analysis was developed as a large R Markdown workflow and is being refactored into a `{targets}` pipeline with modular R functions, YAML parameters, audit outputs, figures, tables, and final reports.

The scientific goal is to analyze lung developmental trajectories, especially alveolar epithelial maturation, using Seurat, Slingshot, tradeSeq, CellChat, and cell composition analyses.

## Core design principle

Keep computation, configuration, outputs, and interpretation separate.

- `R/`: reusable functions only
- `config/`: user-editable parameters, paths, markers, annotations, plotting choices
- `_targets.R`: pipeline dependency graph
- `outputs/`: generated objects, plots, tables, reports, and audit files
- `reports/`: report templates that read completed outputs
- `archive/`: preserved original R Markdown workflow

Do not turn the pipeline back into one giant script.

## Scientific context

The pipeline analyzes pre-processed Seurat objects for fetal and pediatric lung snRNA-seq data. Expected metadata columns include:

- `batch`
- `sex`
- `age`
- `age_days`
- `group`
- `TF_annotation`
- `TF_annotation_lvl1`
- `TF_annotation_lvl2`
- `TF_annotation_lvl3`
- `percent.mt`

The pipeline includes two parallel epithelial/alveolar workflows:

1. `pre` / pre-merge branch: subset epithelial cells from individual fetal and pediatric objects before merging and integration.
2. `post` / post-merge branch: subset epithelial cells from the jointly integrated fetal-pediatric object.

Preserve this distinction. Do not collapse these branches unless the user explicitly asks.

## Rules for modifying `_targets.R`

1. Every major analysis stage should be represented as one or more explicit targets.
2. Target names should be stable, descriptive, and lowercase.
3. Use branch labels consistently: `pre_*` and `post_*`.
4. Targets should call functions from `R/`; do not place long analysis code directly inside `_targets.R`.
5. Reports should be rendered through `tarchetypes::tar_render()` when possible.
6. Heavy objects should be saved as target outputs, not repeatedly regenerated inside reports.
7. Avoid hidden dependencies on files not declared in targets or config.
8. If you add a new output folder, update config and README.

## Rules for modifying R functions

1. Functions must have explicit inputs and outputs.
2. Do not rely on hidden global variables.
3. Do not assume a Seurat object has a specific assay, reduction, or metadata column without checking.
4. Use validation helpers before expensive operations.
5. Return objects rather than only saving files.
6. If a function saves plots or tables, the output path should be passed as an argument or read from config.
7. Emit informative messages for long operations.
8. Prefer small, auditable functions over large multi-stage functions.

## Rules for parameters

1. Analysis thresholds, file paths, dimensions, resolutions, marker lists, and plotting choices should live in YAML files under `config/`.
2. Do not hard-code user-specific absolute paths inside R functions.
3. If a new parameter is introduced, add it to the relevant YAML file and document it in README.
4. Use clear defaults, but avoid silently changing scientific assumptions.
5. If a parameter changes a biological interpretation, add an audit note.

## Rules for auditability

Every major module should generate an audit output that includes:

- module name
- timestamp
- input target or object name
- output target or object name
- number of cells before and after
- number of genes/features when relevant
- metadata columns present
- assays present
- reductions present
- parameter values used
- filtering decisions
- warnings or assumptions
- generated output files

Audit files should be saved under `outputs/audit/` as YAML, CSV, or Markdown when appropriate.

The consolidated audit report should make it easy to understand what changed at each stage without rerunning the full analysis.

## Rules for plots and tables

1. Save plots under `outputs/plots/<module>/`.
2. Save tables under `outputs/tables/<module>/`.
3. Use stable, descriptive filenames.
4. Include branch labels in filenames when applicable, such as `pre_alveolar_umap_by_celltype.png` and `post_alveolar_umap_by_celltype.png`.
5. Do not overwrite important outputs with generic names like `plot.png` or `markers.csv`.
6. Reports should load saved plots/tables or target objects rather than recomputing heavy steps.

## Rules for reports

1. Reports should focus on presentation and interpretation, not heavy computation.
2. Clearly separate computational results from biological interpretation.
3. Do not overstate causal conclusions from observational snRNA-seq data.
4. Include session information and parameter summaries.
5. Include audit summaries for reproducibility.

## Known original-code issues to preserve/fix/document

The original R Markdown workflow had flagged issues that must be fixed or documented:

1. Contaminant gene filtering in the post-merge epithelial approach was written as bare expressions without assignment. Ensure filtering is explicit and assigned.
2. An undefined `curves` variable should be replaced with the correct Slingshot-derived object, such as `getCurves(sds_BT)`, depending on context.
3. An undefined `sce2` variable should be replaced with the correct SingleCellExperiment object, likely `sce_post`, depending on context.
4. Search for additional undefined variables, silent overwrites, or reliance on stale global objects.

## Testing expectations

Before finalizing changes, run lightweight checks when possible:

```r
source("R/00_utils.R")
testthat::test_dir("tests/testthat")
targets::tar_manifest()
targets::tar_validate()
```

Do not require the full raw data to run lightweight tests. Use small mock objects where appropriate.

## What not to do

Do not:

- delete or overwrite the original R Markdown file without preserving it in `archive/`
- commit raw data, large RDS files, or generated target stores
- hard-code local absolute paths in R functions
- mix heavy analysis into report files
- collapse pre-merge and post-merge branches
- silently change biological thresholds or marker definitions
- remove audit outputs to simplify code
- make conclusions stronger than the data support

## Preferred workflow for future agents

1. Read `README.md`, `_targets.R`, and relevant YAML config files first.
2. Inspect target dependencies before editing.
3. Make the smallest safe change.
4. Add or update tests for new helpers.
5. Update audit logic if the change affects analysis outputs.
6. Update README if user-facing behavior changes.
7. Summarize exactly what changed and what still requires manual verification.
