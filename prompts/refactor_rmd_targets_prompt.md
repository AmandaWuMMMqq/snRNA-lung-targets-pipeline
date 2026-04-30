# Claude Code Prompt: Refactor R Markdown snRNA-seq Analysis into a Modular `{targets}` Pipeline

You are working in a GitHub repository that currently contains a large single-file R Markdown analysis pipeline:

- Primary file: `snRNA_analysis_reorganized.Rmd`
- Approximate size: ~2,274 lines
- Analysis type: fetal and pediatric lung snRNA-seq analysis using Seurat, Slingshot, tradeSeq, CellChat, and cell composition analyses

Your task is to refactor this R Markdown workflow into a reproducible, modular, auditable `{targets}`-based analysis pipeline while preserving the scientific logic and outputs of the original analysis.

## High-level goal

Convert the existing monolithic R Markdown file into a clean GitHub repository with:

1. A `{targets}` pipeline for reproducible execution.
2. Modular R functions grouped by analysis stage.
3. User-editable parameter files so I can change thresholds, paths, resolutions, dimensions, and plotting options without editing source code.
4. Automatically generated audit reports after each major module.
5. Automatically generated plots and tables saved in organized output folders.
6. A final Quarto or R Markdown report that reads pipeline outputs and summarizes the analysis.
7. An `AGENTS.md` file that explains how future AI agents should modify, audit, and extend the repository.
8. A complete GitHub-ready repository structure with README, environment notes, and reproducibility instructions.

Do not delete the original R Markdown file. Move or copy it into an archive folder for reference.

---

## Important scientific context

The current analysis investigates developmental trajectories of lung alveolar cells from fetal and pediatric snRNA-seq data.

Input Seurat objects are expected at paths similar to:

```r
.../Sequencing data/20250528/TF_fetal_nuc.rds
.../Sequencing data/20250528/TF_ped_nuc.rds
```

Key metadata columns expected on input objects:

```r
batch
sex
age
age_days
group
TF_annotation
TF_annotation_lvl1
TF_annotation_lvl2
TF_annotation_lvl3
percent.mt
```

The original pipeline includes these major stages:

1. Setup and libraries
2. Data input
3. Individual fetal and pediatric integration
4. Lineage annotation
5. Joint fetal-pediatric integration
6. Epithelial subsetting and re-integration
7. Alveolar lineage analysis
8. Slingshot trajectory inference
9. tradeSeq pseudotime differential expression
10. CellChat cell-cell communication
11. Cell composition analysis
12. Save all final objects

The current Rmd maintains two parallel epithelial/alveolar analysis strategies:

- `pre-merge`: subset from individual fetal/pediatric objects first, then merge and integrate
- `post-merge`: subset from the joint integrated object

Preserve both branches in the new pipeline and make the comparison explicit.

---

## Known issues to fix during refactor

The existing Rmd has flagged code issues. Fix these during modularization and document the fix in an audit note.

1. Contaminant gene filters in the post-merge epithelial approach were written as bare expressions without assignment. These filters were not applied as intended. Refactor into an explicit function and make sure filtered objects are assigned and saved.
2. A section referenced undefined variable `curves`; this should be corrected to use `getCurves(sds_BT)` or the correct Slingshot object-derived curves.
3. A section referenced `sce2`; this should be corrected to `sce_post` or the correct SingleCellExperiment object.
4. Audit for other undefined objects, silent overwrites, inconsistent object names, or chunks that rely on hidden global state.

---

## Required repository structure

Create or update the repository to follow this structure:

```text
.
├── README.md
├── AGENTS.md
├── _targets.R
├── _targets.yaml                  # optional but preferred if useful
├── renv.lock                      # if renv is initialized or already present
├── config/
│   ├── paths.yaml
│   ├── parameters.yaml
│   ├── markers.yaml
│   ├── annotations.yaml
│   └── plotting.yaml
├── R/
│   ├── 00_utils.R
│   ├── 01_io.R
│   ├── 02_qc.R
│   ├── 03_integration.R
│   ├── 04_annotation.R
│   ├── 05_subsetting.R
│   ├── 06_alveolar_analysis.R
│   ├── 07_trajectory_slingshot.R
│   ├── 08_tradeseq.R
│   ├── 09_cellchat.R
│   ├── 10_composition.R
│   ├── 11_plotting.R
│   └── 12_audit.R
├── reports/
│   ├── final_report.qmd
│   ├── module_audit_report.qmd
│   └── templates/
├── scripts/
│   ├── run_pipeline.R
│   ├── run_module.R
│   ├── inspect_targets.R
│   └── clean_outputs.R
├── archive/
│   └── snRNA_analysis_reorganized_original.Rmd
├── outputs/
│   ├── rds/
│   ├── tables/
│   ├── plots/
│   ├── audit/
│   ├── logs/
│   └── reports/
└── tests/
    ├── testthat.R
    └── testthat/
```

If some files already exist, update them instead of overwriting useful content.

---

## Parameterization requirements

Create YAML config files so that I can edit key analysis choices without touching the R source files.

### `config/paths.yaml`

Include at least:

```yaml
input:
  fetal_rds: "PATH/TO/TF_fetal_nuc.rds"
  pediatric_rds: "PATH/TO/TF_ped_nuc.rds"

output:
  root: "outputs"
  rds: "outputs/rds"
  plots: "outputs/plots"
  tables: "outputs/tables"
  audit: "outputs/audit"
  logs: "outputs/logs"
  reports: "outputs/reports"
```

### `config/parameters.yaml`

Include at least:

```yaml
memory:
  future_globals_max_size_gb: 40

seurat:
  dims: 30
  resolution: 0.5
  integration_method: "RPCAIntegration"
  normalization_method: "SCT"
  variable_features_n: 3000

qc:
  min_features: null
  max_percent_mt: null
  min_cells: null

subsetting:
  epithelial_lineage_labels:
    - "Epithelium"
  alveolar_labels:
    - "AT1"
    - "AT2"
    - "LAP"
    - "Transitional"
    - "Bud Tip"
  contaminant_genes:
    - "PTPRC"
    - "COL1A1"
    - "PECAM1"

trajectory:
  start_cluster: "Bud Tip"
  slingshot_reduction: "UMAP"
  n_lineages_to_keep: null
  tradeSeq_nknots: 6
  tradeSeq_cores: 4

cellchat:
  database: "human"
  pathways_of_interest:
    - "BMP"
    - "TGFb"
    - "FGF"

plots:
  width: 8
  height: 6
  dpi: 300
```

### `config/markers.yaml`

Include marker sets for broad lineages and lung epithelial subtypes. At minimum include placeholders for:

```yaml
lineage_markers:
  epithelium:
    - EPCAM
    - KRT8
  mesenchyme:
    - COL1A1
    - DCN
  immune:
    - PTPRC
  endothelium:
    - PECAM1
    - VWF

alveolar_markers:
  AT1:
    - AGER
    - PDPN
    - CAV1
  AT2:
    - SFTPC
    - SFTPA1
    - SFTPA2
    - ABCA3
    - NAPSA
  bud_tip:
    - SOX9
    - ID2
```

### `config/annotations.yaml`

Move manual cluster-to-celltype mappings here where possible. The code should read these mappings rather than hard-coding them deep in scripts.

### `config/plotting.yaml`

Store color palettes, figure dimensions, and plot-specific options here.

---

## Functional module requirements

Each R file in `R/` should define functions, not execute major analysis at import time.

Good pattern:

```r
run_qc <- function(seurat_obj, params, audit_dir = NULL) {
  # validate inputs
  # run QC
  # return object and write audit outputs
}
```

Avoid scripts that rely on many hidden global variables.

Each function should:

1. Take explicit inputs.
2. Read parameters from a `params` list or function arguments.
3. Return explicit outputs.
4. Save plots/tables only through controlled output paths.
5. Emit useful messages.
6. Write an audit summary when appropriate.
7. Fail with informative errors when required metadata or assays are missing.

---

## `{targets}` pipeline requirements

Create `_targets.R` that:

1. Loads required packages.
2. Sources functions from `R/`.
3. Reads YAML config files.
4. Defines targets for each major analysis stage.
5. Saves intermediate RDS objects and audit summaries.
6. Keeps the pre-merge and post-merge epithelial branches separate.
7. Generates final tables, plots, and reports.

Expected target groups:

```r
# config
paths
params
markers
annotations
plotting

# input
fetal_raw
pediatric_raw
input_audit

# individual processing
fetal_integrated
pediatric_integrated
individual_integration_audit

# annotation
fetal_annotated
pediatric_annotated
annotation_audit

# joint integration
joint_integrated
joint_integration_audit

# epithelial branches
pre_epi
post_epi
pre_epi_integrated
post_epi_integrated
epithelial_branch_audit

# alveolar branches
pre_alveolar
post_alveolar
alveolar_audit

# trajectory
pre_slingshot
post_slingshot
trajectory_audit

# tradeSeq
pre_tradeseq
post_tradeseq
tradeseq_audit

# cellchat
pre_cellchat
post_cellchat
cellchat_audit

# composition
composition_tables
composition_plots
composition_audit

# final outputs
final_plots
final_tables
module_audit_report
final_report
```

Use `tar_target()` and where appropriate `tar_render()` from `tarchetypes` for reports.

Use sensible target names that make it easy to rerun only one module.

---

## Audit report requirements

For every major module, generate a machine-readable and human-readable audit output.

Each audit should include:

1. Module name
2. Timestamp
3. Input object names or target names
4. Output object names or target names
5. Number of cells/nuclei before and after
6. Number of genes/features when relevant
7. Metadata columns present
8. Assays present
9. Reductions present
10. Parameter values used
11. Warnings or assumptions
12. Any filtering performed
13. Session/package versions when useful

Save audit outputs as:

```text
outputs/audit/<module_name>_audit.yaml
outputs/audit/<module_name>_audit.csv
outputs/audit/<module_name>_audit.md
```

Create a consolidated audit report:

```text
outputs/reports/module_audit_report.html
```

The audit report should make it easy to answer:

- What was run?
- Which parameters were used?
- How many cells were retained or removed?
- Which outputs were generated?
- Were there warnings or possible issues?
- What changed compared with the previous stage?

---

## Plot and table output requirements

Save plots in organized folders:

```text
outputs/plots/qc/
outputs/plots/integration/
outputs/plots/annotation/
outputs/plots/epithelial/
outputs/plots/alveolar/
outputs/plots/trajectory/
outputs/plots/tradeseq/
outputs/plots/cellchat/
outputs/plots/composition/
```

Save tables in organized folders:

```text
outputs/tables/qc/
outputs/tables/markers/
outputs/tables/trajectory/
outputs/tables/tradeseq/
outputs/tables/cellchat/
outputs/tables/composition/
```

Use stable filenames that include branch labels where needed, e.g.:

```text
pre_alveolar_umap_by_celltype.png
post_alveolar_umap_by_celltype.png
pre_slingshot_lineages.png
post_slingshot_lineages.png
pre_tradeseq_association_results.csv
post_tradeseq_association_results.csv
```

---

## Final report requirements

Create `reports/final_report.qmd` or `reports/final_report.Rmd` that reads only completed target outputs and does not rerun heavy analysis directly inside the report.

The final report should include:

1. Analysis overview
2. Input data summary
3. Individual fetal/pediatric processing summary
4. Joint integration summary
5. Pre-merge vs post-merge epithelial branch comparison
6. Alveolar subset comparison
7. Slingshot trajectory results
8. tradeSeq pseudotime DE results
9. CellChat results focused on BMP, TGFb, FGF, and any other configured pathways
10. Cell composition analysis
11. Key figures
12. Key tables
13. Audit summary
14. Session info

The report should distinguish clearly between:

- computational result
- biological interpretation
- assumption or possible artifact

Do not overstate biological conclusions.

---

## `scripts/` requirements

Create helper scripts:

### `scripts/run_pipeline.R`

Runs the full pipeline:

```r
targets::tar_make()
```

### `scripts/run_module.R`

Allows me to run a selected target or group of targets from the command line. For example:

```bash
Rscript scripts/run_module.R trajectory
Rscript scripts/run_module.R cellchat
Rscript scripts/run_module.R final_report
```

Implement a simple mapping from module names to target names.

### `scripts/inspect_targets.R`

Prints target status, outdated targets, and dependency graph summary.

### `scripts/clean_outputs.R`

Safely cleans selected outputs. It should not delete raw data or archived original Rmd.

---

## README requirements

Update or create `README.md` with:

1. Project overview
2. Repository structure
3. Input data requirements
4. How to configure paths and parameters
5. How to run the full pipeline
6. How to run a single module
7. How to generate final reports
8. How to inspect target status
9. Output organization
10. Troubleshooting common issues
11. Notes on memory requirements, currently expected to be around 40 GB RAM
12. Notes on package setup using `renv` if available
13. A statement that the original Rmd is archived and preserved

---

## AGENTS.md requirements

Create an `AGENTS.md` file with durable instructions for future AI coding agents.

It should explain:

1. The repository purpose
2. The scientific context
3. The pipeline architecture
4. Rules for modifying targets
5. Rules for modifying R functions
6. Rules for adding parameters
7. Rules for auditability
8. Rules for plots and reports
9. Rules for avoiding hidden global state
10. Rules for preserving scientific validity
11. How to test changes
12. What not to do

Make this file practical and specific, not generic.

---

## Testing and validation requirements

Add basic tests using `testthat` where practical.

At minimum test:

1. Config files can be loaded.
2. Required config keys exist.
3. Required metadata columns are checked properly.
4. Audit helper functions produce expected fields.
5. Plot saving helpers create output directories.
6. Branch labels are handled consistently for `pre` and `post` workflows.

Add lightweight validation functions in `R/00_utils.R`, such as:

```r
check_required_metadata <- function(obj, required_cols) { ... }
check_required_assays <- function(obj, required_assays) { ... }
check_output_dirs <- function(paths) { ... }
```

---

## Coding style requirements

Use readable, explicit R code.

Prefer:

```r
run_integration <- function(obj, params) {
  dims <- seq_len(params$seurat$dims)
  resolution <- params$seurat$resolution
  ...
}
```

Avoid:

```r
# implicit globals
obj <- SCTransform(obj)
```

Use consistent naming:

```text
pre_epi
post_epi
pre_alveolar
post_alveolar
pre_slingshot
post_slingshot
```

Do not introduce unnecessary abstractions that make the biological analysis harder to inspect.

---

## GitHub setup requirements

Make the repository GitHub-ready.

Add or update:

```text
.gitignore
README.md
AGENTS.md
LICENSE        # only if one already exists or ask me before choosing one
```

`.gitignore` should exclude heavy/generated outputs and local files, including:

```gitignore
.Rproj.user/
.Rhistory
.RData
.Ruserdata
_targets/
outputs/rds/
outputs/plots/
outputs/tables/
outputs/audit/
outputs/logs/
outputs/reports/
*.rds
*.RDS
*.h5seurat
*.h5ad
.DS_Store
```

Do not commit raw data or large generated RDS files.

Do preserve small config files, source code, reports templates, and documentation.

---

## Implementation strategy

Please proceed in this order:

1. Inspect the existing repository and locate `snRNA_analysis_reorganized.Rmd`.
2. Create a new branch if appropriate, e.g. `refactor-targets-pipeline`.
3. Archive the original Rmd into `archive/` while keeping a copy in place if needed.
4. Extract repeated or logically grouped code into functions under `R/`.
5. Create YAML config files under `config/`.
6. Create `_targets.R` with a minimal working pipeline skeleton first.
7. Add targets module by module, validating each stage.
8. Add audit helpers and module audit outputs.
9. Add plot/table output helpers.
10. Add reports.
11. Add tests.
12. Update README and AGENTS.md.
13. Run lightweight tests and syntax checks.
14. Run `targets::tar_manifest()` and `targets::tar_visnetwork()` or equivalent checks if available.
15. Summarize what was changed, what still needs manual verification, and which targets are expected to be computationally expensive.

---

## Definition of done

This task is complete when:

1. The original Rmd is preserved.
2. The repository has a working `{targets}` structure.
3. Major analysis stages are represented as separate targets.
4. R functions are modular and auditable.
5. Parameters live in YAML config files.
6. Pre-merge and post-merge branches are preserved.
7. Audit reports are generated or scaffolded.
8. Plots and tables are saved to organized output folders.
9. A final report template exists and reads outputs rather than rerunning heavy analysis.
10. README and AGENTS.md explain how to use and maintain the pipeline.
11. The known code issues from the original Rmd are fixed or explicitly documented if not fixable without running data.
12. The codebase is GitHub-ready and avoids committing raw data or large generated outputs.

---

## Final response requested from Claude Code

After implementation, provide a concise summary with:

1. Files created or modified.
2. Pipeline stages implemented.
3. How to run the full pipeline.
4. How to run one module.
5. How to edit parameters.
6. How to find audit outputs.
7. Known limitations or parts requiring manual verification.
8. Any original Rmd code blocks that could not be safely refactored.
