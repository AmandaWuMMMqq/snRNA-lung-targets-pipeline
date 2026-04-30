#!/usr/bin/env Rscript
# Run a selected module or target.
# Usage: Rscript scripts/run_module.R <module_name>
#
# Available module names:
#   config, input, individual_integration, annotation,
#   joint_integration, epithelial, alveolar, trajectory,
#   tradeseq, cellchat, composition, final_plots, final_tables,
#   audit_report, final_report, all

library(targets)

module_map <- list(
  config               = c("paths", "params", "markers", "annotations", "plotting"),
  input                = c("fetal_raw", "pediatric_raw", "input_audit"),
  individual_integration = c("fetal_integrated", "pediatric_integrated",
                              "individual_integration_audit"),
  annotation           = c("fetal_annotated", "pediatric_annotated", "annotation_audit"),
  joint_integration    = c("joint_integrated", "joint_annotated", "joint_integration_audit"),
  epithelial           = c("pre_epi", "post_epi", "epithelial_branch_audit"),
  alveolar             = c("pre_alveolar", "post_alveolar", "alveolar_audit"),
  trajectory           = c("pre_slingshot", "post_slingshot", "trajectory_audit"),
  tradeseq             = c("pre_tradeseq", "post_tradeseq",
                            "pre_tradeseq_results", "post_tradeseq_results",
                            "tradeseq_audit"),
  cellchat             = c("pre_cellchat", "post_cellchat", "cellchat_audit"),
  composition          = c("composition_tables", "composition_plots", "composition_audit"),
  final_plots          = "final_plots",
  final_tables         = "final_tables",
  audit_report         = "module_audit_report",
  final_report         = "final_report"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] == "all") {
  cat("Running full pipeline...\n")
  tar_make()
} else {
  module <- args[1]
  if (!module %in% names(module_map)) {
    cat("Unknown module:", module, "\n")
    cat("Available modules:", paste(names(module_map), collapse = ", "), "\n")
    quit(status = 1)
  }
  targets_to_run <- module_map[[module]]
  cat("Running module:", module, "->", paste(targets_to_run, collapse = ", "), "\n")
  tar_make(names = all_of(targets_to_run))
}
