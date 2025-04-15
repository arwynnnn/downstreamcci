#!/usr/bin/env Rscript

# -------------------------------
# Script: run_from_checkpoint.R
# Purpose: Run the final AUCell calculations and interaction filtering 
#          starting from the previously saved checkpoint objects.
# -------------------------------

# Load required libraries
library(dplyr)
library(furrr)
library(future)
library(purrr)

# Define the checkpoint file prefix
prefix <- "/home/projects2/kam_project/outputs/checkpoints_r5/ckp_"

# Load checkpoint objects
thr_results         <- readRDS(paste0(prefix, "AUCell_thresholds.rds"))
cells_AUC           <- readRDS(paste0(prefix, "AUCell_raw_AUC.rds"))
auc_matrix          <- readRDS(paste0(prefix, "AUC_matrix.rds"))
interactions_df     <- readRDS(paste0(prefix, "interactions_df_with_ids.rds"))
gene_to_go          <- readRDS(paste0(prefix, "gene_to_go.rds"))
valid_gene_go_terms <- readRDS(paste0(prefix, "valid_gene_go_terms.rds"))
go_term_thresholds  <- readRDS(paste0(prefix, "go_term_thresholds.rds"))

# -------------------------------
# Modified Helper Function 
# (Replace the original to limit exported globals to only what is needed.)
# -------------------------------
compute_median_ratio_and_details <- function(gene, cell, params) {
  valid_terms <- params$valid_gene_go_terms[[gene]]
  if (length(valid_terms) == 0) {
    details_df <- data.frame(
      term = character(0),
      auc_val = numeric(0),
      threshold = numeric(0),
      ratio = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(median = 0, details = details_df))
  }
  auc_vals <- params$auc_matrix[valid_terms, cell, drop = TRUE]
  thr_vals <- params$go_term_thresholds[valid_terms]
  ratios <- auc_vals / thr_vals
  med_val <- median(ratios)
  details_df <- data.frame(
    term = valid_terms,
    auc_val = auc_vals,
    threshold = thr_vals,
    ratio = ratios,
    stringsAsFactors = FALSE
  )
  return(list(median = med_val, details = details_df))
}

# -------------------------------
# Package Heavy Globals
# -------------------------------
future_params <- list(
  valid_gene_go_terms = valid_gene_go_terms,
  go_term_thresholds  = go_term_thresholds,
  auc_matrix          = auc_matrix
)

# -------------------------------
# Set up parallel processing with furrr
# -------------------------------
numCores <- 64  # Adjust number of cores as needed
plan(multisession, workers = numCores)

# -------------------------------
# Compute Results for Ligand (Source)
# -------------------------------
message("Computing results for the ligand (source) for each interaction...")
source_results <- future_map2(
  interactions_df$ligand,
  interactions_df$source_cell,
  function(gene, cell) {
    compute_median_ratio_and_details(gene, cell, future_params)
  }
)

# Extract source medians and details
interactions_df$median_ratio_source <- purrr::map_dbl(source_results, "median")
interactions_df$gene_set_results_source <- purrr::map(source_results, "details")

# -------------------------------
# Compute Results for Receptor (Target)
# -------------------------------
message("Computing results for the receptor (target) for each interaction...")
target_results <- future_map2(
  interactions_df$receptor,
  interactions_df$target_cell,
  function(gene, cell) {
    compute_median_ratio_and_details(gene, cell, future_params)
  }
)

# Extract target medians and details
interactions_df$median_ratio_target <- purrr::map_dbl(target_results, "median")
interactions_df$gene_set_results_target <- purrr::map(target_results, "details")

# -------------------------------
# Compute Composite Score & Filter Interactions
# -------------------------------
message("Computing the composite score normalized by distance...")
interactions_df$composite_score <- (interactions_df$median_ratio_source +
  interactions_df$median_ratio_target) /
  (2 * (interactions_df$distance + 1e-06))

message("Filtering interactions with median ratios > 1...")
final_interactions <- interactions_df %>% 
  filter(median_ratio_source > 1, median_ratio_target > 1) %>% 
  arrange(desc(composite_score))

# If a temporary interaction ID exists, remove it
if ("interaction_id" %in% colnames(final_interactions)) {
  final_interactions <- final_interactions %>% select(-interaction_id)
}

# -------------------------------
# Reset future plan and Save Outputs
# -------------------------------
plan(sequential)

output_dir <- "/home/projects2/kam_project/outputs/"
saveRDS(final_interactions, file = paste0(output_dir, "final_interactions.rds"))
saveRDS(interactions_df, file = paste0(output_dir, "full_interactions.rds"))
message("Pipeline completed successfully from checkpoint.")
