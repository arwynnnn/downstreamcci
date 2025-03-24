# Load required libraries
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(data.table)
library(foreach)
library(parallel)
library(doSNOW)   # use doSNOW to allow parallel loops with progress bars
library(CellChat)
library(AUCell)
source("DownstreamCCI_pipeline.R")
source("DownstreamCCI_class.R")
source("DownstreamCCI_visualizations.R")

# Main function that calls the helper functions
runDownstreamCCI <- function(seurat_obj = NULL,
                             cci_network_all = NULL,
                             cellchat_db = NULL,
                             downstreamcci_object = NULL,
                             assay = "RNA",
                             high_exp_threshold = 0.5,
                             interaction_distance = 2,
                             aucMaxRank_top_genes = 0.05,
                             cell_type_col = "cell_type",
                             coordinate_cols = c("X", "Y"),
                             numCores = -1,
                             startFrom = 1) {   # startFrom: 1 to 4 (step number) to re-run from that step
  
  if (is.null(downstreamcci_object)) {
    DownstreamCCI <- list()
  } else {
    DownstreamCCI <- downstreamcci_object
  }
  
  # Print the statement on the amount of cores used
  if (numCores == -1) {
    numCores <- parallel::detectCores() - 1
  }
  message("Using ", numCores, " core(s) for parallel processing.")
  
  # Check for cell type info; if missing, assign default.
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    message("Cell type column '", cell_type_col, "' not found. Assigning all cells the same type ('All').")
    seurat_obj@meta.data[[cell_type_col]] <- "All"
  }
  
  ###########################################################################
  # Step 1: Compute Gene Thresholds & Source/Target Pass Vectors
  ###########################################################################
  if (startFrom <= 1 || !("source_target_pass" %in% names(DownstreamCCI))) {
    res1 <- computeSourceTargetPass(seurat_obj, cci_network_all, assay,
                                    high_exp_threshold, cell_type_col, numCores)
    seurat_obj <- res1$seurat_obj
    DownstreamCCI$source_target_pass <- res1$source_target_pass
  } else {
    message("Step 1: Using existing source/target pass output.")
  }
  
  ###########################################################################
  # Step 2: Compute Spatial Neighbours
  ###########################################################################
  if (startFrom <= 2 || !("neighbours" %in% colnames(seurat_obj@meta.data))) {
    res2 <- computeNeighbours(seurat_obj, coordinate_cols, interaction_distance)
    seurat_obj <- res2$seurat_obj
    DownstreamCCI$neighbours <- res2$neighbours_df
    dist_matrix <- res2$dist_matrix
  } else {
    message("Step 2: Using existing neighbours output.")
    dist_matrix <- as.matrix(dist(as.matrix(seurat_obj@meta.data[, coordinate_cols]), method = "euclidean"))
  }
  
  ###########################################################################
  # Step 3: Find and Annotate Interactions
  ###########################################################################
  if (startFrom <= 3 || !("interactions_annotated" %in% names(DownstreamCCI))) {
    inter_annotated <- findAndAnnotateInteractions(seurat_obj, cci_network_all, cell_type_col, numCores, dist_matrix)
    DownstreamCCI$interactions_annotated <- inter_annotated
  } else {
    message("Step 3: Using existing interactions annotation.")
  }
  
  ###########################################################################
  # Step 4: Calculate AUCell Scores, Filter Interactions, & Final Report
  ###########################################################################
  if (!is.null(cellchat_db)) {
    if (startFrom <= 4 || !("final_interactions" %in% names(DownstreamCCI))) {
      finalRes <- calculateAndFilterInteractions(seurat_obj, DownstreamCCI$interactions_annotated,
                                                 cellchat_db, assay, aucMaxRank_top_genes, numCores)
      final_int <- finalRes$final_interactions
      # Remove rows with NA values to ensure only complete interactions are kept.
      final_int <- final_int[complete.cases(final_int), ]
      DownstreamCCI$final_interactions <- final_int
    } else {
      message("Step 4: Using existing final interactions output.")
    }
  } else {
    message("CellChat database not provided; skipping AUCell-based filtering. Final interactions set to annotated interactions.")
    DownstreamCCI$final_interactions <- DownstreamCCI$interactions_annotated
  }
  
  message("DownstreamCCI pipeline completed successfully.")
  return(list(seurat_obj = seurat_obj, DownstreamCCI = DownstreamCCI))
}
