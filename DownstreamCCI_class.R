library(R6)
library(parallel)
library(doSNOW)
library(foreach)

DownstreamCCI <- R6Class("DownstreamCCI",
 public = list(
   # Stored Seurat object and intermediate results.
   seurat_obj = NULL,
   source_target_pass = NULL,
   neighbours = NULL,
   dist_matrix = NULL,
   interactions_annotated = NULL,
   final_interactions = NULL,
   full_interactions = NULL,
   
   # Stored common parameters (set in Step 1).
   assay = NULL,
   cell_type_col = NULL,
   cci_network_all = NULL,
   
   # Constructor: initialize with the Seurat object only.
   initialize = function(seurat_obj) {
     if (is.null(seurat_obj)) {
       stop("A valid Seurat object must be provided.")
     }
     self$seurat_obj <- seurat_obj
   },
   
   # Step 1: Compute Gene Thresholds & Source/Target Pass Vectors.
   # The three specific parameters (assay, cell_type_col, cci_network_all) are passed in here and then stored.
   computeSourceTargetPass = function(assay, cell_type_col, cci_network_all,
                                      high_exp_threshold = 0.5, numCores = -1,
                                      force = FALSE) {
     if (!force && !is.null(self$source_target_pass)) {
       message("Step 1 already computed. Use force = TRUE to re-run.")
       return(invisible(self))
     }
     # Save the common parameters in the object.
     self$assay <- assay
     self$cell_type_col <- cell_type_col
     self$cci_network_all <- cci_network_all
     
     # If numCores is -1, compute available cores.
     if (numCores == -1) {
       numCores <- parallel::detectCores() - 1
     }
     
     res1 <- computeSourceTargetPass(
       seurat_obj = self$seurat_obj,
       cci_network_all = self$cci_network_all,
       assay = self$assay,
       high_exp_threshold = high_exp_threshold,
       cell_type_col = self$cell_type_col,
       numCores = numCores
     )
     self$seurat_obj <- res1$seurat_obj
     self$source_target_pass <- res1$source_target_pass
     invisible(self)
   },
   
   # Step 2: Compute Spatial Neighbours and Annotate Interactions.
   computeNeighboursAndAnnotateInteractions = function(coordinate_cols = c("X", "Y"),
                                                       interaction_distance = 2,
                                                       numCores = -1,
                                                       force = FALSE) {
     # Check that Step 1 has been run.
     if (is.null(self$source_target_pass)) {
       stop("Error: Step 1 (computeSourceTargetPass) must be run before computing neighbours and annotating interactions.")
     }
     
     # If numCores is -1, compute available cores.
     if (numCores == -1) {
       numCores <- parallel::detectCores() - 1
     }
     
     # Call the merged helper function (which integrates your raw code for Steps 2 and 3)
     res2 <- computeNeighboursAndAnnotateInteractions(
       seurat_obj = self$seurat_obj,
       coordinate_cols = coordinate_cols,
       interaction_distance = interaction_distance,
       cci_network_all = self$cci_network_all,
       cell_type_col = self$cell_type_col,
       numCores = numCores
     )
     
     # Update the object with returned values.
     self$seurat_obj <- res2$seurat_obj
     self$neighbours <- res2$neighbours_df
     self$dist_matrix <- res2$dist_matrix
     self$interactions_annotated <- res2$interactions_annotated
     
     invisible(self)
   },
   
   # Step 3: Calculate AUCell Scores, Filter Interactions, & Final Report.
   # This step uses the stored assay.
   calculateAndFilterInteractions = function(aucMaxRank_top_genes = 0.05,
                                             collection = "C2",
                                             numCores = -1,
                                             force = FALSE,
                                             pathway_col="receptor") {
     # Check that the merged Step 2 & 3 has been run.
     if (is.null(self$interactions_annotated)) {
       stop("Error: Step 2 (computeNeighboursAndAnnotateInteractions) must be run before Step 3.")
     }
     if (numCores == -1) {
       numCores <- parallel::detectCores() - 1
     }
     if (!force && !is.null(self$final_interactions)) {
       message("Step 3 already computed. Use force = TRUE to re-run.")
       return(invisible(self))
     }
     finalRes <- calculateAndFilterInteractions(
       seurat_obj = self$seurat_obj,
       interactions_df = self$interactions_annotated,
       collection = collection,
       assay = self$assay,
       aucMaxRank_top_genes = aucMaxRank_top_genes,
       numCores = numCores,
       pathway_col=pathway_col)
     self$final_interactions <- finalRes$final_interactions
     self$full_interactions <- finalRes$full_interactions
     invisible(self)
   }
 )
)
