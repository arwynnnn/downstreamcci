
###############################################################
# Helper Functions for DownstreamCCI Pipeline
# This file contains the helper functions for each step of the pipeline.
# Make sure the required libraries (doSNOW, foreach, etc.) are loaded.
###############################################################

library(doSNOW)   # for parallel loops with progress bars
library(BiocParallel)
library(foreach)
library(msigdbr)
library(dplyr)
library(stringr)
library(AUCell)
library(future)
library(furrr)
library(purrr)

env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)

###############################
# Step 1: Compute Gene Thresholds & Source/Target Pass Vectors
###############################
computeSourceTargetPass <- function(seurat_obj, cci_network_all, assay,
                                    high_exp_threshold, cell_type_col, numCores) {
  message("Step 1: Computing gene thresholds and source/target pass vectors...")
  
  cell_types <- seurat_obj@meta.data[[cell_type_col]]
  unique_cell_types <- unique(cell_types)
  
  # Extract unique genes from the ligand and receptor columns
  unique_genes <- unique(c(cci_network_all$ligand, cci_network_all$receptor))
  # Split composite gene names (separated by ':' or '_')
  all_genes_needed <- unique(unlist(lapply(unique_genes, function(g) {
    if (grepl("[:_]", g)) {
      strsplit(g, "[:_]")[[1]]
    } else {
      g
    }
  })))
  # Retain only genes that are in the assay data
  all_genes_needed <- intersect(all_genes_needed, rownames(seurat_obj[[assay]]@data))
  
  message("  Found ", length(all_genes_needed), " genes present in assay '", assay, "'.")
  
  # Compute thresholds per gene & cell type (using nonzero expression values)
  thresholds <- list()
  message(paste0("  Computing expression thresholds (", high_exp_threshold, " quantile) for each gene..."))
  pb_thresh <- txtProgressBar(min = 0, max = length(all_genes_needed), style = 3)
  for (i in seq_along(all_genes_needed)) {
    g <- all_genes_needed[i]
    expr <- as.numeric(seurat_obj[[assay]]@data[g, ])
    thresholds[[g]] <- sapply(unique_cell_types, function(ct) {
      expr_ct <- expr[cell_types == ct]
      expr_ct_nonzero <- expr_ct[expr_ct != 0]
      if (length(expr_ct_nonzero) > 0) {
        quantile(expr_ct_nonzero, probs = high_exp_threshold, na.rm = TRUE)
      } else {
        NA
      }
    })
    names(thresholds[[g]]) <- unique_cell_types
    setTxtProgressBar(pb_thresh, i)
  }
  close(pb_thresh)
  
  # Helper: compute pass vector for a gene (or composite gene)
  get_pass_vector <- function(gene, seurat_obj, thresholds, cell_types) {
    if (grepl("[:_]", gene)) {
      subunits <- unlist(strsplit(gene, "[:_]"))
      missing <- setdiff(subunits, names(thresholds))
      if (length(missing) > 0) {
        warning(paste("Subunit(s)", paste(missing, collapse = ", "),
                      "of composite gene", gene, "not found; skipping."))
        return(NULL)
      }
      pass_list <- lapply(subunits, function(sub) {
        expr <- as.numeric(seurat_obj[[assay]]@data[sub, ])
        thr <- thresholds[[sub]][cell_types]
        expr >= thr
      })
      pass <- Reduce("&", pass_list)
      return(pass)
    } else {
      if (!(gene %in% names(thresholds))) {
        warning(paste("Gene", gene, "not found in the Seurat object; skipping."))
        return(NULL)
      }
      expr <- as.numeric(seurat_obj[[assay]]@data[gene, ])
      thr <- thresholds[[gene]][cell_types]
      return(expr >= thr)
    }
  }
  
  nCells <- ncol(seurat_obj[[assay]]@data)
  # Setup parallel processing with progress bar using doSNOW
  message("  Determining source and target genes above expression thresholds...")
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  pb_foreach <- txtProgressBar(min = 0, max = length(unique_genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb_foreach, n)
  opts <- list(progress = progress)
  
  res_list <- foreach(gene = unique_genes, .packages = "stats", .options.snow = opts) %dopar% {
    pass <- get_pass_vector(gene, seurat_obj, thresholds, cell_types)
    if (is.null(pass)) {
      return(list(source = rep("", nCells), target = rep("", nCells)))
    }
    source_vec <- rep("", nCells)
    target_vec <- rep("", nCells)
    if (gene %in% cci_network_all$ligand) {
      source_vec[pass] <- gene
    }
    if (gene %in% cci_network_all$receptor) {
      target_vec[pass] <- gene
    }
    list(source = source_vec, target = target_vec)
  }
  close(pb_foreach)
  stopCluster(cl)
  
  # Aggregate results per cell
  source_pass_vec <- sapply(1:nCells, function(i) {
    genes <- sapply(res_list, function(res) res$source[i])
    genes <- genes[genes != ""]
    paste(genes, collapse = ";")
  })
  target_pass_vec <- sapply(1:nCells, function(i) {
    genes <- sapply(res_list, function(res) res$target[i])
    genes <- genes[genes != ""]
    paste(genes, collapse = ";")
  })
  
  # Add new metadata columns to the Seurat object
  seurat_obj[["source_pass"]] <- source_pass_vec
  seurat_obj[["target_pass"]] <- target_pass_vec
  source_target_df <- data.frame(cell = colnames(seurat_obj[[assay]]@data),
                                 source_pass = source_pass_vec,
                                 target_pass = target_pass_vec,
                                 stringsAsFactors = FALSE)
  return(list(seurat_obj = seurat_obj, source_target_pass = source_target_df))
}

###############################
# Step 2: Compute Spatial Neighbours and Annotate Interactions
###############################
computeNeighboursAndAnnotateInteractions <- function(seurat_obj, coordinate_cols, interaction_distance, 
                                                     cci_network_all, cell_type_col, numCores) {
  # ---- Step 2A: Compute Spatial Neighbours ----
  message("Step 2: Computing spatial neighbours for each cell based on coordinates...")
  meta_data <- seurat_obj@meta.data
  if (!all(coordinate_cols %in% colnames(meta_data))) {
    stop("Coordinate columns not found in seurat_obj metadata.")
  }
  coords <- as.matrix(meta_data[, coordinate_cols])
  dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  
  nCells <- nrow(meta_data)
  neighbours_vec <- vector("character", nCells)
  message("  Determining neighbours for each cell within distance ", interaction_distance, "...")
  pb_nb <- txtProgressBar(min = 0, max = nCells, style = 3)
  for (i in 1:nCells) {
    neighbour_indices <- which(dist_matrix[i, ] <= interaction_distance & dist_matrix[i, ] > 0)
    if (length(neighbour_indices) > 0) {
      neighbours_vec[i] <- paste(rownames(meta_data)[neighbour_indices], collapse = ";")
    } else {
      neighbours_vec[i] <- ""
    }
    setTxtProgressBar(pb_nb, i)
  }
  close(pb_nb)
  
  seurat_obj[["neighbours"]] <- neighbours_vec
  neighbours_df <- data.frame(cell = rownames(meta_data),
                              neighbours = neighbours_vec,
                              stringsAsFactors = FALSE)
  
  # ---- Step 2B: Annotate Interactions ----
  meta_data <- seurat_obj@meta.data
  cell_barcodes <- rownames(meta_data)
  
  message("  Finding and annotating interactions between cell...")
  # Set up parallel processing
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  pb_inter <- txtProgressBar(min = 0, max = length(cell_barcodes), style = 3)
  progress <- function(n) setTxtProgressBar(pb_inter, n)
  opts <- list(progress = progress)
  
  edge_results <- foreach(cell = cell_barcodes, .combine = "rbind", 
    .packages = "stats", .options.snow = opts) %dopar% {
      source_cells <- character()
      target_cells <- character()
      interaction_names <- character()
      ligands <- character()
      receptors <- character()
      
      src_genes_str <- meta_data[cell, "source_pass"]
      if (is.na(src_genes_str) || src_genes_str == "") return(NULL)
      src_genes <- unlist(strsplit(src_genes_str, ";"))
      
      nbrs_str <- meta_data[cell, "neighbours"]
      if (is.na(nbrs_str) || nbrs_str == "") return(NULL)
      neighbour_cells <- unlist(strsplit(nbrs_str, ";"))
      
      for (src_gene in src_genes) {
        interactions <- cci_network_all[cci_network_all$ligand == src_gene, ]
        if (nrow(interactions) == 0) next
        for (i in seq_len(nrow(interactions))) {
          target_gene <- interactions$receptor[i]
          for (nbr in neighbour_cells) {
            nbr_target_str <- meta_data[nbr, "target_pass"]
            if (is.na(nbr_target_str) || nbr_target_str == "") next
            nbr_target_genes <- unlist(strsplit(nbr_target_str, ";"))
            if (target_gene %in% nbr_target_genes) {
              source_cells <- c(source_cells, cell)
              target_cells <- c(target_cells, nbr)
              interaction_names <- c(interaction_names, interactions$interaction_name[i])
              ligands <- c(ligands, src_gene)
              receptors <- c(receptors, target_gene)
            }
          }
        }
      }
      
      if (length(source_cells) > 0) {
        data.frame(source_cell = source_cells,
                   target_cell = target_cells,
                   interaction_name = interaction_names,
                   ligand = ligands,
                   receptor = receptors,
                   stringsAsFactors = FALSE)
      } else {
        NULL
      }
    }
  close(pb_inter)
  stopCluster(cl)
  
  if (is.null(edge_results) || nrow(edge_results) == 0) {
    warning("No interactions found.")
    interactions_annotated <- data.frame()
  } else {
    # Merge with pathway info from cci_network_all based on interaction_name.
    mapping_df <- unique(cci_network_all[, c("interaction_name", "pathway_name")])
    inter_res <- merge(edge_results, mapping_df, by = "interaction_name", all.x = TRUE)
    
    # Annotate with distance and cell-type info.
    inter_res$distance <- mapply(function(src, tgt) {
      if (src %in% rownames(dist_matrix) && tgt %in% rownames(dist_matrix)) {
        return(dist_matrix[src, tgt])
      } else {
        return(NA)
      }
    }, inter_res$source_cell, inter_res$target_cell)
    
    inter_res$source_cell_type <- meta_data[inter_res$source_cell, cell_type_col]
    inter_res$target_cell_type <- meta_data[inter_res$target_cell, cell_type_col]
    
    interactions_annotated <- inter_res
  }
  
  return(list(seurat_obj = seurat_obj, 
              neighbours_df = neighbours_df, 
              dist_matrix = dist_matrix, 
              interactions_annotated = interactions_annotated))
}


###############################
# Step 3: Calculate AUCell Scores, Filter Interactions, & Final Report
###############################
calculateAndFilterInteractions <- function(seurat_obj, interactions_df, collection, assay, aucMaxRank_top_genes, numCores, pathway_col) {
  message("Step 3: Calculating enrichment and computing final scores...")
  
  # --- Build Candidate and Gene-to-GO Mappings ---
  # Retrieve GO:BP gene sets from msigdbr.
  msigdb_df <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = "GO:BP")
  
  # Build candidate mapping: each GO term (gs_name) maps to its unique gene set.
  candidate_mapping <- msigdb_df %>%
    group_by(gs_name) %>%
    reframe(genes = list(unique(gene_symbol))) %>%
    ungroup() %>%
    { setNames(.$genes, .$gs_name) }
  
  # Restrict candidate mapping to only those GO terms that include at least one gene from our interactions.
  relevant_genes <- unique(c(interactions_df$ligand, interactions_df$receptor))
  candidate_mapping <- candidate_mapping[sapply(candidate_mapping, function(genes) length(intersect(genes, relevant_genes)) > 0)]
  
  # Build gene-to-GO mapping: for each gene, list all associated GO terms.
  gene_to_go <- msigdb_df %>% 
    group_by(gene_symbol) %>% 
    summarize(go_terms = list(unique(gs_name))) %>% 
    ungroup() %>% 
    { setNames(.$go_terms, .$gene_symbol) }
  
  # --- Use ALL cells for AUCell ranking (global context) ---
  expr_dense_full <- as.matrix(seurat_obj[[assay]]@data)
  expr_dense <- expr_dense_full  # Using all cells ensures more robust rankings
  
  message("  Building AUCell rankings for assay '", assay, "'...")
  cells_rankings <- AUCell_buildRankings(expr_dense, plotStats = FALSE, verbose = TRUE)
  
  aucMaxRank_value <- ceiling(aucMaxRank_top_genes * nrow(expr_dense))
  message(paste0("  Calculating AUC scores using top ", aucMaxRank_top_genes * 100, "% genes..."))
  if(numCores == -1) { 
    numCores <- parallel::detectCores() - 1
  }
  cells_AUC <- AUCell_calcAUC(candidate_mapping, cells_rankings, aucMaxRank = aucMaxRank_value, nCores = numCores, verbose = TRUE)
  
  message("  Computing thresholds for enrichment scores...")
  thr_results <- AUCell_exploreThresholds(cells_AUC, plotHist = FALSE, assign = TRUE, nCores = numCores, verbose = TRUE)
  auc_matrix <- getAUC(cells_AUC)

  
  message("  Computing composite scores and final matrix...")
  # Add interaction_id temporarily to join with AUC dat
  
  # For Linux, use multicore.
  plan(multisession, workers=numCores-2)
  
  # Add an interaction ID.
  interactions_df <- interactions_df %>% mutate(interaction_id = row_number())
  
  # Precompute valid GO terms for each gene:
  # Only GO terms in auc_matrix are retained.
  message("   Precomputing valid GO terms for each gene...")
  valid_gene_go_terms <- lapply(gene_to_go, function(go_terms) {
    if (is.null(go_terms) || length(go_terms) == 0) 
      return(character(0))
    intersect(go_terms, rownames(auc_matrix))
  })
  
  # Precompute thresholds for each GO term in the auc_matrix.
  message("   Precomputing thresholds for each GO term...")
  go_term_thresholds <- sapply(rownames(auc_matrix), function(term) {
    if (!is.null(thr_results[[term]]) && !is.null(thr_results[[term]]$aucThr$selected)) {
      thr_results[[term]]$aucThr$selected
    } else {
      1
    }
  }, USE.NAMES = TRUE)

  # CHECKPOINT
  prefix <- "/home/projects2/kam_project/outputs/checkpoints_r5/ckp_"
  saveRDS(thr_results, paste0(prefix, "AUCell_thresholds.rds"))
  saveRDS(cells_AUC, paste0(prefix, "AUCell_raw_AUC.rds"))
  saveRDS(auc_matrix, paste0(prefix, "AUC_matrix.rds"))
  saveRDS(interactions_df, paste0(prefix, "interactions_df_with_ids.rds"))
  saveRDS(gene_to_go, paste0(prefix, "gene_to_go.rds"))
  saveRDS(valid_gene_go_terms, paste0(prefix, "valid_gene_go_terms.rds"))
  saveRDS(go_term_thresholds, paste0(prefix, "go_term_thresholds.rds"))

  future_params <- list(
  valid_gene_go_terms = valid_gene_go_terms,
  go_term_thresholds = go_term_thresholds,
  auc_matrix = auc_matrix)
  
  # Define a function that computes both the median ratio and details for a given gene in a given cell.
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
  
  # Compute results for the ligand (source) for each interaction in parallel.
  message("   Computing results for the ligand (source) for each interaction...")
  source_results <- future_map2(
  interactions_df$ligand,
  interactions_df$source_cell,
  function(gene, cell) {
    compute_median_ratio_and_details(gene, cell, future_params)
    }
  )
  
  # Extract the median ratio and details for source.
  message("   Extracting the median ratio and details for source...")
  interactions_df$median_ratio_source <- map_dbl(source_results, "median")
  interactions_df$gene_set_results_source <- map(source_results, "details")
  
  # Compute results for the receptor (target) for each interaction in parallel.
  message("   Computing results for the receptor (targer) for each interaction...")
  target_results <- future_map2(
    interactions_df$receptor,
    interactions_df$target_cell,
    function(gene, cell) {
      compute_median_ratio_and_details(gene, cell, future_params)
    }
  )
  
  # Extract the median ratio and details for target.
  message("   Extracting the median ratio and details for target...")
  interactions_df$median_ratio_target <- map_dbl(target_results, "median")
  interactions_df$gene_set_results_target <- map(target_results, "details")
  
  # Compute the composite score normalized by distance.
  message("   Computing the composite score normalized by distance...")
  interactions_df$composite_score <- (interactions_df$median_ratio_source +
    interactions_df$median_ratio_target) / (2 * (interactions_df$distance + 1e-06))
  
  # Filter interactions to keep only those with both median ratios exceeding 1 and sort by composite score.
  message("   Filtering interactions...")
  ranked_interactions <- interactions_df %>% 
    filter(median_ratio_source > 1, median_ratio_target > 1) %>% 
    arrange(desc(composite_score))
  
  # Remove the temporary interaction_id column before returning.
  message("   Cleaning the image...")
  final_interactions <- ranked_interactions %>% select(-interaction_id)
  interaction_results <- interactions_df %>% select(-interaction_id)
  
  # Return the final objects.
  message("Pipeline completed.")
  list(
    final_interactions = final_interactions,
    full_interactions = interaction_results,
    auc_matrix = auc_matrix
  )
}
