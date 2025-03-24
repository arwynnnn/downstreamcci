subset_seurat <- function(seurat_obj, xmin, xmax, ymin, ymax, celltypes = NULL, cell_type_col = "cell_type") {

  md <- seurat_obj@meta.data
  
  # Subset based on X and Y coordinate ranges
  cells_to_keep <- rownames(md)[md$X >= xmin & md$X <= xmax & md$Y >= ymin & md$Y <= ymax]
  
  # Further subset by cell types if specified
  if (!is.null(celltypes)) {
    cells_to_keep <- intersect(cells_to_keep, rownames(md)[md[[cell_type_col]] %in% celltypes])
  }
  subset(seurat_obj, cells = cells_to_keep)
}


printInteractionNumbers <- function(dcc, top_n = 10) {
  library(dplyr)
  
  counts <- dcc$full_interactions %>%
    group_by(interaction_name) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  print(head(counts, top_n))
}


plotEnrichmentVsDistance <- function(dcc, selected_interaction, enrichment_metric = "percentage") {
  library(dplyr)
  library(ggplot2)
  
  # Filter full_interactions for only the selected interaction
  sub_int <- dcc$full_interactions %>%
    filter(interaction_name == selected_interaction)
  
  if(nrow(sub_int) == 0) {
    stop("No interactions found for the specified interaction name.")
  }
  
  # Aggregate enrichment metric per target cell
  cell_enrichment <- if(enrichment_metric == "percentage") {
    sub_int %>%
      group_by(target_cell) %>%
      summarise(enrichment = mean(enriched_percentage, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(y_label = "% Enriched Gene Sets")
  } else if(enrichment_metric == "ratio") {
    sub_int %>%
      group_by(target_cell) %>%
      summarise(enrichment = mean(median_ratio, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(y_label = "Median Enrichment Ratio")
  } else {
    stop("enrichment_metric must be either 'percentage' or 'ratio'.")
  }
  
  # Ensure the target cells are available in the distance matrix
  available_cells <- intersect(cell_enrichment$target_cell, rownames(dcc$dist_matrix))
  if(length(available_cells) == 0) {
    stop("No target cells found in the distance matrix.")
  }
  
  # Compute, for each available cell, the distance to the closest other cell (with the same interaction)
  min_distances <- sapply(available_cells, function(cell) {
    dists <- dcc$dist_matrix[cell, available_cells]
    dists <- dists[dists > 0]  # remove self-distance (0)
    if(length(dists) == 0) return(NA) else return(min(dists, na.rm = TRUE))
  })
  
  dist_df <- data.frame(target_cell = available_cells,
                        min_distance = min_distances,
                        stringsAsFactors = FALSE)
  
  # Merge enrichment and distance data
  plot_df <- merge(cell_enrichment, dist_df, by = "target_cell")
  
  # Plot enrichment vs. minimum distance
  p <- ggplot(plot_df, aes(x = min_distance, y = enrichment)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Distance to Closest Interacting Partner", 
         y = unique(plot_df$y_label),
         title = paste("Enrichment vs Distance for", selected_interaction))
  
  print(p)
}
