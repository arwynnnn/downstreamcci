library(ggplot2)


plotCellChatNetwork <- function(dcc, ...) {
  # Check that the final interactions exist in the DownstreamCCI object.
  if (is.null(dcc$final_interactions)) {
    stop("The DownstreamCCI object does not contain 'final_interactions'. Run Step 3 before (calculateAndFilterInteractions)")
  }
  
  # Extract and clean the final interactions data frame.
  ranked_interactions <- dcc$final_interactions[complete.cases(dcc$final_interactions), ]
  
  # Aggregate scores by source and target cell types.
  celltype_summary <- ranked_interactions %>%
    dplyr::group_by(source_cell_type, target_cell_type) %>%
    dplyr::summarise(
      mean_score = mean(score, na.rm = TRUE),
      n_interactions = dplyr::n(),
      .groups = "drop"
    )
  
  # Create a score matrix with cell types as row and column names.
  cell_types <- unique(c(as.character(celltype_summary$source_cell_type),
                         as.character(celltype_summary$target_cell_type)))
  score_matrix <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
                         dimnames = list(cell_types, cell_types))
  
  for (i in seq_len(nrow(celltype_summary))) {
    src <- as.character(celltype_summary$source_cell_type[i])
    tgt <- as.character(celltype_summary$target_cell_type[i])
    score_matrix[src, tgt] <- celltype_summary$mean_score[i]
  }
  
  # Plot the network using netVisual_circle from CellChat.
  netPlot <- CellChat::netVisual_circle(score_matrix, weight.scale = TRUE, ...)
  return(netPlot)
}


visualize_spatial_feature <- function(seurat_obj, feature, xcoord="X", ycoord="Y", assay = "nn_bulk", pt.size = 1) {
  
  spatial_coords <- seurat_obj@meta.data[, c(xcoord, ycoord)]
  expression_vals <- as.data.frame(t(seurat_obj[[assay]]@data[feature, , drop = FALSE]))
  merged_df <- cbind(spatial_coords, expression = expression_vals[[feature]])
  
  p <- ggplot(merged_df, aes_string(x = xcoord, y = ycoord, color = "expression")) +
    geom_point(size = pt.size) +
    scale_color_viridis_c() +
    theme_minimal() +
    labs(title = paste("Spatial Expression of", feature),
         color = feature) +
    scale_y_reverse()
  
  print(p)
}

visualize_spatial_celltypes <- function(seurat_obj, cell_type_col = "cell_type", xcoord="X", ycoord="Y", celltypes=NULL, pt.size = 1) {
  
  merged_df <- seurat_obj@meta.data
  if (!is.null(celltypes)) {
    merged_df <- merged_df[merged_df[[cell_type_col]] %in% celltypes, ]}
  
  p <- ggplot(merged_df, aes_string(x = xcoord, y = ycoord, color = cell_type_col)) +
    geom_point(size = pt.size) +
    theme_minimal() +
    labs(title = "Spatial Distribution of Cell Types", color = cell_type_col) +
    scale_y_reverse()
  
  print(p)
}


