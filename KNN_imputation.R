compute_spatial_neighbors <- function(obj, max_distance = 500, max_neighbors = 32, meta.x = "x", meta.y = "y") {
  spatial_data <- obj@meta.data[, c(meta.x, meta.y)]
  rownames(spatial_data) <- colnames(obj)
  
  # Compute pairwise Euclidean distances
  dist_matrix <- as.matrix(dist(spatial_data))
  nearest_neighbors <- lapply(seq_len(nrow(dist_matrix)), function(i) {
    neighbors <- which(dist_matrix[i, ] <= max_distance & dist_matrix[i, ] > 0) 
    if (length(neighbors) > max_neighbors) {
      neighbors <- neighbors[order(dist_matrix[i, neighbors])][1:max_neighbors]  # Take closest `max_neighbors`
    }
    return(neighbors)
  })
  names(nearest_neighbors) <- colnames(obj)  
  return(nearest_neighbors)
}

pseudobulk_neighbors <- function(
    obj, 
    nn_count = 30, 
    nn_cutoff = 1/2,
    npcs = 15,
    FUN = rowMeans,
    markers = NULL,
    as = c("seurat", "bioc"),
    return_seurat = TRUE,
    normalize_data = TRUE,
    assay = "counts",
    pb = TRUE,
    meta.x = "x",
    meta.y = "y",
    mc.cores = 1) {
  
  # Compute spatial neighbors
  spatial_neighbors <- compute_spatial_neighbors(obj, meta.x=meta.x, meta.y=meta.y)
  
  # Build SNN graph
  if (!is.null(obj@graphs$RNA_snn)) {
    snn_graph <- obj@graphs$RNA_snn
  } else {
    snn_graph <- build_snn(obj, nn_count, npcs)
  }
  
  # Extract expression data
  expr <- extract_expression(obj, assay, normalize_data)
  if (is.null(markers)) markers <- rownames(expr)
  expr <- expr[markers, ]
  cols <- colnames(obj)
  
  num_cells <- ncol(snn_graph)
  message("Pseudo-bulking each cell with its ", nn_count, " neighbors")
  pseudobulked_expr <- mclapply(seq_len(num_cells), function(cell_id) {
    snn_neighbors <- which(snn_graph[, cell_id] > nn_cutoff)
    spatial_neighbors_cell <- spatial_neighbors[[cell_id]]
    
    # Prioritize neighbors present in both SNN and spatial
    common_neighbors <- intersect(snn_neighbors, spatial_neighbors_cell)
    
    if (length(common_neighbors) >= nn_count) {
      neighbors <- common_neighbors[seq_len(nn_count)]
    } else {
      remaining <- setdiff(snn_neighbors, common_neighbors)
      remaining <- remaining[seq_len(min(length(remaining), nn_count - length(common_neighbors)))]
      neighbors <- c(common_neighbors, remaining)
    }
    
    if (length(neighbors) == 1) neighbors <- c(neighbors, neighbors)
    pseudobulked <- FUN(expr[markers, neighbors])
    return(pseudobulked)
  }, mc.cores = mc.cores)
  
  rm(expr)
  pseudobulked_expr <- do.call(cbind, pseudobulked_expr)
  colnames(pseudobulked_expr) <- cols
  
  if (!return_seurat) return(pseudobulked_expr)
  
  message("Adding pseudo-bulked expression data to assay 'nnbulk'")
  pseudobulked_expr <- Seurat::CreateAssayObject(data = pseudobulked_expr, key = "nnbulk_")
  obj[["nnbulk"]] <- pseudobulked_expr
  Seurat::DefaultAssay(obj) <- "nnbulk"
  
  return(obj)
}


extract_expression <- function(obj, assay, normalize_data) { 
  if (class(obj)[1] == "Seurat") { 
    assay <- ifelse(normalize_data, "SCT", "RNA") 
    layer <- ifelse(normalize_data, "data", "counts") 
    expr <- Seurat::GetAssayData(obj, assay = assay, layer = layer) 
  } else if (class(obj)[1] == "SingleCellExperiment") { 
    expr <- SummarizedExperiment::assay(obj, assay = assay) 
  } else { 
    stop("Unsupported object class") 
  } 
  return(expr) 
}

build_snn <- function(obj, nn_count, npcs) { 
  snn_graph <- switch(class(obj)[1], 
                      "Seurat" = build_snn_seurat(obj, nn_count, npcs), 
                      "SingleCellExperiment" = build_snn_sce(obj, nn_count, npcs), 
                      stop("Unsupported object class") 
  ) 
  return(snn_graph) 
}

build_snn_seurat <- function(seu_obj, nn_count, npcs) { 
  if (!"SCT_snn" %in% SeuratObject::Graphs(seu_obj)) { 
    message("Building SNN graph") 
    seu_obj <- Seurat::FindNeighbors( 
      seu_obj,  
      reduction = "pca",  
      dims = npcs, 
      k.param = nn_count, 
      compute.SNN = TRUE) 
  } 
  snn_graph <- SeuratObject::Graphs(seu_obj, slot = "SCT_snn") 
  snn_graph <- as(snn_graph, "dgCMatrix") 
  return(snn_graph) 
}