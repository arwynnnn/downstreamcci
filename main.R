#setwd("~/#DTU/KU/nichenet")
source("KNN_imputation.R")
source("DownstreamCCI.R")
source("DownstreamCCI_utils.R")
library(Seurat)
library(SeuratObject)

### LOAD DATA ###

# expression
HumanTonsil_spatial <- read.csv("/home/projects2/kam_project/data/HumanTonsil_spatial.csv")
HumanTonsil_expression = as.data.frame(fread("/home/projects2/kam_project/data/HumanTonsil_expression.csv.gz"))
rownames(HumanTonsil_expression) <- HumanTonsil_expression$V1
HumanTonsil_expression <- HumanTonsil_expression[, -1]

# cell-cell interaction prior knowledge
CellChatDB <- CellChatDB.human
cci_network_all <- CellChatDB$interaction[CellChatDB$interaction$annotation == "Cell-Cell Contact", ]


### PREPROCESSING ###

# load into seurat
seurat_obj <- CreateSeuratObject(counts = HumanTonsil_expression)
rownames(HumanTonsil_spatial) <- HumanTonsil_spatial$NAME
seurat_obj <- AddMetaData(seurat_obj, metadata = HumanTonsil_spatial)

#view the sample and filter based on cell types
visualize_spatial_celltypes(seurat_obj, pt.size=2.5, celltypes=c("B_germinal_center", "B_naive"))
#seurat_obj_min <- subset_seurat(seurat_obj, xmin=750, xmax=1500, ymin=1100, ymax=2250, celltypes=c("B_germinal_center", "B_naive"))
seurat_obj_min <- subset_seurat(seurat_obj, xmin=750, xmax=1000, ymin=1000, ymax=1250)
visualize_spatial_celltypes(seurat_obj_min, pt.size=2.5)

# data normalisation and imputation
seurat_obj_min <- SCTransform(seurat_obj_min)
seurat_obj_min <- RunPCA(seurat_obj_min, features = VariableFeatures(seurat_obj_min))
seurat_obj_min <- pseudobulk_neighbors(seurat_obj_min, meta.x = "X", meta.y = "Y")


### RUN CCI ANALYSIS ###

# Create a DownstreamCCI container with the Seurat object.
dcc <- DownstreamCCI$new(seurat_obj_min)

# Step 1: Compute Gene Thresholds & Source/Target Pass Vectors.
dcc$computeSourceTargetPass(
  assay = "nnbulk",
  cell_type_col = "cell_type",
  cci_network_all = cci_network_all,
  high_exp_threshold = 0.8,
  numCores = -1)

# Step 2: Compute Spatial Neighbours and Annotate Interactions.
dcc$computeNeighboursAndAnnotateInteractions(coordinate_cols = c("X", "Y"), interaction_distance = 5, numCores = -1)

# Step 3: Calculate AUCell Scores, Filter Interactions
dcc$calculateAndFilterInteractions(aucMaxRank_top_genes = 0.05, collection="C2", pathway_col="receptor", numCores = 1)

# save the file
saveRDS(dcc, file = "/home/projects2/kam_project/downstreamcci/outputs/dcc_full_sample.rds")

# enrichment vs distance
#printInteractionNumbers(dcc)
#plotEnrichmentVsDistance(dcc, "CD99_CD99", enrichment_metric = "ratio")

# Visualize the network using the final interactions stored in the object.
#plotCellChatNetwork(
#  dcc=dcc,
#  title.name = "Cell-Cell Communication Network",
#  arrow.width = 2,
#  arrow.size = 0.6)

