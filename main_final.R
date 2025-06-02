
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
seurat_obj_min <- subset_seurat(seurat_obj, xmin=1000, xmax=1500, ymin=1250, ymax=2000)

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
  high_exp_threshold = 0.25,
  numCores = 16)

# Step 2: Compute Spatial Neighbours and Annotate Interactions.
dcc$computeNeighboursAndAnnotateInteractions(coordinate_cols = c("X", "Y"), interaction_distance = 100, numCores = 16)

# Step 3: Calculate AUCell Scores, Filter Interactions
dcc$calculateAndFilterInteractions(aucMaxRank_top_genes = 0.20, collection="C5", pathway_col="receptor", numCores = 16)

# save the file
saveRDS(dcc, file = "/home/projects2/kam_project/outputs2/dcc_full.rds")
