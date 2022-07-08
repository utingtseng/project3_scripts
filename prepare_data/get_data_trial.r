library(GEOquery)
library(Seurat)
library(tibble)
library(magrittr)

#This is the metadata
gse <- getGEO('GSE99933')


# Read in data
data <- data.table::fread("~/Documents/neuroblastoma/rna_seq_data/GSE99933_E12.5_counts.txt", data.table = F)

# Gene names are in the first column so we need to move them to rownames
data <- data %>% 
  column_to_rownames("V1")

# Create Seurat
seurat <- CreateSeuratObject(counts = data)

head(seurat@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#normalize the data
seurat <- NormalizeData(seurat)

seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
DimPlot(seurat, reduction = "pca")
