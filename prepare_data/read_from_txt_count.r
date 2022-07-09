library(GEOquery)
library(Seurat)
library(tibble)
library(magrittr)

#This is the metadata
gse <- getGEO('GSE99933')


# Read in data
data <- data.table::fread("~/Documents/neuroblastoma/rna_seq_data/GSE99933_E12.5_counts.txt", data.table = F)

# Gene names are in the first column so we need to move them to rownames
row.names(data) <- data$V1
#remove the first column as we don't need it anymore
data$V1 <- NULL

# Create Seurat
seurat <- CreateSeuratObject(counts = data)

head(seurat@meta.data, 5)


#subset the seurat object QC
seurat_filtered <- subset(seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#normalize, scale and find variable: steps before pca
seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered, selection.method = "vst", nfeatures = 2000)

seurat_filtered <- RunPCA(seurat_filtered, features = VariableFeatures(object = seurat_filtered))
DimPlot(seurat_filtered, reduction = "pca")
