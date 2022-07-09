#For all datasets in GSE147821
library(dplyr)
library(Seurat)
library(patchwork)

#get data location
dirs <- list.dirs(path= "", recursive=F, full.names = F)
#create seurat objects for each of the files
for (x in dir){
  name <- gsub(".raw_feature_bc_matrix.h5", "", x)
  cts <- Read10X_h5(paste0(dir_path,x))
  assign(name, CreateSeuratObject(counts = cts))
}

#type ls() to see all objects
#merge datasers
merged_seurat <- merge(GSM4446535_10X_19_001, y= ls()[2:3],
                       add.cell.ids = ls()[1:3],
                       project = "NB") #all the other objects
#QC
View(merged_seurat@meta.data)
#create a sample column by copying the rowname of the metadata
merged_seurat$sample <-  rownames(merged_seurat@meta.data)
#split the column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = "sample", into= c("first item", "second item"), sep = "/")

#calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurate, pattern = "^MT-")

#filtering (detailed would be provided in the research paper)
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)