#getting familiar with r, following tutorial on youtube
library(GEOquery)
library(dplyr)
library(tidyverse)
library(tidyr)

data_path <- ("~/Documents/neuroblastoma/rna_seq_data")
#load expression data
dat <- read.csv(file.path(data_path, "GSE99933_E12.5_counts.txt"), sep='\t')
#get dimension of the dataframe
dim(dat)

#get metadata
gse<- getGEO(GEO="GSE99933", GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#select columns of interest with a pipeline
metadata.modified <- metadata %>%
  select(1,10,11)%>%
  #rename the column name
  rename(celltype = characteristics_ch1.1)%>%
  #change the content of the string
  #gsub globally substitute cell type to nothing
  mutate(celltype = gsub("cell type:", "", celltype))%>%
  rename(genotype= characteristics_ch1)%>%
  mutate(genotype = gsub("genotype", "", genotype))%>%
  #\\ is used to escape the []
  mutate(title = gsub(" \\[single cell RNA-seq\\]","", title))

#look at the data, difficult to add metadata to the wide format (samples as columns)
head(dat)

#adding a column named genes (essentially a copy from rownames)
dat <- data.frame(genes = row.names(dat), dat)
rownames(dat)<- NULL

#reshaping the data
#somehting is wrong here
data.long <- dat%>%
  #-gene means I exclude the gene column
  gather("sample", "FPKM", -genes)%>%
  head()

#join dataframes
data.long <- data.long%>%
  left_join(., metadata.modified, by = c("sample"= "title"))%>%
  head()

#explore data
data.long %>%
  filter(genes == "Mrpl15")%>%
  group_by(genes, sample)%>%
  summarize(mean_FPKM = mean(FPKM))%>%
  head()

#something is off??
