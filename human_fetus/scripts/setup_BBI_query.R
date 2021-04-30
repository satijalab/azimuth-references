#!/usr/bin/env Rscript

# read matrix and metadata
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(stringr)
library(magrittr)
ref <- readRDS(args[1])
samples <- str_split(args[2], ' ')[[1]]
df_cell <- readRDS(args[3])
df_gene <- readRDS(args[4])
rownames(df_cell) <- df_cell$sample

objs <- list()
# load and subset
for (i in samples) {
  mat <- readRDS(i)
  mat <- mat[, sample(1:ncol(mat), size=min(ncol(mat),5000))]
  objs <- c(objs,mat)
}

# merge and make object
print(all(sapply(objs,function(o){all(rownames(o)==rownames(objs[[1]]))})))
merged.obj <- do.call(cbind, objs)
print(all(colnames(merged.obj)%in% rownames(df_cell)))
df_cell <- df_cell[colnames(merged.obj),]
merged.obj <- merged.obj[, which(!is.na(df_cell$Organ_cell_lineage))]
df_cell <- df_cell[which(!is.na(df_cell$Organ_cell_lineage)), ]
merged.obj <- CreateSeuratObject(merged.obj, meta.data = df_cell)

# exclude reference cells and save
merged.obj <- subset(merged.obj, cells=Cells(ref),invert=T)
merged.obj <- subset(merged.obj,features=which(rowSums(merged.obj[['RNA']]@counts)>0))
saveRDS(merged.obj,'out/bbi_query.rds')



