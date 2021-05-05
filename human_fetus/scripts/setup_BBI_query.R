#!/usr/bin/env Rscript

print("HI")
# read matrix and metadata
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(stringr)
library(magrittr)

ref <- readRDS(file = args[1])
samples <- args[2:15]
df_cell <- readRDS(file = args[16])
df_gene <- readRDS(file = args[17])
rownames(x = df_cell) <- df_cell$sample

objs <- list()
# load and subset
set.seed(seed = 42)
for (i in samples) {
  mat <- readRDS(file = i)
  celltypes <- df_cell[colnames(mat), 'Organ_cell_lineage']
  sample <- c()
  for (j in unique(celltypes)) {
    print(j)
    ct.sample <- which(celltypes == j)
    print(min(length(x = ct.sample), 500))
    sample <- c(sample, sample(ct.sample, size = min(length(x = ct.sample), 1000)))
  }
  mat <- mat[, sample]
  objs <- c(objs, mat)
}

# merge and make object
print(all(sapply(objs, function(o){all(rownames(x = o) == rownames(x = objs[[1]]))})))
merged.obj <- do.call(cbind, objs)
print(all(colnames(x = merged.obj) %in% rownames(x = df_cell)))

df_cell <- df_cell[colnames(x = merged.obj), ]
merged.obj <- merged.obj[, which(x = !is.na(x = df_cell$Organ_cell_lineage))]
df_cell <- df_cell[which(x = !is.na(x = df_cell$Organ_cell_lineage)), ]
merged.obj <- CreateSeuratObject(counts = merged.obj, meta.data = df_cell)

# exclude reference cells and save
merged.obj <- subset(merged.obj, cells = Cells(x = ref), invert = TRUE)

merged.obj <- subset(merged.obj, features = which(x = rowSums(GetAssayData(object = merged.obj[['RNA']], slot = "counts")) > 0))
saveRDS(object = merged.obj, file = args[18])
