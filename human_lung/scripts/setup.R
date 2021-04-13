#!/usr/bin/env Rscript
library(Seurat)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

mat <- fread(file = args[1], header = TRUE, sep = ",")
mat <- as(object = as.matrix(x = mat, rownames = 1), Class = 'dgCMatrix')
meta <- read.csv(file = args[2], row.names = 1)

ob <- CreateSeuratObject(counts = mat, meta.data = meta)
ob <- SCTransform(object = ob)
ob <- RunPCA(object = ob, verbose = FALSE)
ob <- RunUMAP(object = ob, dims = 1:30, return.model = TRUE)

saveRDS(object = ob, file = args[3])
