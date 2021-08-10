#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(Matrix)

mat <- readMM(file = "data/kidney_demo_expression.mtx")
cells <- read.csv(file = "data/kidney_demo_cells.csv", header = FALSE)
features <- read.csv(file = "data/kidney_demo_features.csv", header = FALSE)
meta <- read.csv(file = "data/kidney_demo_metadata.csv", row.names = 1)
meta <- meta[, c("Project", "Experiment", "celltype", "compartment", "broad_celltype")]
cells.keep <- read.csv("data/stewart_cells.csv", row.names = 1)
rownames(x = mat) <- cells$V1
colnames(mat) <- features$V1
mat <- mat[cells.keep[, 1], ]
mat <- t(x = mat)

ob <- CreateSeuratObject(counts = mat, meta.data = meta)
saveRDS(object = ob, file = args[1])