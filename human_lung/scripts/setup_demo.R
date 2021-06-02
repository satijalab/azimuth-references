#!/usr/bin/env Rscript
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

load(file = args[1])
meta <- read.table(file = args[2], header = T, sep = "\t", row.names = 1)
meta <- meta[, c("ID", "location", "celltype")]

ob <- CreateSeuratObject(counts = raw_counts, meta.data = meta)
saveRDS(object = ob, file = args[3])
