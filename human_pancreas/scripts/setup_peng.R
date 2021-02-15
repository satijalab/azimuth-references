#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)
Convert(args[1], dest = "h5seurat", overwrite = TRUE)
h5Seurat.file <- paste0(
  "data/",
  tools::file_path_sans_ext(x = basename(path = args[1])),
  ".h5seurat"
)

ob <- LoadH5Seurat(file = h5Seurat.file)
ob <- subset(x = ob, subset = Patient == "N7")
for (meta in setdiff(x = colnames(x = ob[[]]), y = 'Cell_type')) {
  ob[[meta]] <- NULL
}
saveRDS(object = ob, file = args[2])
