#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(file = args[1])
fullref <- readRDS(file = args[2])
fullref <- subset(x = fullref, cells = Cells(x = ref))
fullref[['umap']] <- ref[['refUMAP']]
Key(object = fullref[['umap']]) <- "umap_"
DefaultAssay(object = fullref[['umap']]) <- "RNA"

DefaultAssay(object = fullref) <- "RNA"
fullref <- NormalizeData(object = fullref)
fullref <- DietSeurat(
  object = fullref,
  dimreducs = "umap",
  assays = "RNA"
)
for (i in colnames(x = fullref[[]])) {
  fullref[[i]] <- NULL
}
fullref <- AddMetaData(object = fullref, metadata = ref[[]])
Misc(object = fullref[['umap']], slot = "model") <- NULL

fullref <- RenameCells(object = fullref, new.names = paste0("cell", 1:ncol(x = fullref)))

for (i in colnames(x = fullref[[]])) {
  if (is.factor(x = fullref[[i, drop = TRUE]])) {
    fullref[[i]] <- as.character(x = fullref[[i, drop = TRUE]])
  }
}

SaveH5Seurat(object = fullref, file = args[3], overwrite = TRUE)
Convert(args[3], dest = "h5ad", overwrite = TRUE)


