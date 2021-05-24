#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(file = args[1])
fullref <- LoadH5Seurat(file = args[2], assays = "ADT", reductions = FALSE)
fullref <- subset(x = fullref, cells = Cells(x = ref))
dat_3p <- ReadMtx(mtx = args[3], cells = args[4], features = args[5])

fullref[["RNA"]] <- CreateAssayObject(counts = dat_3p[, Cells(x = ref)])
DefaultAssay(object = fullref) <- "RNA"
fullref <- NormalizeData(object = fullref)
adt.mat <- GetAssayData(object = fullref[["ADT"]], slot = "data")
rownames(x = adt.mat) <- paste0("ADT-", rownames(x = adt.mat))

combined.rna.adt <- CreateAssayObject(data = rbind(
    GetAssayData(object = fullref[["RNA"]], slot = "data"),
    adt.mat)
)
fullref[["RNA"]] <- combined.rna.adt

fullref[['umap']] <- ref[['refUMAP']]
Key(object = fullref[['umap']]) <- "umap_"
DefaultAssay(object = fullref[['umap']]) <- "RNA"
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

SaveH5Seurat(object = fullref, file = args[6], overwrite = TRUE)
Convert(args[6], dest = "h5ad", overwrite = TRUE)


