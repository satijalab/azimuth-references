#!/usr/bin/env Rscript
library(Seurat)
library(feather)
args <- commandArgs(trailingOnly = TRUE)

dat <- readRDS(file = args[1])
metadata <- as.data.frame(x = read_feather(path = args[2]))
rownames(x = metadata) <- metadata$sample_id
metadata <- metadata[colnames(x = dat), ]

ob <- CreateSeuratObject(counts = dat, meta.data = metadata)
ob.list <- SplitObject(object = ob, split.by = "donor_id")
ob.list <- lapply(X = ob.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 3000)
ob.list <- PrepSCTIntegration(object.list = ob.list, anchor.features = features)
ob.list <- lapply(X= ob.list, FUN = RunPCA, verbose = FALSE)

anchors <- FindIntegrationAnchors(
  object.list = ob.list,
  anchor.features = features,
  normalization.method = "SCT",
  reduction = "rpca",
  reference = which(names(ob.list) %in% c("F008", "M008")),
  dims = 1:30
)

ob.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
ob.int <- RunPCA(ob.int, verbose = FALSE)
ob.int <- RunUMAP(ob.int, dims = 1:30, reduction = "pca")
saveRDS(object = ob.int, file = args[3])
