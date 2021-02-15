#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

panc8 <- readRDS(file = args[1])
xin <- readRDS(file = args[2])
ob.list <- c(panc8, xin)
features <- SelectIntegrationFeatures(object.list = ob.list)
ob.list <- PrepSCTIntegration(object.list = ob.list, anchor.features = features)
ob.list <- lapply(X = ob.list, FUN = RunPCA, features = features, verbose = FALSE)
ref <- unname(obj = which(x = sapply(X = ob.list, FUN = function(x) x$dataset[1] == "indrop1")))
anchors <- FindIntegrationAnchors(
  object.list = ob.list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  reference = ref,
  dims = 1:30
)
panc8_xin <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  normalization.method = "SCT",
)
panc8_xin <- RunPCA(object = panc8_xin)
panc8_xin <- RunUMAP(object = panc8_xin, dims = 1:30, return.model = TRUE)
panc8_xin <- FindNeighbors(object = panc8_xin, dims = 1:30)
panc8_xin <- FindClusters(object = panc8_xin, resolution = c(0.8, 2))
saveRDS(object = panc8_xin, file = args[3])
