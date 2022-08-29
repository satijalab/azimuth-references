#!/usr/bin/env Rscript

# Script to build human tonsil reference 1.0.0

# Rscript {input.script} {input.counts} {input.annotations} {input.dr} {output.ref} {output.idx} {output.obj}

library(Seurat)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
obj <- readRDS(args[1])
annotations <- read.csv(args[2])
rownames(annotations) <- annotations$X
annotations$X <- NULL
azimuth.reference <- args[3]
azimuth.index <- args[4]
full.reference <- args[5]

# Add annotations
obj <- AddMetaData(obj, annotations)

# Rename cell types
obj$celltype.l1 <- obj$annotation_figure_1_renamed
obj$celltype.l2 <- obj$annotation_20220215
obj$celltype.l1 <- as.factor(obj$celltype.l1)
obj$celltype.l2 <- as.factor(obj$celltype.l2)

# Normalize with SCT
obj <- SCTransform(obj, vst.flavor = "v2")

# Build an SNN graph based on their harmony integrated space
obj <- FindNeighbors(obj, reduction = "harmony")

# RunSPCA on the SCT residuals supervised by the harmony SNN graph
obj <- RunSPCA(obj, graph = "RNA_snn", npcs = 50, reduction.name = "spca")

# Build UMAP model from harmony dimension reduction
obj <- RunUMAP(obj, dims=1:50, reduction="harmony", reduction.name = "harmony.umap", return.model = TRUE)

saveRDS(obj, full.reference)

azimuth.obj <- AzimuthReference(
  object = obj,
  refUMAP = "harmony.umap",
  plotref = "harmony.umap",
  refDR = "spca",
  refAssay = "SCT",
  metadata = c("celltype.l1", "celltype.l2"),
  dims = 1:50,
  reference.version = "1.0.0"
)
ValidateAzimuthReference(object = azimuth.obj)
SaveAnnoyIndex(
  object = azimuth.obj[["refdr.annoy.neighbors"]],
  file = azimuth.index
)
saveRDS(
  object = azimuth.obj,
  file = azimuth.reference
)
