#!/usr/bin/env Rscript

# Script to build human tonsil reference 2.0.0

# Rscript {input.script} {input.counts} {input.annotations} {input.dr} {output.ref} {output.idx} {output.obj}

library(Seurat)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
obj_original <- readRDS(args[1])
annotations <- read.csv(args[2], row.names = 1)
azimuth.reference <- args[3]
azimuth.index <- args[4]
full.reference <- args[5]

# Exclude multiome cells
obj_original <- subset(obj_original, subset = assay == "3P")

# Rebuild Object from RNA Assay 
obj <- CreateSeuratObject(obj_original[["RNA"]])
obj[["harmony"]] <- obj_original[["harmony"]]

# Add annotations
obj <- AddMetaData(obj, annotations)

# Exclude unknown cells
obj <- subset(obj, subset = celltype.l2 != "unknown")

# Normalize with SCT
obj <- SCTransform(obj, vst.flavor = "v2")

# Build an SNN graph based on their harmony integrated space
obj <- FindNeighbors(obj, reduction = "harmony")

# RunSPCA on the SCT residuals supervised by the harmony SNN graph
obj <- RunSPCA(obj, graph = "RNA_snn", npcs = 50, reduction.name = "spca")

# Build UMAP model from harmony dimension reduction
obj <- RunUMAP(obj, dims=1:50, reduction="harmony", reduction.name = "harmony.umap", return.model = TRUE)

saveRDS(obj, full.reference)

# Build Azimuth Reference
azimuth.obj <- AzimuthReference(
  object = obj,
  refUMAP = "harmony.umap",
  plotref = "harmony.umap",
  refDR = "spca",
  refAssay = "SCT",
  metadata = c("celltype.l1", "celltype.l2"),
  dims = 1:50,
  reference.version = "2.0.0"
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


