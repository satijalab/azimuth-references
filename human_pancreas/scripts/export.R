#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

ref.dir <- "reference/"
ob.dir <- "seurat_objects/"
ref <- readRDS(file = args[1])
annotations <- readRDS(file = args[2])
Idents(object = ref) <- annotations

if ("remove" %in% levels(x = ref)) {
  ref <- subset(x = ref, idents = "remove", invert = TRUE)
  ref <- RunPCA(object = ref, verbose = FALSE)
}
ref$annotation.l1 <- Idents(object = ref)
ref <- RunUMAP(object = ref, dims = 1:30, return.model = TRUE)
full.ref <- ref

colormap <- list(annotation.l1 = CreateColorMap(object = ref, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "integrated",
  metadata = c("annotation.l1"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = full.ref, file = file.path(ob.dir, "fullref.Rds"))

