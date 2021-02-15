#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)
args <- commandArgs(trailingOnly = TRUE)

ref.dir <- "reference/"
ob.dir <- "seurat_objects/"
ref <- readRDS(file = args[1])

Idents(object = ref) <- "subclass_label"
ref$class <- ref$class_label
ref$cluster <- ref$cluster_label
ref$subclass <- ref$subclass_label
ref$cross_species_cluster <- ref$cross_species_cluster_label

full.ref <- ref
ref <- subset(x = ref, downsample = 2000)
ref <- RunUMAP(object = ref, reduction = "pca", dims = 1:30, return.model = TRUE)

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "integrated",
  metadata = c("class", "cluster", "subclass", "cross_species_cluster"),
  dims = 1:50,
  k.param = 31,
  reference.version = "1.0.0"
)
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = full.ref, file = file.path(ob.dir, "fullref.Rds"))
