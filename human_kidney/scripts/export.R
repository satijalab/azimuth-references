#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(Azimuth)

ob <- readRDS(file = args[1])

ref <- AzimuthReference(
  object = ob, 
  refDR = 'pca',
  metadata = c('annotation.l3', 'annotation.l2', 'annotation.l1', 'celltype'),
  dims = 1:100,
  reference.version = "1.0.0"
)

saveRDS(object = ref, file = args[2])
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = args[3])

