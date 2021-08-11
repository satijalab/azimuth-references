#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(Azimuth)

ob <- readRDS(file = args[1])

# small naming changes requested by investigators
old.names <- c("Monocyte","Tissue-resident Macrophage","Natural Killer / Natural Killer T")
new.names <- c("Non-classical monocyte", "M2 Macrophage", "Natural Killer T")

l3 <- as.character(ob$annotation.l3)
l2 <- as.character(ob$annotation.l2)

for(i in 1:3) {
  l3[l3==old.names[i]] <- new.names[i]
  l2[l2==old.names[i]] <- new.names[i]
}

ob$annotation.l2 <- factor(l2,levels = sort(unique(l2)))
ob$annotation.l3 <- factor(l3,levels = sort(unique(l3)))

ref <- AzimuthReference(
  object = ob, 
  refDR = 'pca',
  metadata = c('annotation.l3', 'annotation.l2', 'annotation.l1'),
  dims = 1:100,
  reference.version = "1.0.0"
)

saveRDS(object = ref, file = args[2])
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = args[3])

