#!/usr/bin/env Rscript

# read in intermediate object
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(stringr)
library(magrittr)
library(Azimuth)

print("Reading object")
bbi <- readRDS(file = args[1])
print("Done reading in")
df_cell <- readRDS(file = args[2])
df_gene <- readRDS(file = args[3])
rownames(x = df_cell) <- df_cell$sample

# annotate based on existing metadata
annotation.l1 <-
  str_match(
    bbi[['Organ_cell_lineage',drop=T]],
    '(?:(?:[A-Za-z]*)-(.*))|(.*)'
  )[,c(2,3)] %>% 
  t %>% 
  as.vector %>% 
  Filter(f = function(a) !is.na(a))
annotation.l2 <- bbi[['Organ_cell_lineage', drop = TRUE]]
organ <- str_match(
    bbi[['Organ_cell_lineage', drop = TRUE]],
    '(?:([A-Za-z]*)-(?:.*))|(.*)'
)[, 2] 
bbi[['annotation.l1']] <- annotation.l1
bbi[['annotation.l2']] <- annotation.l2
bbi[['organ']] <- organ

# subset plotref
set.seed(seed = 42)
plotting.cells <- sample(x = Cells(x = bbi), size = 10**5)
plotref <- subset(x = bbi[["umap"]], cells = plotting.cells)

# make azimuth reference
ref <- AzimuthReference(
    object = bbi,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    plotref = plotref,
    metadata = c('organ','annotation.l1','annotation.l2'), 
    dims = 1:100,
    reference.version = "1.0.0"
)

# save
saveRDS(object = ref, file = args[4])
SaveAnnoyIndex(object = ref[['refdr.annoy.neighbors']], file = args[5])
saveRDS(object = bbi, file = args[6])

