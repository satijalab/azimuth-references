#!/usr/bin/env Rscript

# read matrix and metadata
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)

mat <- readRDS(file = args[1])
df_cell <- readRDS(file = args[2])
df_gene <- readRDS(file = args[3])
rownames(x = df_cell) <- df_cell$sample
rownames(x = mat) <- df_gene[rownames(x = mat), ]$gene_short_name

# prepare object
bbi <- CreateSeuratObject(counts = mat, meta.data = df_cell[colnames(x = mat), ])
Idents(object = bbi) <- 'Organ_cell_lineage'

# normalize
bbi <- SCTransform(
  object = bbi,
  method = 'glmGamPoi', 
  clip.range = c(-10000, 10),
  do.correct.umi = FALSE,
  conserve.memory = TRUE
)

# dimensional reduction
bbi <- RunPCA(object = bbi, npcs = 100)
bbi <- RunUMAP(object = bbi, dims = 1:100, return.model = TRUE)
saveRDS(object = bbi, file = args[4])
