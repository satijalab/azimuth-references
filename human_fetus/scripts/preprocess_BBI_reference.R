#!/usr/bin/env Rscript

# read matrix and metadata
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(stringr)
library(magrittr)
mat <- readRDS(args[1])
df_cell <- readRDS(args[2])
df_gene <- readRDS(args[3])
rownames(df_cell) <- df_cell$sample

# prepare object
bbi <- CreateSeuratObject(mat, meta.data = df_cell[colnames(mat), ])
Idents(bbi) <- 'Organ_cell_lineage'

# normalize
bbi <- SCTransform(
  bbi,
  method = 'glmGamPoi', 
  clip.range = c(-10000, 10),
  do.correct.umi = F,
  conserve.memory = T
)

# dimensional reduction
bbi <- RunPCA(bbi, npcs=100)
bbi <- RunUMAP(bbi, dims=1:100, return.model=T)
saveRDS(bbi,'out/ref_intermediate.Rds',compress=F)
