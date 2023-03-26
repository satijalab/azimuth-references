#!/usr/bin/env Rscript

# Script to build human liver demo

library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
demo.path <- args[1] 
demo.output <- args[2]

demo <- Read10X(demo.path, gene.column = 1)
demo.obj <- CreateSeuratObject(counts = demo, min.cells = 3, min.features = 200)
demo.obj <- subset(demo.obj, subset = nFeature_RNA < 10000 & nCount_RNA < 60000)
demo.obj <- subset(demo.obj, downsample = 50000)
saveRDS(demo.obj, demo.output)
