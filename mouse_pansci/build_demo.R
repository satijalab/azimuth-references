#!/usr/bin/env Rscript

# Script to build human liver demo

# Demo object can be downloaded at https://cellxgene.cziscience.com/collections/fe0e718d-2ee9-42cc-894b-0b490f437dfd

library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
demo.path <- args[1] 
demo.output <- args[2]

demo <- readRDS(demo.path, gene.column = 1)
demo <- subset(demo, cells = sample(Cells(demo), size = 20000))

saveRDS(demo.obj, demo.output)
