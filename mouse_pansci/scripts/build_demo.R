#!/usr/bin/env Rscript

# Script to build demo

library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
demo.path <- args[1] 
demo.output <- args[2]

demo <- readRDS(demo.path, gene.column = 1)
demo <- subset(demo, cells = sample(Cells(demo), size = 20000))

saveRDS(demo.obj, demo.output)
