#!/usr/bin/env Rscript

# Script to build human heart demo

library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
demo.path <- args[1]
demo.output <- args[2]

print("got args")
Convert(demo.path, dest = "h5seurat", overwrite = TRUE)
print("converted")
demo <- LoadH5Seurat(gsub("h5ad", "h5seurat", demo.path))


demo <- subset(demo, downsample = 50000)
saveRDS(demo, demo.output)
