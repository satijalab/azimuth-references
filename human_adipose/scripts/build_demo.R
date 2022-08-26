#!/usr/bin/env Rscript

# Script to build human adipose demo

library(Seurat)

args = commandArgs(trailingOnly=TRUE)
demo.dir <- args[1]
demo.output <- args[2]

files <- list.files(demo.dir, pattern = "*.csv", full.names = TRUE)
tables <- lapply(files, read.table, sep=",", header = TRUE, row.names = 1)
names(tables) <- list.files(demo.dir, pattern = "*.csv",full.names=FALSE)
seurat_objs <- list()
for (i in 1:length(tables)){
  seurat_objs[i] <- CreateSeuratObject(tables[[i]])
}

demo.merge <- merge(seurat_objs[[1]], seurat_objs[2:length(seurat_objs)])
demo <- subset(demo.merge, downsample = 50000)

saveRDS(demo, demo.output)