#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

panc8 <- readRDS(file = args[1])
panc.list <- SplitObject(object = panc8, split.by = "dataset")
panc.list <- lapply(X = panc.list, FUN = SCTransform)
saveRDS(object = panc.list, file = args[2])

