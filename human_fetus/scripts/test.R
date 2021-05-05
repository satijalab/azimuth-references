#!/usr/bin/env Rscript

# read matrix and metadata
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(stringr)
library(magrittr)
# ref <- readRDS(args[1])
print(args[2])
samples <- str_split(args[2], ' ')[[1]]
print(samples)
df_cell <- readRDS(args[3])
df_gene <- readRDS(args[4])
rownames(df_cell) <- df_cell$sample

objs <- list()
# load and subset
for (i in samples) {
  print(i)
}