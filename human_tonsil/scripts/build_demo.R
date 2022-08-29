#!/usr/bin/env Rscript

# Script to build demo dataset for human tonsil reference 1.0.0

library(Seurat)

args = commandArgs(trailingOnly=TRUE)
demo.counts.dir <- args[1]
demo.cell.calls.dir <- args[2]
demo.output <- args[3]

demo.counts.grep <- paste0("count/outs/raw")
counts.paths <- file.path(demo.counts.dir, list.dirs(demo.counts.dir, full.names = FALSE)[grepl(demo.counts.grep, list.dirs(demo.counts.dir, full.names = FALSE))])
cell.calls.paths <- file.path(demo.cell.calls.dir, list.files(demo.cell.calls.dir))

seurat.objs <- list()
for (i in 1:length(counts.paths)) {
    raw.counts <- counts.paths[[i]]
    cell.calls.path <- cell.calls.paths[[i]]
    cell.calls <- read.csv(cell.calls.path)[, 2]
    counts <- Read10X(raw.counts)
    filtered.counts <- counts[["Gene Expression"]][, cell.calls]
    seurat.objs[[i]] <- CreateSeuratObject(filtered.counts)
}

obj <- merge(seurat.objs[[1]], seurat.objs[2:length(counts.paths)])
saveRDS(obj, demo.output)
