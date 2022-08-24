# Demo 

library(Seurat)

args = commandArgs(trailingOnly=TRUE)
demo.dir <- args[1]
demo.output <- args[2]

Convert(demo.path, dest = "h5seurat", overwrite = TRUE)
demo <- LoadH5Seurat(gsub("h5ad", "h5seurat", demo.path))

demo <- subset(demo, downsample = 50000)
saveRDS(demo, demo.output)
