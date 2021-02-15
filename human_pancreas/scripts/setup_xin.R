#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
data.dir <- "data/xin"

objs <- list()
for (i in 1:12) {
  donor <- paste0("Donor_", i)
  donor.dir <- file.path(data.dir, donor)
  if (!dir.exists(paths = donor.dir)) {
    files <- list.files(path = data.dir, pattern = paste0(donor, "_"))
    system2(command = "mkdir", args = , donor.dir)
    for (f in files) {
      system2(command = "mv", args = paste0("-f ", c(file.path(data.dir, f), file.path(donor.dir, Seurat:::ExtractField(string = f, field = 4)))))
      system2(command = "gunzip", args = file.path(donor.dir, Seurat:::ExtractField(string = f, field = 4)))
    }
  }
  counts <- Read10X(data.dir = donor.dir)
  objs[[i]] <- CreateSeuratObject(counts = counts) # already QC'd
  objs[[i]] <- SCTransform(object = objs[[i]])
  objs[[i]]$dataset <- donor
}
saveRDS(object = objs, file = args[1]) 
