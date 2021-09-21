#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
data.dir <- "data/nhlbi"

objs <- list()

donors <- c("GSM3396161","GSM3396162","GSM3396163",
           "GSM3396164","GSM3396165","GSM3396166",
           "GSM3396167","GSM3396168","GSM3396169",
           "GSM3396170","GSM3396171","GSM3396172",
           "GSM3396173","GSM3396174","GSM3396175",
           "GSM3396176","GSM3396177","GSM3396178",
           "GSM3396179","GSM3396180","GSM3396181",
           "GSM3396182","GSM3396183","GSM3396184","GSM3396185")

for (i in 1:length(donors)) {
  donor <- donors[[i]]
  donor.dir <- file.path(data.dir, donor)
  if (!dir.exists(paths = donor.dir)) {
    files <- list.files(path = data.dir, pattern = paste0(donor, "_"))
    system2(command = "mkdir", args = , donor.dir)
    for (f in files) {
      system2(command = "mv", args = paste0("-f ", c(file.path(data.dir, f), 
                                                     file.path(donor.dir, gsub("_[a-zA-Z0-9]*", "", Seurat:::ExtractField(string = f, field = 2:3))))))
      system2(command = "gunzip", args = file.path(donor.dir, 
                                                   gsub("_[a-zA-Z0-9]*", "", Seurat:::ExtractField(string = f, field = 2:3))))
    }
  }
  counts <- Read10X(data.dir = donor.dir)
  objs[[i]] <- CreateSeuratObject(counts = counts)
  objs[[i]][["percent.mt"]] <- PercentageFeatureSet(objs[[i]], pattern = "^MT-")
  objs[[i]] <- subset(objs[[i]], subset = nFeature_RNA > 500  & percent.mt < 8) #use same qc metrics as original paper
  objs[[i]]$dataset <- donor
}
saveRDS(object = objs, file = args[1]) 

