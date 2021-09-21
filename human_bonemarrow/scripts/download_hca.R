if (!requireNamespace(package = "SeuratDisk", quietly = TRUE)) {
  if (!requireNamespace(package = "remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("mojaveazure/seurat-disk")
}
if (!requireNamespace(package = "curl", quietly = TRUE)) {
  install.packages("curl")
}
library(Seurat)
library(SeuratDisk)
library(curl)
library(future)
library(future.apply)
plan(strategy = "multicore", workers = 12)
options(future.globals.maxSize = Inf)

args <- commandArgs(trailingOnly = TRUE)

md <- read.table(args[1], sep = "\t", header = FALSE)
info <- list()
for (i in seq_len(nrow(md))) {
  info[[i]] <- list(
    url = md$V17[i],
    donor = md$V38[i],
    sex = md$V37[i],
    age = md$V43[i],
    sample = md$V52[i]
  )
}

obj.list <- future_lapply(
  X = info,
  FUN = function(i) {
    fname <- paste0("data/hca/", i$sample, ".loom")
    curl_download(url = i$url, destfile = fname, quiet = FALSE)
    x <- Connect(filename = fname, mode = "r")
    x.seurat <- as.Seurat(x = x)
    x$close()
    x.seurat <- x.seurat[, x.seurat$emptydrops_IsCell == 1]
    x.seurat$sex <- i$sex
    x.seurat$age <- i$age
    x.seurat$donor <- i$donor
    x.seurat$sample <- i$sample
    return(x.seurat)
  }
)

saveRDS(object = obj.list, args[2])