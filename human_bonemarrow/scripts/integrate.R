library(Seurat)
library(future)
library(future.apply)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = Inf)

set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)

nhlbi <- readRDS(args[1])
hca <- readRDS(args[2])
mpal <- readRDS(args[3])

# split hca by donor only
hca <- merge(x = hca[[1]], y = hca[2:length(hca)])
hca <- hca[, hca$nCount_RNA > 500]
hca <- SplitObject(object = hca, split.by = "donor")

nhlbi <- lapply(nhlbi, function(x) {
  x$orig.ident <- x$dataset
  x
})
hca <- lapply(hca, function(x) {
  x$orig.ident <- x$sample
  x
})
process_objects <- function(x) {
  x$percent.mt <- PercentageFeatureSet(object = x, pattern = "^MT-")
  x <- x[, x$percent.mt < 5]
  x <- SCTransform(x, variable.features.n = 3000)
  x <- RunPCA(x)
  return(x)
}

# sct integrate
DefaultAssay(mpal) <- "RNA"
obj.list <- SplitObject(mpal, split.by = "tissue")
obj.list <- future_lapply(obj.list, process_objects)
hca <- future_lapply(hca, process_objects)
nhlbi <- future_lapply(nhlbi, process_objects)
obj.list <- c(obj.list, nhlbi, hca)

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  reduction = "rpca",
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30
)

integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)

saveRDS(object = integrated, file = args[4])
