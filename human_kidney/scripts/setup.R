#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

#################### Packages ##################################################
library(Seurat)
library(dplyr)
library(Matrix)
library(Azimuth)

#################### Helper Fxns ###############################################
pagoda_batch_correction <- function(
  batch,
  countMatrix,
  object
) {
  counts <- t(x = countMatrix) # cell x gene now
  cell.depth <- colSums(x = countMatrix)
  gene.av <- (colSums(x = counts) + length(x = unique(x = batch)))/
    (sum(cell.depth) + length(x = unique(x = batch)))
  gene.sum.perbatch <- t(
    x = sapply(X = levels(x = batch), function(b) {
      print(x = b)
      colSums(x = counts[which(batch == b), ])
    }))
  depth.per.batch <- as.numeric(x = tapply(X = cell.depth, INDEX = batch, FUN = sum))
  # normalize each batch
  gene.sum.perbatch <- t(
    x = log(x = gene.sum.perbatch + 1) - 
      log(x = depth.per.batch + 1)
  )
  # "center" each gene
  gene.sum.perbatch <- 
    exp(x = gene.sum.perbatch - log(x = gene.av))
  gene.id <- rep(1:counts@Dim[2], diff(counts@p))
  # divide counts by what was computed above
  counts@x <- as.numeric(
    x = counts@x / gene.sum.perbatch[cbind(gene.id, as.integer(x = batch)[counts@i + 1])]
  )
  depthScale <- 1000
  trim <- 10
  if (trim == 0) { # we clip genes in ScaleData, but we never clip cells
    counts <- counts / as.numeric(x = cell.depth)
    pagoda2:::inplaceWinsorizeSparseCols(sY = counts, n = trim)
    counts <- counts * as.numeric(x = cell.depth)
    cell.depth <- round(x = Matrix::rowSums(x = counts))
  }
  counts <- counts/as.numeric(x = cell.depth/depthScale)
  counts@x <- as.numeric(x = log(x = counts@x + 1))
  obj2 <- CreateSeuratObject(
    counts = t(x = counts), 
    meta.data = object[[]] 
  ) 
  obj2 <- SCTransform(object = obj2, clip.range = c(-10, 10))
  obj2 <- RunPCA(object = obj2) 
  obj2 <- RunUMAP(object = obj2, dims = 1:50, return.model = TRUE)
  return(obj2)
}

LoadData <- function(mtx, features, cells, metadata, study) {
  mat <- readMM(file = mtx)
  features <- read.table(file = features)
  rownames(x = mat) <- features[, 1]
  cells <- read.table(file = cells)
  colnames(x = mat) <- cells[, 1]
  metadata <- read.csv(file = metadata, row.names = 1)
  ob <- CreateSeuratObject(counts = mat, meta.data = metadata)
  ob$study <- study
  return(ob)
}

################## Processing ##################################################
kpmp <- LoadData(
  mtx = args[1],
  features = args[2],
  cells = args[3],
  metadata = args[4],
  study = "KPMP"
)
kpmp_pagoda <- pagoda_batch_correction(
  batch = as.factor(x = kpmp[['batch', drop = TRUE]]),
  countMatrix = GetAssayData(object = kpmp[['RNA']], slot = "counts"),
  object = kpmp
) 

lake <- LoadData(
  mtx = args[5],
  features = args[6],
  cells = args[7],
  metadata = args[8],
  study = "Lake"
)
lake_pagoda <- pagoda_batch_correction(
  batch = as.factor(x = lake[['library', drop = TRUE]]),
  countMatrix = GetAssayData(object = lake[['RNA']], slot = "counts"),
  object = lake
) 

## merge kpmp and lake kidney
Idents(object = kpmp_pagoda) <- Idents(object = kpmp)
Idents(object = lake_pagoda) <- 'annotation.l3'
feats <- SelectIntegrationFeatures(
  object.list = list(kpmp_pagoda, lake_pagoda), 
  nfeatures = 5000
)
kidney.merged <- merge(x = kpmp_pagoda, y = lake_pagoda)
VariableFeatures(object = kidney.merged) <- feats
kidney.merged <- RunPCA(object = kidney.merged, npcs = 100)
kidney.merged <- RunUMAP(object = kidney.merged, dims = 1:100)

merged.nuconly <- subset(
  x = kidney.merged, 
  subset = study == 'Lake'
)
merged.nuconly <- RunUMAP(
  object = merged.nuconly, 
  dims = 1:100, 
  return.model = TRUE
)
merged.nuconly[['SCT']]@SCTModel.list$model1.1 <- NULL

annotation.l3 <- lake$subclass.full
annotation.l3 <- sapply(
  X = annotation.l3, 
  FUN = function(x) {
    gsub(" Cell", "", as.character(x = x))
  }
)
merged.nuconly$celltype <- annotation.l3

# convert annotations to readable levels
levels.l2 <- unlist(lapply(sort(unique(merged.nuconly$annotation.l2)),function(x) names(which.max(table(annotation.l3,merged.nuconly$annotation.l2)[,x]))))
names(levels.l2) <- sort(unique(merged.nuconly$annotation.l2))

translate.l2 <- data.frame(levels.l2)
translate.l2["DCT",1] <- 'Distal Convoluted Tubule'
annotation.l2 <- translate.l2[merged.nuconly$annotation.l2,1]

new.l1 <- c("Ascending Thin Limb","Connecting Tubule","Distal Convoluted Tubule","Descending Thin Limb","Endothelial", "Fibroblast", "Intercalated", "Immune","Schwann", "Papillary Tip Epithelial",
            "Principal","Parietal Epithelial", "Podocyte", "Proximal Tubule", "Thick Ascending Limb", "Vascular Smooth Muscle / Pericyte")
translate.l1 <- data.frame(new.l1,row.names =sort(unique(merged.nuconly$annotation.l1)))
annotation.l1 <- translate.l1[merged.nuconly$annotation.l1,1]
merged.nuconly$annotation.l1 <- factor(annotation.l1,levels = sort(unique(annotation.l1)))
merged.nuconly$annotation.l2 <- factor(annotation.l2,levels = sort(unique(annotation.l2)))
merged.nuconly$annotation.l3 <- factor(annotation.l3,levels = sort(unique(annotation.l3)))

saveRDS(object = merged.nuconly, file = args[9])
