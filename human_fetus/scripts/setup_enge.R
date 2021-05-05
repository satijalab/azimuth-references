#!/usr/bin/env Rscript

library(Seurat)
library(magrittr)

# helper function for annotation
annotate <- function(obj, curr, new) {
  if (!is.list(curr)) curr <- list(curr)
  curr <- lapply(curr,function(vec){if (is.numeric(vec)) as.character(vec) else vec})
  
  stopifnot(length(curr)==length(new))
  new <- rep(new, times = sapply(curr, FUN = length))
  obj <- RenameIdents(obj, setNames(as.list(new), nm = unlist(curr)))
  return(obj)
}

# read files
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
files <- grep('GSM', list.files(path = path), value = TRUE)
cells <- NULL
for (f in files) {
  cells <- cbind(cells, read.table(file = file.path(path, f))[, 2])
}
rownames(x = cells) <- read.table(file = file.path(path, f))[, 1]
colnames(x = cells) <- paste0("cell_", 1:ncol(x = cells))
cells <- as.matrix(x = cells)
cells <- cells[1:(nrow(x = cells)-7), ] # not gene features

# construct object
obj <- CreateSeuratObject(counts = cells) %>% SCTransform %>% RunPCA %>% RunUMAP(dims = 1:50)
obj <- obj %>% FindNeighbors(dims = 1:50) %>% FindClusters(resolution = 0.8)

# annotate
curr <- list(c(12,3,5,4,21,6,0,18),17,c(2,24,13,23,25,8,20),c(16,1,19,14),c(9,11,7,10,22),c(15))
new <- c('alpha','delta','beta','acinar','ductal','stellate')
obj <- annotate(obj, curr, new)
# reannotate some clusters to capture endothelial
obj.subset <- subset(obj, idents = c('stellate', 'delta'))
obj.subset <- obj.subset %>% NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA %>% RunUMAP(dims = 1:10)
obj.subset <- obj.subset %>% FindNeighbors(dims = 1:50) %>% FindClusters(resolution = 2.5)
idents <- setNames(as.vector(Idents(obj)), Cells(obj))
idents[unlist(CellsByIdentities(obj.subset)[c('5','10')])] <- 'endothelial'
Idents(object = obj) <- idents
obj[['celltype']] <- Idents(object = obj)
DefaultAssay(object = obj)<-'RNA'
obj <- DietSeurat(obj, assays = 'RNA')

# save
saveRDS(object = obj, file = args[2])
