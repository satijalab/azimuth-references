#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)

# Helper fxn
annotate <- function(obj, curr, new) {
  if (!is.list(x = curr)) curr <- list(curr)
  curr <- lapply(X = curr, FUN = function(vec) {
    if (is.numeric(x = vec)) as.character(x = vec) else vec
  })
  stopifnot(length(x = curr) == length(x = new))
  new <- rep(x = new, times = sapply(X = curr, FUN = length))
  obj <- RenameIdents(obj, setNames(as.list(x = new), nm = unlist(x = curr)))
  return(obj)
}

args <- commandArgs(trailingOnly = TRUE)
ref.dir <- "reference/"
ob.dir <- "seurat_objects/"

ob <- readRDS(file = args[1])
ob[['annotation.l2']] <- ob[['free_annotation']]
Idents(object = ob) <- 'annotation.l2'

l2.id <- list(
  c('Basophil/Mast 1','Basophil/Mast 2'),
  c('Ciliated','Proximal Ciliated'),
  c('Signaling Alveolar Epithelial Type 2','Alveolar Epithelial Type 2'),
  c('Bronchial Vessel 1','Bronchial Vessel 2'),
  c('Capillary Intermediate 1','Capillary Intermediate 2'),
  c('Proximal Basal','Proliferating Basal','Differentiating Basal'),
  c('TREM2+ Dendritic','IGSF21+ Dendritic','Myeloid Dendritic Type 1',
    'Myeloid Dendritic Type 2','EREG+ Dendritic'),
  c('Nonclassical Monocyte','Intermediate Monocyte'),
  c('Classical Monocyte','OLR1+ Classical Monocyte'),
  c('Airway Smooth Muscle','Vascular Smooth Muscle'),
  c('Mesothelial','Lipofibroblast'),
  c('Alveolar Fibroblast','Adventitial Fibroblast'),
  c('Myofibroblast','Fibromyocyte'),
  c('CD8+ Naive T','CD8+ Memory/Effector T'),
  c('CD4+ Memory/Effector T','CD4+ Naive T'),
  c('Neuroendocrine','Ionocyte')
)

l1.id <- c(
  'Basophil/Mast','Ciliated','Alveolar Epithelial Type 2','Bronchial Vessel',
  'Capillary Intermediate','Basal','Dendritic','CD16+ Monocyte','CD14+ Monocyte',
  'Smooth Muscle','Mesothelial','Fibroblast','Myofibroblast',
  'CD8 T','CD4 T','Ionocyte'
)

ob <- annotate(ob, l2.id, l1.id)
ob[['annotation.l1']] <- Idents(object = ob)

ref <- AzimuthReference(
  object = ob,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1", "annotation.l2"),
  dims = 1:50,
  k.param = 31,
  reference.version = "1.0.0"
)
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = ob, file = file.path(ob.dir, "fullref.Rds"))
