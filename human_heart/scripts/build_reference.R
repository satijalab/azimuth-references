#!/usr/bin/env Rscript

# Script to build human heart reference 1.0.0

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
remotes::install_github("satijalab/seurat", "feat/dictionary")
library(Seurat)
library(SeuratDisk)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
science.path <- args[1]
nature.path <- args[2]
naturecardio.path <- args[3]
annotations.path <- args[4]
ref.path <- args[5]
annoy.path <- args[6]
full.obj.path <- args[7]

# Get Data
science <- Read10X(science.path, cell.column = 2)
science.meta <- read.table(paste0(science.path, "/GSE165838_CARE_RNA_metadata.txt.gz"), row.names = 1)
science.obj <- CreateSeuratObject(science, min.cells = 3, min.features = 200, meta.data = science.meta, project = "Science") # min.cells = 3, min.features = 200

Convert(nature.path, dest = "h5seurat", overwrite = TRUE)
nature.obj <- LoadH5Seurat(gsub("h5ad", "h5seurat", nature.path))
nature.obj@meta.data$celltype.l1 <- nature.obj@meta.data$cell_type
nature.obj@meta.data$celltype.l2 <- nature.obj@meta.data$cell_states  
nature.obj <- subset(nature.obj, subset = celltype.l1 == "NotAssigned", invert = TRUE)
nature.obj <- subset(nature.obj, subset = celltype.l1 == "doublets", invert = TRUE)
nature.obj <- subset(nature.obj, subset = source == "CD45+", invert = TRUE) 
nature.obj <- subset(nature.obj, subset = percent_mito < .1)

nature_cardio.obj <- load(naturecardio.path, verbose = TRUE)
nature_cardio.obj <- RefMerge

# Integration
objs <- c(science.obj, nature.obj, nature_cardio.obj)
names(objs) <- c("science", "nature", "nature_cardio")
atoms.list <- list()

for (i in 1:length(objs)) {
  object <- objs[[i]]
  object$dataset <- names(objs[i])
  object <- RenameCells(object = object, add.cell.id = names(objs[i]))
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  atoms.i <- LeverageScoreSampling(object = object, num.cells = as.numeric(10000))
  atoms.list[[i]] <- atoms.i
}

names(atoms.list) <- c("science", "nature", "nature_cardio")
atoms.list <- unlist(atoms.list)

for (i in 1:length(atoms.list)) {
  atoms.list[[i]][["percent.mt"]] <- PercentageFeatureSet(atoms.list[[i]], pattern = "^MT-")
  atoms.list[[i]][["percent.rib"]] <- PercentageFeatureSet(atoms.list[[i]], pattern = "^RP[SL][[:digit:]]")
}

atoms.list[['nature']] <- SplitObject(atoms.list[['nature']], split.by = "source")
atoms.list[['nature_cardio']] <- SplitObject(atoms.list[['nature_cardio']], split.by = "tech")
atoms.list <- unlist(atoms.list)

ref_list <- c("science", "nature.Nuclei", "nature_cardio.SN")
ref_num <- c()
for (n in ref_list){
  num <- which(names(atoms.list) == n)
  ref_num <- c(ref_num, num)
}

for (i in 1:length(atoms.list)) {
  atoms.list[[i]] <- SCTransform(atoms.list[[i]], verbose = FALSE)
}

for (n in 1:length(atoms.list)){
  ribo.genes <- grep(pattern = "^[Rr][Pp][LlSs]", x = rownames(x = atoms.list[[n]]@assays$RNA), value = TRUE)
  mt.genes <- grep(pattern = "^[Mm][Tt]", x = rownames(x = atoms.list[[n]]@assays$RNA), value = TRUE)
  rm.genes <- c(ribo.genes, mt.genes)
  var.genes <- VariableFeatures(atoms.list[[n]])
  genes <- setdiff(var.genes, rm.genes)
  VariableFeatures(atoms.list[[n]]) <- genes
}

integration.features <- SelectIntegrationFeatures(object.list = atoms.list ,nfeatures = 3000)
atoms.list <- PrepSCTIntegration(object.list = atoms.list,anchor.features = integration.features)
atoms.anchors <- FindIntegrationAnchors(object.list = atoms.list, normalization.method = "SCT",anchor.features = integration.features, reference = ref_num)
atoms.merge <- IntegrateData(anchorset = atoms.anchors,normalization.method = "SCT",dims = 1:50) 
atoms.merge <- RunPCA(atoms.merge, verbose = FALSE)
atoms.merge <- RunUMAP(atoms.merge, reduction = "pca", dims = 1:30, return.model = TRUE)

# Expand to full objects
integrated_objects <- list()
for (i in 1:length(objs)) {
  object <- objs[[i]]
  dataset_name <- names(objs[i])
  object$dataset <- dataset_name
  object <- RenameCells(object = object, add.cell.id = dataset_name)
  object <- NormalizeData(object)
  object <- IntegrateSketchEmbeddings(object = object, 
                                      atom.sketch.object = atoms.merge, 
                                      atom.sketch.reduction = "pca", 
                                      features = integration.features)
  
  integrated_objects[[i]] <- object
  rm(object)
}

obj.merge <- merge(integrated_objects[[1]], integrated_objects[2:length(integrated_objects)], merge.dr = "pca")
obj.merge <- RunUMAP(obj.merge, reduction = "pca", dims = 1:30)

# Annotate 
annotations <- read.csv(annotations.path, row.names = 1)
obj.merge <- AddMetaData(obj.merge, annotations)
obj.merge <- subset(obj.merge, subset = celltype.l1 != "NA")
saveRDS(obj.merge, full.obj.path)

# Azimuth
DefaultAssay(obj.merge) <- "SCT" 
obj.merge <- SCTransform(obj.merge)
obj.umap <- RunUMAP(obj.merge, reduction = "pca", dims = 1:50, return.model = TRUE)
obj.azimuth <- AzimuthReference(obj.umap, refUMAP='umap', refDR='pca', plotref='umap', refAssay = 'SCT', metadata = c("celltype.l1", "celltype.l2"))
saveRDS(object = obj.azimuth, file =  ref.path, compress=F)
SaveAnnoyIndex(object = obj.azimuth[["refdr.annoy.neighbors"]], file = annoy.path)
