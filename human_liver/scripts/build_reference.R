#!/usr/bin/env Rscript

# Script to build human liver reference 1.0.0
library(remotes)
library(Seurat)
library(SeuratDisk)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
aizarani.path <- args[1]
macparland.path <- args[2]
payen.path <- args[3]
ramachandran.path <- args[4]
zhang.path <- args[5]
annotations.path <- args[6]
ref.path <- args[7]
idx.path <- args[8]
obj.path <- args[9]

# Get Data
aizarani <- read.table(aizarani.path)
aizarani.obj <- CreateSeuratObject(aizarani)
aizarani.obj <- subset(aizarani.obj, subset = nCount_RNA > 1000)

macparland <- read.table(macparland.path, header = T, row.names=1, sep=",", as.is=T)
macparland.obj <- CreateSeuratObject(macparland)
macparland.obj[["percent.mt"]] <- PercentageFeatureSet(macparland.obj, pattern = "^MT-")
macparland.obj <- subset(macparland.obj, subset = nCount_RNA >= 1500 & percent.mt <= 50)
folders <- list.files(payen.path, full.names = T)
objs <- list()
for (i in 1:length(folders)){
  list <- list.files(folders[i])
  print(list)
  data <- Read10X(folders[i])
  objs[i] <- CreateSeuratObject(data)
}
payen.obj <- merge(objs[[1]], objs[2:length(objs)], add.cell.ids = c("GSM2", "GSM3", "GSM4", "GSM5", "GSM6", "GSM7", "GSM8", "GSM9"))
payen.obj[["percent.mt"]] <- PercentageFeatureSet(payen.obj, pattern = "^MT-")
payen.obj <- subset(payen.obj, subset = nCount_RNA > 2000 & percent.mt < 25)

folders <- list.files(ramachandran.path, full.names = T)
objs <- list()
for (i in 1:length(folders)){
  list <- list.files(folders[i])
  print(list)
  data <- Read10X(folders[i])
  objs[i] <- CreateSeuratObject(data)
}
ramachandran.obj <- merge(objs[[1]], objs[2:length(objs)])
ramachandran.obj[["percent.mt"]] <- PercentageFeatureSet(ramachandran.obj, pattern = "^MT-")
ramachandran.obj <- subset(ramachandran.obj, subset = nFeature_RNA >300 & percent.mt < 20)

adj_18 <- read.table(paste0(zhang.path, "/GSM4116579_ICC_18_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)
adj_23 <- read.table(paste0(zhang.path, "/GSM4116582_ICC_23_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)
adj_25 <- read.table(paste0(zhang.path, "/GSM4116586_ICC_25_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)

tables <- list(adj_18, adj_23, adj_25) 
names(tables) <- c("adjacent_18", "adjacent_23", "adjacent_25")
objs <- list()
for (i in 1:length(tables)){
  objs[i] <- CreateSeuratObject(tables[[i]])
}
zhang.obj <- merge(objs[[1]], objs[2:length(objs)])
zhang.obj[["percent.mt"]] <- PercentageFeatureSet(zhang.obj, pattern = "^MT-")
zhang.obj <- subset(zhang.obj, subset = nCount_RNA > 2001 & nFeature_RNA < 6001 & nFeature_RNA >501 & percent.mt < 20)
# Sketching 
objs <- c(aizarani.obj, macparland.obj, payen.obj, ramachandran.obj, zhang.obj)
names(objs) <- c("aizarani_v1", "macparland_v1", "payen_v2", "Ramachandran_healthy", "zhang_healthy")

atoms.list <- list()
for (i in 1:length(objs)) {
  object <- objs[[i]]
  dataset_name <- names(objs[i])
  object$dataset <- dataset_name
  object <- RenameCells(object = object, add.cell.id = dataset_name)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 3000)
  atoms.i <- LeverageScoreSampling(object = object, num.cells = 10000)
  atoms.list[[i]] <- atoms.i
}

# Integration 
for (i in 1:length(atoms.list)) {
  atoms.list[[i]] <- SCTransform(atoms.list[[i]], vst.flavor =2, verbose = T)
}
integration.features <- SelectIntegrationFeatures(object.list = atoms.list ,nfeatures = 3000)
atoms.list <- PrepSCTIntegration(object.list = atoms.list,anchor.features = integration.features)
atoms.anchors <- FindIntegrationAnchors(object.list = atoms.list, normalization.method = "SCT",anchor.features = integration.features)
atoms.merge <- IntegrateData(anchorset = atoms.anchors,normalization.method = "SCT",dims = 1:50) 
atoms.merge <- RunPCA(atoms.merge, verbose = FALSE)
atoms.merge <- RunUMAP(atoms.merge, reduction = "pca", dims = 1:30, return.model = TRUE)

# Expand to full object
integrated_objects <- list()
for (i in 1:length(objs)) {
  object <- objs[[i]]
  object$dataset <- names(objs[i])
  object <- RenameCells(object = object, add.cell.id = names(objs[i]))
  object <- NormalizeData(object)
  object <- IntegrateSketchEmbeddings(object = object, 
                                      atom.sketch.object = atoms.merge, 
                                      atom.sketch.reduction = "pca", 
                                      features = integration.features)
  integrated_objects[[i]] <- object
  rm(object)
}
obj.merge <- merge(integrated_objects[[1]], integrated_objects[2:length(integrated_objects)], merge.dr = "pca")
DefaultAssay(obj.merge) <- 'RNA'
obj.merge[['SCT']] <- CreateAssayObject(data = obj.merge[['RNA']]@data)
obj.merge[['SCT']] <- as(object = obj.merge[['SCT']] , Class = 'SCTAssay')
obj.merge[['SCT']]@SCTModel.list$ref.model <- atoms.list[[Seurat:::SampleIntegrationOrder(atoms.merge@tools$Integration@sample.tree)[1]]][['SCT']]@SCTModel.list$model1

# Annotate 
annotations <- read.csv(annotations.path, row.names = 1)
obj.merge <- AddMetaData(obj.merge, annotations)
obj.merge <- subset(obj.merge, subset = celltype.l1 != "NA")
obj.umap <- RunUMAP(obj.merge, reduction = "pca", dims = 1:30, return.model = TRUE)
saveRDS(obj.umap, obj.path)

# Azimuth
obj.azimuth <- AzimuthReference(obj.umap, refUMAP='umap', refDR='pca', plotref='umap', refAssay = 'SCT', metadata = c("celltype.l1", "celltype.l2"))
saveRDS(object = obj.azimuth, file =  ref.path, compress=F)
SaveAnnoyIndex(object = obj.azimuth[["refdr.annoy.neighbors"]], file = idx.path)