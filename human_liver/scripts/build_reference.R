#!/usr/bin/env Rscript

# Script to build human liver reference 1.0.0
library(remotes)
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
#remotes::install_github("satijalab/seurat", "feat/dictionary", lib = "/home/mollag/R/x86_64-pc-linux-gnu-library/4.2")
#myPaths <- .libPaths()
#myPaths <- c(myPaths, "/home/mollag/R/x86_64-pc-linux-gnu-library/4.2")
#.libPaths("/home/mollag/R/x86_64-pc-linux-gnu-library/4.2")
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

print("aizarani done")
macparland <- read.table(macparland.path, header = T, row.names=1, sep=",", as.is=T)
macparland.obj <- CreateSeuratObject(macparland)
print("macparland done")
folders <- list.files(payen.path, full.names = T)
objs <- list()
for (i in 1:length(folders)){
  list <- list.files(folders[i])
  print(list)
  data <- Read10X(folders[i])
  objs[i] <- CreateSeuratObject(data)
}
payen.obj <- merge(objs[[1]], objs[2:length(objs)], add.cell.ids = c("GSM2", "GSM3", "GSM4", "GSM5", "GSM6", "GSM7", "GSM8", "GSM9"))
print("payen done")

folders <- list.files(ramachandran.path, full.names = T)
objs <- list()
for (i in 1:length(folders)){
  list <- list.files(folders[i])
  print(list)
  data <- Read10X(folders[i])
  objs[i] <- CreateSeuratObject(data)
}
ramachandran.obj <- merge(objs[[1]], objs[2:length(objs)])
print("ramachandran done")
adj_18 <- read.table(paste0(zhang.path, "/GSM4116579_ICC_18_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)
adj_23 <- read.table(paste0(zhang.path, "/GSM4116582_ICC_23_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)
adj_25 <- read.table(paste0(zhang.path, "/GSM4116586_ICC_25_Adjacent_UMI.csv.gz"), sep=",", header = TRUE, row.names = 1)

tables <- list(adj_18, adj_23, adj_25) 
names(tables) <- c("adjacent_18", "adjacent_23", "adjacent_25")
objs <- list()
for (i in 1:length(tables)){
  print("step1")
  objs[i] <- CreateSeuratObject(tables[[i]])
  print("step2")
}
print("about to merge")
zhang.obj <- merge(objs[[1]], objs[2:length(objs)])
print("zhang done")
# Sketching 
objs <- c(aizarani.obj, payen.obj, macparland.obj, ramachandran.obj, zhang.obj)
names(objs) <- c("aizarani_v1", "macparland_v1", "payen_v2", "Ramachandran_healthy", "zhang_healthy")
print("obj list made")

atoms.list <- list()
for (i in 1:length(objs)) {
  print("step1")
  object <- objs[[i]]
  print(object)
  print("step2")
  dataset_name <- names(objs[i])
  object$dataset <- dataset_name
  print(object)
  print("step3")
  object <- RenameCells(object = object, add.cell.id = dataset_name)
  print(object)
  head(colnames(object))
  head(colnames(object[["RNA"]]@counts))
  print("step4")
  object <- NormalizeData(object)
  print("step5")
  object <- FindVariableFeatures(object)
  print("step6")
  atoms.i <- LeverageScoreSampling(object = object, num.cells = 10000)
  print("step7")
  atoms.list[[i]] <- atoms.i
}
 
# Integration 
for (i in 1:length(atoms.list)) {
  atoms.list[[i]] <- SCTransform(atoms.list[[i]], vst.flavor =2, verbose = T)
}
print("SCT done")
integration.features <- SelectIntegrationFeatures(object.list = atoms.list ,nfeatures = 3000)
print("found int features")
atoms.list <- PrepSCTIntegration(object.list = atoms.list,anchor.features = integration.features)
print("prepped sct")
atoms.anchors <- FindIntegrationAnchors(object.list = atoms.list, normalization.method = "SCT",anchor.features = integration.features)
print("found int anchors")
atoms.merge <- IntegrateData(anchorset = atoms.anchors,normalization.method = "SCT",dims = 1:50) 
print("integrated data")
atoms.merge <- RunPCA(atoms.merge, verbose = FALSE)
atoms.merge <- RunUMAP(atoms.merge, reduction = "pca", dims = 1:30, return.model = TRUE)
print("sketch integration done")

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
  print("added to list")
  rm(object)
}
print("scale up done")
obj.merge <- merge(integrated_objects[[1]], integrated_objects[2:length(integrated_objects)], merge.dr = "pca")
print("merged")
DefaultAssay(obj.merge) <- 'RNA'
print("default assay changed")
#obj.merge[['integrated']] <- NULL
obj.merge[['SCT']] <- CreateAssayObject(data = obj.merge[['RNA']]@data)
print("assay object created")
obj.merge[['SCT']] <- as(object = obj.merge[['SCT']] , Class = 'SCTAssay')
obj.merge[['SCT']]@SCTModel.list$ref.model <- atoms.list[[Seurat:::SampleIntegrationOrder(obj.merge@tools$Integration@sample.tree)[1]]][['SCT']]@SCTModel.list$model1
print("scale up done")

# Annotate 
annotations <- read.csv(annotations.path, row.names = 1)
obj.merge <- AddMetaData(obj.merge, annotations)
obj.merge <- subset(obj.merge, subset = celltype.l1 != "NA")
saveRDS(obj.merge, full.obj.path)

# Azimuth
obj.umap <- RunUMAP(obj.merge, reduction = "pca", dims = 1:30, return.model = TRUE)
obj.azimuth <- AzimuthReference(obj.umap, refUMAP='umap', refDR='pca', plotref='umap', refAssay = 'SCT', metadata = c("celltype.l1", "celltype.l2"))
saveRDS(object = obj.azimuth, file =  ref.path, compress=F)
SaveAnnoyIndex(object = obj.azimuth[["refdr.annoy.neighbors"]], file = annoy.path)