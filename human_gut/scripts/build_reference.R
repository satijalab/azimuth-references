#!/usr/bin/env Rscript

# Script to build human liver reference 1.0.0
library(remotes)
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
#remotes::install_github("satijalab/seurat", "feat/dictionary", lib = "/home/mollag/R/x86_64-pc-linux-gnu-library/4.2")
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(data.table)
library(patchwork)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
data.path <- args[1]
#annotations.path <- args[6]
#ref.path <- args[7]
#idx.path <- args[8]
obj.path <- args[2]

# Get Data
dataset_name = "Nature_2021"
Convert(data.path, dest = "h5seurat", overwrite = TRUE)
print("conversion done ")
obj <- LoadH5Seurat(gsub("h5ad", "h5seurat", data.path))
print("loading done")
obj = subset(obj, subset = total_counts <= 60000)
obj <- RenameCells(object = obj, add.cell.id = dataset_name)
# i removed the step of normalizing data here since you do it again for the sketch 
# save the whole QCed dataset
#saveRDS(obj, file = "Seurat_obj_Gut_atlas_Nature_2021.rds")


# split and sketch based on sources: "Pediatric healthy" "Pediatric Crohn Disease" "Healthy adult" "fetal""
Idents(obj) = "Diagnosis"
new.list = SplitObject(obj)
# and save each subset into a seurat obj
data_name_list = c("Pediatric", 
                   "Pediatric_Crohn",
                   "Adult", 
                   "Fetal")

for (i in 1:length(new.list)) {
  saveRDS(new.list[[i]], file = paste0("data/subset_Seurat_obj_", 
                                       data_name_list[i], 
                                       "_Gut_atlas_Nature_2021.rds"))
}
# empty list to store the sketched subsets
atoms.list = list()  

# run sketching
for (i in 1:length(new.list)) {
  if(i == 2){
    # skipping Crohn-disease subset 
    next()
  }
  object <- new.list[[i]]
  object$dataset <- data_name_list[i]
  object <- RenameCells(object = object, add.cell.id = data_name_list[i])
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
  # calculate leverage score and sample 5k/10k cells based on leverage score
  if(i == 1){
    atoms.i <- LeverageScoreSampling(object = object, num.cells = 5000)
  } else if (i %in% c(3, 4)) {
    atoms.i <- LeverageScoreSampling(object = object, num.cells = 10000)
  }
  atoms.list[[i]] <- atoms.i
  rm(object)
}
atoms.list[[2]] = NULL # essential, or there will be a NULL object in atoms.list[[2]]

#saveRDS(atoms.list, file = "subset_sketched_Seurat_obj_5k_10k_10k__Gut_atlas_Nature_2021.rds")

# run standard multiCCA integration
#for (i in 1:length(atoms.list)) { - you shouldnt need to do this bc you already ran find variable features, just going to add these parameters to the sketched part
#  atoms.list[[i]] <- FindVariableFeatures(atoms.list[[i]], selection.method = "vst", nfeatures = 3000)
#}

# run find anchors
features <- SelectIntegrationFeatures(object.list = atoms.list, nfeatures = 3000)
atoms.anchors <- FindIntegrationAnchors(object.list = atoms.list, anchor.features = features, k.anchor = 10, reference = c(2, 3))  # use dataset 2 and 3 as reference and map dataset 1 onto them
atoms.combined <- IntegrateData(anchorset = atoms.anchors)

# Run the standard workflow for visualization and clustering
DefaultAssay(atoms.combined) <- "integrated"
atoms.combined <- ScaleData(atoms.combined, verbose = FALSE)
atoms.combined <- RunPCA(atoms.combined, npcs = 100, verbose = FALSE)

# use the top 100 PCs 
atoms.combined <- RunUMAP(atoms.combined, reduction = "pca", dims = 1:100)

#saveRDS(atoms.combined, file = "data/v2_Sketched_integrated_Seurat_obj_5k_10k_10k_ref23_anchor10_Gut_atlas_Nature_2021.rds")

data_name_list = c("Pediatric", 
                   "Adult", 
                   "Fetal")

# empty list to store the scaled-up datasets
integrated_objects <- list()
for (i in 1:3) {
  dataset_name <- data_name_list[i]
  # load in Seurat object / basic preprocessing
  object <- readRDS(file = paste0("data/subset_Seurat_obj_", 
                                  data_name_list[i], 
                                  "_Gut_atlas_Nature_2021.rds"))
  object$dataset <- dataset_name
  object <- RenameCells(object = object, add.cell.id = dataset_name)
  object <- NormalizeData(object)
  
  # Integrate all cells into the same space as the atoms
  object <- IntegrateSketchEmbeddings(object = object, 
                                      atom.sketch.object = atoms.combined, 
                                      atom.sketch.reduction = "pca",
                                      features = features)
  integrated_objects[[i]] <- object
  rm(object)
}

# integrate the 3 scaled-up datasets
obj.merge <- merge(integrated_objects[[1]], integrated_objects[2:length(integrated_objects)], merge.dr = "pca")
obj.merge <- RunUMAP(obj.merge, reduction = "pca", dims = 1:100)

# 
obj.merge = subset(obj.merge, subset = nFeature_RNA > 0)

saveRDS(obj.merge, obj.path)


# Annotate 
#annotations <- read.csv(annotations.path, row.names = 1)
#obj.merge <- AddMetaData(obj.merge, annotations)
#obj.merge <- subset(obj.merge, subset = celltype.l1 != "NA")
#saveRDS(obj.merge, full.obj.path)

# Azimuth
#obj.umap <- RunUMAP(obj.merge, reduction = "pca", dims = 1:30, return.model = TRUE)
#obj.azimuth <- AzimuthReference(obj.umap, refUMAP='umap', refDR='pca', plotref='umap', refAssay = 'SCT', metadata = c("celltype.l1", "celltype.l2"))
#saveRDS(object = obj.azimuth, file =  ref.path, compress=F)
#SaveAnnoyIndex(object = obj.azimuth[["refdr.annoy.neighbors"]], file = annoy.path)