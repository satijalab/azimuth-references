#!/usr/bin/env Rscript
library(Seurat)
library(SeuratObject)
library(Azimuth)
library(Matrix)
# used to read parquet files which efficiently compress CSVs
library(arrow)

args = commandArgs(trailingOnly=TRUE)
counts.path <- args[1]
annotations.path <- args[2]
dr.path <- args[3]
ref.path <- args[4]
annoy.path <- args[5]
full.obj.path <- args[6]

mtx <- readRDS(counts.path)
obj <- CreateSeuratObject(counts = mtx)

# load annotations
annotations <- read_parquet(annotations.path)
annotations <- as.data.frame(annotations)[,c("ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level")]
rownames(annotations) <- Cells(obj)
obj <- AddMetaData(obj, metadata = annotations)
print(head(obj))

# load in the scANVI latent space which contains 30 dimensions
latent.space <- as.matrix(read_parquet(dr.path))
rownames(latent.space) <- Cells(obj)
scanvi.dr <- CreateDimReducObject(embeddings = as.matrix(latent.space), key = "SCANVI")
obj[["scanvi"]] <- scanvi.dr

# find neighbors based on scANVI latent space
obj <- FindNeighbors(obj, reduction = "scanvi")

# Run SCTransform on the raw counts
obj <- SCTransform(obj, method = "glmGamPoi")

# run sPCA
obj <- RunSPCA(object = obj, assay = "SCT", graph = "RNA_snn")

# Force RunUMAP to run with n.epochs to prevent RunUMAP from running for 0 epochs
# Related: https://scanpy.discourse.group/t/umap-incorrectly-installed/663
obj <- RunUMAP(obj, dims = 1:30, reduction = "scanvi", n.epochs = 200, return.model = TRUE)

# save the full size object to perform marker identification
saveRDS(obj, file = full.obj.path)

annotations <- c("ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level")
ref <- AzimuthReference(obj,
                        refUMAP = 'umap',
                        refDR = 'spca',
                        dims = 1:50, # use 50 dimensions from the sPCA dimensional reduction
                        plotref = 'umap',
                        reference.version = '2.0.1',
                        metadata = annotations)

saveRDS(object = ref, file = ref.path, compress = F)
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]],
               file = annoy.path)
