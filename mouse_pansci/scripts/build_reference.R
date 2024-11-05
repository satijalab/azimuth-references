#!/usr/bin/env Rscript

# Script to build mouse pansci reference 1.0.0

# Rscript {input.script} {input.counts} {input.genenames} {input.dr} {output.ref} {output.idx} {output.obj}

library(Seurat)
library(Azimuth)

args = commandArgs(trailingOnly=TRUE)
obj_original <- readRDS(args[1])
genenames <- read.csv(args[2], row.names = 1)
azimuth.reference <- args[3]
azimuth.index <- args[4]
full.reference <- args[5]

# Convert gene IDs to gene symbols based on annotations file 
pansci.mtx <- LayerData(obj_original, assay="RNA", layer="counts")
gene_names <- read.csv(genenames)
rownames(pansci.mtx) <- gene_names$gene_name
subset.pansci.mtx <- pansci.mtx[grep("^ENSMUSG", rownames(pansci.mtx), invert = TRUE), ]

# Rebuild object from RNA assay and converted gene names
pansci.assay <- CreateAssayObject(counts = subset.pansci.mtx)
obj <- CreateSeuratObject(counts=pansci.assay, assay="RNA")
obj <- AddMetaData(obj, obj_original[[]])

# Normalize wtih SCT
obj <- SCTransform(obj, variable.features.n = 5000)

# Run PCA on the SCT residuals 
obj <- RunPCA(obj, assay = "SCT", reduction.name = "pca", npcs = 100)

# Build UMAP model 
obj <- RunUMAP(obj, dims = 1:100, return.model = T)

saveRDS(obj, full.reference)

# Build Azimuth Reference
azimuth.obj <- AzimuthReference(obj,
                        refUMAP='umap',
                        refDR='pca',
                        refAssay = 'SCT',
                        dims = 1:100,
                        plotref='umap',
                        metadata = c("Main_cell_type","Sub_cell_type", "Sub_cell_type_organ"),
                        reference.version = "1.0.0"
)

ValidateAzimuthReference(object = azimuth.obj)

SaveAnnoyIndex(
  object = azimuth.obj[["refdr.annoy.neighbors"]],
  file = azimuth.index
)
saveRDS(
  object = azimuth.obj,
  file = azimuth.reference
)