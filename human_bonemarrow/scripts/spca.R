library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

integrated <- readRDS(file = args[1])
annot <- read.table(file = args[2], sep = ",", row.names = 1, header = TRUE)

integrated <- integrated[, rownames(annot)]
integrated <- AddMetaData(integrated, metadata = annot)

# compute PCA, UMAP and graph based on integrated PCA
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:50, return.model = TRUE, n.neighbors = 30L)
integrated <- FindNeighbors(integrated, dims = 1:50, reduction = "pca")

# sctransform
DefaultAssay(integrated) <- "RNA"
integrated <- SCTransform(
  object = integrated,
  residual.features = rownames(integrated[['integrated']]),
  conserve.memory = TRUE
)

# spca
integrated <- RunSPCA(
  object = integrated,
  graph = "integrated_snn",
  assay = "SCT"
)

saveRDS(integrated, args[3])
