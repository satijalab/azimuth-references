library(Azimuth)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
integrated <- readRDS(args[1])

ob.ref <- AzimuthReference(
  object = integrated,
  refUMAP = "umap",
  refDR = "spca",
  refAssay = "SCT",
  metadata = c("celltype.l1", "celltype.l2"),
  dims = 1:50,
  reference.version = "1.0.0"
)

ValidateAzimuthReference(object = ob.ref)
SaveAnnoyIndex(
  object = ob.ref[["refdr.annoy.neighbors"]],
  file = args[2]
)
saveRDS(
  object = ob.ref,
  file = args[3]
)
