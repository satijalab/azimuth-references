

# Convert to seurat object
# sceasy::convertFormat("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.h5ad", from="anndata", 
# to="seurat", outFile = "/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

a <- readRDS("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

a <- NormalizeData(a)
a <- FindVariableFeatures(a, nfeatures = 5000)
a <- ScaleData(a)
a <- RunPCA(a, npcs = 100)
a <- RunUMAP(a, dims = 1:100)

##### Post processing ####
a <- readRDS("/brahms/mollag/azimuth/pansci/pansci_processed_umap_100_dims.rds")

DimPlot(a, group.by = "Main_cell_type", label = T, repel = T, alpha = 0.1) + NoLegend()
DimPlot(a, group.by = "Lineage", label = T, repel = T, alpha = 0.1) + NoLegend()
DimPlot(a, group.by = "Organ_name", label = T, repel = T, alpha = 0.1) + NoLegend()
saveRDS(obj, full.reference)



# Build Azimuth Reference
azimuth.obj <- AzimuthReference(
  object = obj,
  refUMAP = "umap",
  plotref = "umap",
  refDR = "pca",
  refAssay = "RNA",
  metadata = c("celltype.l1", "celltype.l2"),
  dims = 1:100,
  reference.version = "2.0.0"
)

########### SCTransform VERSION #########
a <- readRDS("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

a <- SCTransform(a)
a <- FindVariableFeatures(a, nfeatures = 5000)
a <- ScaleData(a)
a <- RunPCA(a, npcs = 100)
a <- RunUMAP(a, dims = 1:100)



# they actually include "sub_umap" coordinates in the metadata? not sure if they computed the umap themselves? 
# Skylar looked into this and they dont look the best
