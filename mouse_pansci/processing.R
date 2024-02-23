# Script to process pansci downsampled dataset

# Convert to seurat object
# sceasy::convertFormat("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.h5ad", from="anndata", 
# to="seurat", outFile = "/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

a <- readRDS("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

a <- NormalizeData(a)
a <- FindVariableFeatures(a, nfeatures = 5000)
a <- ScaleData(a)
a <- RunPCA(a, npcs = 100)
a <- RunUMAP(a, dims = 1:100, return.model = T)
saveRDS(a, "/brahms/mollag/azimuth/pansci/pansci_processed_umap_100_dims.rds")

##### Post processing ####
a <- readRDS("/brahms/mollag/azimuth/pansci/pansci_processed_umap_100_dims.rds")

DimPlot(a, group.by = "Main_cell_type", label = T, repel = T, alpha = 0.1) + NoLegend()
DimPlot(a, group.by = "Lineage", label = T, repel = T, alpha = 0.1) + NoLegend()
DimPlot(a, group.by = "Organ_name", label = T, repel = T, alpha = 0.1) + NoLegend()




########### SCTransform VERSION #########
b <- readRDS("/brahms/mollag/azimuth/pansci/20240222_PanSci_all_cells_adata_sampled_by_subtype.rds")

b <- SCTransform(b, variable.features.n = 5000)
b <- ScaleData(b)
a <- RunPCA(a, npcs = 100)
a <- RunUMAP(a, dims = 1:100)

