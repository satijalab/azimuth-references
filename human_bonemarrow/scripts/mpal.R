library(Seurat)


args <- commandArgs(trailingOnly = TRUE)

file.list <- list.files("data/mpal/", full.names = TRUE)
obj.list <- lapply(X = file.list, FUN = readRDS)
names(obj.list) <- sapply(file.list, function(x) substr(x, start = 23, stop = nchar(file.list[[1]])-4))

drop_na <- function(x) {
  na_row <- apply(X = x, MARGIN = 1, FUN = function(x) all(is.na(x)))
  return(x[!na_row, ])
}

bmmc_d1t1 <- CreateSeuratObject(counts = obj.list$scRNA_BMMC_D1T1, names.delim = ":")
bmmc_d1t1[['ADT']] <- CreateAssayObject(counts = drop_na(obj.list$scADT_BMMC_D1T1))

bmmc_d1t2 <- CreateSeuratObject(counts = obj.list$scRNA_BMMC_D1T2, names.delim = ":")
bmmc_d1t2[['ADT']] <- CreateAssayObject(counts = drop_na(obj.list$scADT_BMMC_D1T2))

cd34_d2t1 <- CreateSeuratObject(counts = obj.list$scRNA_CD34_D2T1, names.delim = ":")
counts.use <- drop_na(obj.list$scADT_CD34_D2T1)
cd34_d2t1 <- cd34_d2t1[, colnames(counts.use)]
cd34_d2t1[['ADT']] <- CreateAssayObject(counts = counts.use)

cd34_d3t1 <- CreateSeuratObject(counts = obj.list$scRNA_CD34_D3T1, names.delim = ":")
cd34_d3t1[['ADT']] <- CreateAssayObject(counts = drop_na(obj.list$scADT_CD34_D3T1))

obj <- list(bmmc_d1t1, bmmc_d1t2, cd34_d2t1, cd34_d3t1) #pbmc_d4t1, pbmc_d4t2) 

so <- merge(x = obj[[1]], y = obj[2:length(obj)])

so <- SCTransform(so)
so <- RunPCA(so)
so <- RunUMAP(so, reduction = "pca", dims = 1:30, return.model = TRUE)

DefaultAssay(so) <- "ADT"
so <- NormalizeData(so, normalization.method = "CLR")
so <- ScaleData(so, features = rownames(so))
so <- RunPCA(so, features = rownames(so), npcs = 20, reduction.name = "pca.adt")
so <- RunUMAP(so, reduction = 'pca.adt', dims = 1:20, reduction.name = "umap.adt", return.model = TRUE)
DimPlot(so, reduction = "umap.adt")

so$tissue <- sapply(strsplit(so$orig.ident, "_"), `[[`, 1)
saveRDS(so, file = args[1])
