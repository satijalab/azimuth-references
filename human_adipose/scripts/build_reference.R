#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(yaml)
library(tidyverse)
library(stringr)
library(Azimuth)
library(DoubletFinder)

args = commandArgs(trailingOnly=TRUE)
rna.path <- args[1]
nuc.path <- args[2]
annotations.path <- args[3]
ref.path <- args[4]
annoy.path <- args[5]
full.obj.path <- args[6]

tsv_to_matrix <- function(file, sparse = TRUE, int = FALSE) {
  if (int) {
    default <- col_integer()
  } else {
    default <- col_double()
  }
  dat <- read_tsv(file, progress = FALSE, col_types = cols("GENE" = col_character(), .default = col_double()))
  dat <- as.matrix(column_to_rownames(data.frame(dat, check.names = FALSE), "GENE"))
  if (sparse) {
    dat <- as(dat, "dgCMatrix")
  }
  return(dat)
}
params <- yaml.load_file("scripts/pipeconfig.yaml")

###################### Prep #################################
for (i in 1:length(params$data)) {
  dge <- params$data[[i]]
  name <- names(params$data)[i]
  data <- tsv_to_matrix(dge$path, int = TRUE)
  if (length(dge) > 1) {
    metadata <- data.frame(dge[names(dge) != "path"])
    metadata <- metadata[rep(1, ncol(data)), ]
    row.names(metadata) <- colnames(data)
    new_seurat <- CreateSeuratObject(data, project = name, meta.data = metadata)}
  if (i == 1) {
    seurat <- new_seurat
  } else {
    seurat <- merge(seurat, new_seurat)
  }
}

rna.obj <- CreateSeuratObject(
  seurat[["RNA"]]@counts,
  project = params$project_name,
  min.cells = 2,
  meta.data = seurat@meta.data
)

rna.obj@meta.data[c('random', 'Number')] <- str_split_fixed(rna.obj@meta.data$individual, 's', 2)
rna.obj@meta.data <- subset(rna.obj@meta.data, select = -c(random, individual))
rna.obj@meta.data$species <- "Hs"
rna.list = SplitObject(rna.obj, split.by = "Number") 

nuc.data <- Read10X(data.dir = nuc.path)
nuc.obj <- CreateSeuratObject(counts = nuc.data, min.cells = 3, min.features = 200)
nuc.obj@meta.data[c('Patient', 'ID')] <- str_split_fixed(rownames(nuc.obj@meta.data), '-', 2)
nuc.obj@meta.data[c('Type', 'Fat', 'Number')] <- str_split_fixed(nuc.obj@meta.data$Patient, '_', 3)
nuc.obj@meta.data <- subset(nuc.obj@meta.data, select = -c(Patient, ID))

nuc.list <- SplitObject(nuc.obj, split.by = "Number")

# QC 

all.list <- c(nuc.list, rna.list)
all.list <- lapply(all.list, function (x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^[Mm][Tt]-")
  x[["percent.rib"]] <- PercentageFeatureSet(x, pattern = "^RP[SL][[:digit:]]")
  x
})

for (i in seq_len(length(all.list))){
  all.list[[i]] <- subset(all.list[[i]], subset = nCount_RNA > 400 & percent.mt < 10)
}

var_genes <- params$variable_genes
all.list <- lapply(X = all.list, FUN = SCTransform, vars.to.regress = "percent.mt")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for (i in seq_len(length(all.list))) {
  all.list[[i]] <- CellCycleScoring(all.list[[i]], s.features = s.genes, g2m.features = g2m.genes)
}

all.list <- lapply(X = all.list, FUN = SCTransform, variable.features.n = var_genes,vars.to.regress = c("percent.mt", "percent.rib", "S.Score", "G2M.Score"))

######################## Integration #############################

ref_list <- c("01", "12", "13", "254", "255", "256", "266") 
ref_num <- c()
for (n in ref_list){
  num <- which(names(all.list) == n)
  ref_num <- c(ref_num, num)
}

for (n in 1:length(all.list)){
  ribo.genes <- grep(pattern = "^[Rr][Pp][LlSs]", x = rownames(x = all.list[[n]]@assays$RNA), value = TRUE)
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = all.list[[n]]@assays$RNA), value = TRUE)
  rm.genes <- c(ribo.genes, mito.genes)
  var.genes <- VariableFeatures(all.list[[n]])
  genes <- setdiff(var.genes, rm.genes)
  VariableFeatures(all.list[[n]]) <- genes
}

# Reference Based Integration w/ RPCA 
integration.features <- SelectIntegrationFeatures(object.list = all.list,nfeatures = params$integrate$nfeatures)
all.list <- lapply(all.list,Seurat::RunPCA,features = integration.features,verbose = FALSE)
all.list <- PrepSCTIntegration(object.list = all.list,anchor.features = integration.features)
all.anchors <- FindIntegrationAnchors(object.list = all.list,reduction = "rpca",normalization.method = "SCT",anchor.features = integration.features, reference = ref_num)
all.combined <- IntegrateData(anchorset = all.anchors,normalization.method = "SCT",dims = 1:params$intdims, k.weight = 100)

all.combined <- RunPCA(all.combined, verbose = FALSE)
all.combined <- RunUMAP(all.combined, reduction = "pca", dims = 1:40)

########################## Annotation #########################################
annotations <- read.csv(annotations.path, row.names = 1)
all.combined <- AddMetaData(all.combined, annotations)
all.combined <- subset(all.combined, subset = celltype.l1 != "NA")
saveRDS(all.combined, full.obj.path)

################## Azimuth #####################################################
DefaultAssay(all.combined) <- "integrated"
human_azimuth <- RunUMAP(all.combined, reduction = "pca", dims = 1:50, return.model = TRUE)
human_ref <- AzimuthReference(human_azimuth, refUMAP='umap', refDR='pca', plotref='umap', refAssay = 'integrated', metadata = c("celltype.l1", "celltype.l2"))

saveRDS(object = human_ref, file =  ref.path, compress=F)
SaveAnnoyIndex(object = human_ref[["refdr.annoy.neighbors"]], file = annoy.path)
