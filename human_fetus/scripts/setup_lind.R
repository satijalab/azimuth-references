#https://pubmed.ncbi.nlm.nih.gov/29449449/
# mkdir -p data/lind_data
# wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102596&format=file' -O data/lind_data/
# mv data/lind_data/'index.html?acc=GSE102596&format=file' lind.tar
# tar -xvf data/lind_data/lind.tar
# gunzip data/lind_data/GSM2741551_count-table-human16w.tsv.gz
# mv data/lind_data/GSM2741551_count-table-human16w.tsv data/lind_data/counts.tsv

args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(magrittr)
counts <- read.csv(file = args[1], sep = '\t', row.names = 1)
lind <- CreateSeuratObject(counts = counts)

#process
lind <- lind %>% NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA %>% RunUMAP(dims = 1:50)
lind <- lind %>% FindNeighbors(dims = 1:50) %>% FindClusters(resolution = 2)

# # annotate based on gene markers given in Fig. 7
# FeaturePlot(lind, features=toupper(c('plvap','gng11','tie1')))
# FeaturePlot(lind, features=toupper(c('ccl3','srgn')))
# FeaturePlot(lind, features=toupper(c('top2a','hist1h4c')))
# FeaturePlot(lind, features=toupper(c('lypd1','six1','dapl1','cited1')))
# FeaturePlot(lind, features=toupper(c('cdh6','emx2','lrp2')))
# FeaturePlot(lind, features=toupper(c('mal','aldh1a1','krt18','gata3')))
# FeaturePlot(lind, features=toupper(c('ren','mgp','tcf21','aldh1a2')))
# FeaturePlot(lind, features=toupper(c('vcam1','fabp2','tagln','meis1')))
# FeaturePlot(lind, features=toupper(c('lum','sfrp2','dcn')))

# helper function for annotation
annotate <- function(obj, curr, new) {
  if (!is.list(curr)) curr <- list(curr)
  curr <- lapply(curr,function(vec){if (is.numeric(vec)) as.character(vec) else vec})
  
  stopifnot(length(curr)==length(new))
  new <- rep(new, times = sapply(curr, FUN = length))
  obj <- RenameIdents(obj, setNames(as.list(new), nm = unlist(curr)))
  return(obj)
}
curr <- list(2,c(10,20),18,c(16,17),c(4,9,5,14),
          c(7,19,21),c(0,1,8,3,11,13,6,12,15))
new <- c('Developing vasculature','Immune','Erythroid','Cell-cycling','Nephron progenitors',
          'Differentiating nephron cell-types','Interstitial lineage')
lind <- annotate(lind, curr, new) 
lind[['celltype']] <- Idents(lind)

# save 
saveRDS(object = lind, file = args[2]) 
