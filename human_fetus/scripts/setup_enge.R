# mkdir -p data/enge_data
# wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81547&format=file' -O data/enge_data/enge.tar
# tar -xvf data/enge_data/enge.tar

# parse args
args <- commandArgs(trailingOnly = TRUE)
library(optparse)
library(stringr)
option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to directory with data", metavar="character"),
  make_option(c("-f", "--counts-file"), type="character", default=NULL,
              help="path to file with counts data", metavar="character"),
  make_option(c("-m", "--metadata-file"), type="character", default=NULL,
              help="path to metadata file", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

library(Seurat)
library(magrittr)
library(readr)
# helper function for annotation
annotate <- function(obj, curr, new) {
  if (!is.list(curr)) curr <- list(curr)
  curr <- lapply(curr,function(vec){if (is.numeric(vec)) as.character(vec) else vec})
  
  stopifnot(length(curr)==length(new))
  new <- rep(new, times = sapply(curr, FUN = length))
  obj <- RenameIdents(obj, setNames(as.list(new), nm = unlist(curr)))
  return(obj)
}

# read files
path <- opt$path
if (str_sub(path,nchar(path)) != "/") {
  path <- paste0(path,'/')
}
files <- grep('GSM', system(paste0("ls ",path),intern=T), value=T)
cells <- NULL
for (f in files) {
  cells <- cbind(cells,(readr::read_delim(paste0(path,f),delim='\t',col_names = F) %>% as.data.frame)[,2])
}
rownames(cells) <- (readr::read_delim(paste0(path,f),delim='\t',col_names = F) %>% as.data.frame)[,1]
colnames(cells) <- paste0("cell_",1:ncol(cells))
cells <- as.matrix(cells)
cells <- cells[1:(nrow(cells)-7),] # not gene features

# construct object
obj <- CreateSeuratObject(counts=cells) %>% SCTransform %>% RunPCA %>% RunUMAP(dims=1:50)
obj <- obj %>% FindNeighbors(dims=1:50) %>% FindClusters(resolution=0.8)

# annotate
curr <- list(c(12,3,5,4,21,6,0,18),17,c(2,24,13,23,25,8,20),c(16,1,19,14),c(9,11,7,10,22),c(15))
new <- c('alpha','delta','beta','acinar','ductal','stellate')
obj <- annotate(obj,curr,new)
# reannotate some clusters to capture endothelial
obj.subset <- subset(obj, idents = c('stellate','delta'))
obj.subset <- obj.subset %>% NormalizeData %>% FindVariableFeatures %>% ScaleData %>% RunPCA %>% RunUMAP(dims=1:10)
obj.subset <- obj.subset %>% FindNeighbors(dims=1:50) %>% FindClusters(resolution=2.5)
idents <- setNames(as.vector(Idents(obj)), Cells(obj))
idents[unlist(CellsByIdentities(obj.subset)[c('5','10')])] <- 'endothelial'
Idents(obj)<-idents
obj[['celltype']] <- Idents(obj)
DefaultAssay(obj)<-'RNA'
obj <- DietSeurat(obj,assays = 'RNA')

# save
saveRDS(obj,"out/enge.rds")
