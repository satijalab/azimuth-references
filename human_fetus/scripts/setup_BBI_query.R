#!/usr/bin/env Rscript

# parse args
args <- commandArgs(trailingOnly = TRUE)
library(optparse)
library(stringr)
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="count matrix file(s)", metavar="character"),
  make_option(c("-c", "--cell-metadata"), type="character", default=NULL,
              help="cell metadata RDS file path", metavar="character"),
  make_option(c("-g", "--gene-metadata"), type="character", default=NULL, 
              help="gene metadata RDS file path", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="reference file to exclude cells", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
opt$input <- str_split(opt$input,' ')[[1]]

print(opt$input)
# read matrix and metadata
library(Seurat)
library(stringr)
library(magrittr)
ref <- readRDS(opt$ref)
df_gene <- readRDS(opt$`gene-metadata`)
df_cell <- readRDS(opt$`cell-metadata`)
rownames(df_cell) <- df_cell$sample

objs <- list()
# load and subset
for (i in opt$input) {
  mat <- readRDS(i)
  mat <- mat[, sample(1:ncol(mat), size=min(ncol(mat),5000))]
  objs <- c(objs,mat)
}

# merge and make object
print(all(sapply(objs,function(o){all(rownames(o)==rownames(objs[[1]]))})))
merged.obj <- do.call(cbind, objs)
print(all(colnames(merged.obj)%in% rownames(df_cell)))
df_cell <- df_cell[colnames(merged.obj),]
merged.obj <- merged.obj[, which(!is.na(df_cell$Organ_cell_lineage))]
df_cell <- df_cell[which(!is.na(df_cell$Organ_cell_lineage)), ]
merged.obj <- CreateSeuratObject(merged.obj, meta.data = df_cell)

# exclude reference cells and save
merged.obj <- subset(merged.obj, cells=Cells(ref),invert=T)
merged.obj <- subset(merged.obj,features=which(rowSums(merged.obj[['RNA']]@counts)>0))
saveRDS(merged.obj,'out/bbi_query.rds')



