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
              help="gene metadata RDS file path", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
opt$input <- str_split(opt$input,' ')[[1]]

# read matrix and metadata
library(Seurat)
library(stringr)
library(magrittr)
mat <- readRDS(opt$input)
df_gene <- readRDS(opt$`gene-metadata`)
df_cell <- readRDS(opt$`cell-metadata`)
rownames(df_cell) <- df_cell$sample

# prepare object
bbi <- CreateSeuratObject(mat, meta.data = df_cell[colnames(mat), ])
Idents(bbi) <- 'Organ_cell_lineage'

# normalize
bbi <- SCTransform(
  bbi,
  method = 'glmGamPoi', 
  clip.range = c(-10000, 10),
  do.correct.umi = F,
  conserve.memory = T
)

# dimensional reduction
bbi <- RunPCA(bbi, npcs=100)
bbi <- RunUMAP(bbi, dims=1:100, return.model=T)
saveRDS(bbi,'out/ref_intermediate.Rds',compress=F)
