#!/usr/bin/env Rscript

# parse args
args <- commandArgs(trailingOnly = TRUE)
library(optparse)
library(stringr)
option_list = list(
  make_option(c("-f", "--fullref"), type="character", default=NULL, 
              help="full object RDS file", metavar="character"),
  make_option(c("-g", "--dfgene"), type="character", default=NULL,
              help="gene metadata RDS file path", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

library(Seurat)
fullref <- readRDS(opt$fullref)
DefaultAssay(fullref) <- 'RNA'
fullref <- DietSeurat(
  fullref,
  assays = 'RNA',
  dimreducs = 'umap'
)
df_gene <- readRDS(opt$dfgene)

# add gene names
new.names <- df_gene[rownames(fullref),]$gene_short_name
notdup <- !duplicated(new.names)
new.names <- new.names[notdup]
fullref <- subset(fullref, features=rownames(fullref)[notdup])
rownames(fullref[['RNA']]@meta.features) <- new.names
rownames(fullref[['RNA']]@data) <- new.names
rownames(fullref[['RNA']]@counts) <- new.names

# save
saveRDS(fullref,'out/fullref_vitessce.Rds')