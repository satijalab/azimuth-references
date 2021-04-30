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

# read in intermediate object
library(Seurat)
library(stringr)
library(magrittr)
print("Reading object")
bbi <- readRDS(opt$input)
print("Done reading in")
df_gene <- readRDS(opt$`gene-metadata`)
df_cell <- readRDS(opt$`cell-metadata`)
rownames(df_cell) <- df_cell$sample

# annotate based on existing metadata
annotation.l1 <-
  str_match(
    bbi[['Organ_cell_lineage',drop=T]],
    '(?:(?:[A-Za-z]*)-(.*))|(.*)'
  )[,c(2,3)] %>% 
  t %>% 
  as.vector %>% 
  Filter(f=function(a)!is.na(a))
annotation.l2 <- bbi[['Organ_cell_lineage',drop=T]]
organ <- 
  str_match(
    bbi[['Organ_cell_lineage',drop=T]],
    '(?:([A-Za-z]*)-(?:.*))|(.*)'
  )[,2] 
bbi[['annotation.l1']] <- annotation.l1
bbi[['annotation.l2']] <- annotation.l2
bbi[['organ']] <- organ

# make azimuth reference
library(Azimuth)
ref <- 
  AzimuthReference(
    bbi,
    metadata = c('organ','annotation.l1','annotation.l2'), 
    refDR='pca',
    dims=1:100
  )

# subset plotref
ref@tools$AzimuthReference@plotref <- 
  subset(ref@tools$AzimuthReference@plotref, 
         cells = sample(Cells(ref), 10**5))

# convert gene names
old.names <- rownames(ref)
new.names <- df_gene[rownames(ref),]$gene_short_name
rownames(ref[['refAssay']]@data) <- new.names
rownames(ref[['refAssay']]@meta.features) <- new.names
rownames(ref[['refDR']]@feature.loadings) <- new.names
ref[['refAssay']]@SCTModel.list$refmodel@feature.attributes <- 
  ref[['refAssay']]@SCTModel.list$refmodel@feature.attributes[old.names,]
rownames(ref[['refAssay']]@SCTModel.list$refmodel@feature.attributes) <- new.names

# save
saveRDS(ref, 'out/ref.Rds')
SaveAnnoyIndex(ref[['refdr.annoy.neighbors']], 'out/idx.annoy')
saveRDS(bbi, 'out/fullref.Rds')

