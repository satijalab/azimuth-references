# MAIN link: https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/

# mkdir -p data/asp_data
# wget 'https://data.mendeley.com/datasets/mbvhhf8m62/2' -O data/asp_data
# rm data/asp_data/2 # junk
# unzip data/asp_data/Developmental_heart_scRNA-seq.zip
# gunzip data/asp_data/share_files/*
# mv data/asp_data/share_files/all_cells_count_matrix_filtered.tsv data/asp_data/counts_filtered.tsv
# mv data/asp_data/share_files/all_cells_meta_data_filtered.tsv data/asp_data/md_filtered.tsv


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
counts <- read.csv(opt$`counts-file`,sep='\t',row.names = 1)
metadata <- read.csv(opt$`metadata-file`,sep='\t')
metadata <- metadata[grep('\\(|\\)|/|\\&',metadata$X,invert=T),] # remove bad lines
rownames(metadata) <- metadata$X; metadata$X <- NULL
asp <- CreateSeuratObject(counts, meta.data = metadata)
saveRDS(asp, 'out/asp.rds')