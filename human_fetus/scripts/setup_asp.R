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
library(Seurat)

counts <- read.csv(file = args[1], sep = '\t', row.names = 1)
metadata <- read.csv(file = args[2], sep = '\t')
metadata <- metadata[grep('\\(|\\)|/|\\&', metadata$X, invert = TRUE), ] # remove bad lines
rownames(x = metadata) <- metadata$X
metadata$X <- NULL
asp <- CreateSeuratObject(counts = counts, meta.data = metadata)
saveRDS(object = asp, file = args[3])
