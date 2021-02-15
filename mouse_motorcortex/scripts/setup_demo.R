#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

dat <- data.table::fread(input = args[1])
cell.names <- dat[, 1]
dat[, 1] <- NULL

meta <- read.csv(file = args[2], row.names = 1)
mop.cells <- rownames(meta)[which(meta$region_label == "MOp")]
mop.idx <- pmatch(x = mop.cells, table = cell.names$sample_name)

dat <- dat[mop.idx, ]
cell.names <- cell.names[mop.idx, ]
gc(verbose = FALSE)
dat <- as(object = t(x = as.matrix(x = dat)), Class = "dgCMatrix")
colnames(x = dat) <- cell.names$sample_name
meta <- meta[cell.names$sample_name, ]

ob <- CreateSeuratObject(counts = dat, meta.data = meta)
ob$subclass_label[ob$subclass_label == ""] <- NA
ob$subclass_color[ob$subclass_color == ""] <- NA
ob$subclass_order[ob$subclass_order == ""] <- NA

ob$class <- ob$class_label
ob$cluster <- ob$cluster_label
ob$subclass <- ob$subclass_label
ob$donor_sex <- ob$donor_sex_label
ob$region <- ob$region_label
ob$cortical_layer <- ob$cortical_layer_label
ob$external_donor_name <- ob$external_donor_name_label

metadata.keep <- c("class", "subclass", "cluster", "donor_sex", "region", "cortical_layer", "external_donor_name")
for (i in names(x = ob[[]])) {
  if (i %in% metadata.keep) {
    next
  } else {
    ob[[i]] <- NULL
  }
}
ob <- subset(x = ob, cells = names(which(is.na(ob$subclass))), invert = TRUE)
Idents(object = ob) <- "subclass"
ob <- subset(x = ob, downsample = 300)

saveRDS(ob, file = args[3])
