#!/usr/bin/env Rscript
library(Seurat)
library(Matrix)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

dat <- data.table::fread(input = args[1])
cell.names <- dat[, 1]
dat[, 1] <- NULL
dat <- as(object = t(x = as.matrix(x = dat)), Class = "dgCMatrix")
colnames(x = dat) <- cell.names$sample_name
meta <- read.csv(file = args[2], row.names = 1)

ob <- CreateSeuratObject(counts = dat, meta.data = meta)
Idents(object = ob) <- "region_label"
ob <- subset(x = ob, idents = c("M1lm", "M1ul"))
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

saveRDS(ob, file = args[3])
