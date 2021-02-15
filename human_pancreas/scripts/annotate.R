#!/usr/bin/env Rscript
library(Matrix)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
pi <- readRDS(file = args[1])

Idents(object = pi) <- "integrated_snn_res.2"
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(0, 5, 7, 8, 9, 15, 27)), value = "alpha")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(1, 2, 4, 11, 13, 16, 19, 21, 26)), value = "beta")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(10, 18, 24, 30)), value = "acinar")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(3, 17, 25)), value = "ductal")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(6, 29)), value = "delta")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(14)), value = "gamma")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(34)), value = "epsilon")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(12)), value = "activated_stellate")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(31)), value = "schwann")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(20)), value = "endothelial")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(22)), value = "quiescent_stellate")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(28, 33)), value = "immune")
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(32)), value = "cycling")
# 23 mostly defined by mitochondrial
pi <- SetIdent(object = pi, cells = WhichCells(object = pi, idents = c(23)), value = "remove")
annotations <- Idents(object = pi)

saveRDS(annotations, file = args[2])

