library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(Nebulosa)

#load Seurat Object
mir223Seurat <- readRDS(".../Data/s451_7_seurat_ABX_vs_21-05-12.RDS")
mir223Seurat@active.ident = factor(mir223Seurat$seurat_clusters)

# Rename Clusters
mir223Seurat$Chosen_Cluster <- paste("C", mir223Seurat$seurat_clusters)
mir223Seurat$populations <- paste(mir223Seurat$Chosen_Cluster)
mir223Seurat@active.ident = factor(mir223Seurat$Chosen_Cluster)
mir223Seurat <- RenameIdents (mir223Seurat, 'C 0' = "B Cells", 'C 1' = "T Cells", 'C 2' = "Granulocytes", 'C 3' = "Fibroblasts", 'C 4' = "MF & DC", 'C 5' = "T Cells", 'C 6' = "AT2", 'C 7' = "NK Cells", 'C 8' = "Endothelial Cells", 'C 9' = "SMC", 'C 10' = "AM", 'C 11' = "Pericytes", 'C 12' = "Endothelial Cells", 'C 13' = "T Cells", 'C 14' = "Mesothelial Cells", 'C 15' = "T Cells", 'C 16' = "Ciliated Cells", 'C 17' = "ly. Endothelial Cells", 'C 18' = "AT1")
mir223Seurat$populations <- paste(mir223Seurat@active.ident)

UMAPPlot(mir223Seurat, group.by = "orig.ident", shuffle=T, cols = c("#ABAD59", "#CC5165", "#90BFF9")) + (coord_fixed(ratio=1)) 
