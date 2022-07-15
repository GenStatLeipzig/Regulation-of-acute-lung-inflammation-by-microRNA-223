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

#Marker Dotplot (Fig. S4A)

Celltypeorder = c('AM', 'MF & DC', 'Granulocytes', 'NK Cells', 'T Cells', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Endothelial Cells', 'ly. Endothelial Cells', 'Fibroblasts', 'SMC', 'Pericytes', 'Mesothelial Cells')
my_levels <- c(Celltypeorder)
mir223Seurat@active.ident <- factor(x = mir223Seurat@active.ident, levels = my_levels)
p1<- DotPlot(mir223Seurat, features = c("Marco","Adgre1", "Flt3", "S100a8", "Siglecf", "Ncr1", "Cd3e", "Cd79a", "Epcam", "Akap5", "Lamp3", "Foxj1", "Cdh5", "Mmrn1", "Inmt", "Acta2", "Cox4i2", "Msln"))
p2 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
p2


#P Val. Diff expr. (Fig. S4B)
#for function reverselog_trans see: https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2

library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#see script "Figure 5C" to create GOI_Sig_AB
p3 <- ggplot(data=GOI_Sig_AB, aes(x=celltype, y=p_val_adj_tissue, label=gene)) + geom_point() + geom_text_repel() + theme_minimal() +  scale_y_continuous(trans=reverselog_trans(10))
p3

#Violin Plots of Granulocytes (Fig. S5A-S5C)
mir223Seurat$celltype_group <- paste(mir223Seurat$populations, mir223Seurat$orig.ident, sep = "_")
Idents(mir223Seurat) <- "celltype_group"
subset_PNM_BA <- subset(x = mir223Seurat, idents = c("Granulocytes_A_wt_d39wt", "Granulocytes_B_mir_d39wt"))
UMAPPlot(subset_PNM_BA)

VlnPlot(subset_PNM_BA, features = "Spi1", split.by = "orig.ident", cols = c("#ABAD59","#CC5165"))
VlnPlot(subset_PNM_BA, features = "Cxcl2", split.by = "orig.ident", cols = c("#ABAD59","#CC5165"))
VlnPlot(subset_PNM_BA, features = "Stat3", split.by = "orig.ident", cols = c("#ABAD59","#CC5165"))
VlnPlot(subset_PNM_BA, features = "Mmp8", split.by = "orig.ident", cols = c("#ABAD59","#CC5165"))
VlnPlot(subset_PNM_BA, features = "Nfkbiz", split.by = "orig.ident", cols = c("#ABAD59","#CC5165"))
