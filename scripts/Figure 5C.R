require(toolboxH)
library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(ggrepel)

#load diff. expr. Table
wilcoxtest = read_excel2(".../Data/s459_1_differential_expression_noDoublet.xlsx")



#define your genes of interest (GOI)
GOI = c("Sftpc", "Lcn2", "Fth1", "Cybb", "Ctsd", "Ctsb", "Mxd1", "Ywhah", "Il1b", "Nfil3", "Nfkbia", "Lyz2", "Cxcl3", "Cxcl2", "Parp9", "Zfp36", "Lst1", "Cyba", "Spi1", "Prdx5", "Cd52", "Tpt1", "Fau", "Atp1b1", "Hmgb2", "Il6ra", "Ifitm2", "Cd177", "Plac8", "Chil3")
GOI
Goitable = wilcoxtest[gene %in% GOI]
uniqueN(Goitable$gene)

# Filter for significant Genes
GOI_Sig_AB <- dplyr::filter(Goitable, c((Goitable$p_val_adj_tissue_sign=="TRUE") & (Goitable$contrast=="A vs. B")))
GOI_Sig_AB$Experssionvalue <- pmax(GOI_Sig_AB$meanexp_A, GOI_Sig_AB$meanexp_B)

range(GOI_Sig_AB$Experssionvalue)
range(GOI_Sig_AB$avg_log2FC)

order2 = c('B Cells', 'T Cells', 'NK Cells', 'Mesothelial Cells' , 'Fibroblasts', 'SMC', 'Pericytes','Endothelial Cells', 'ly. Endothelial Cells', 'AT1', 'AT2', 'Ciliated Cells', 'AM', 'MF & DC', 'Granulocytes')
GOI_Sig_AB$celltype <- factor(GOI_Sig_AB$celltype, levels = c(order2)) 
GOI_Sig_AB <- GOI_Sig_AB[order(GOI_Sig_AB$celltype),]
GOI_Sig_AB

ggplot(GOI_Sig_AB, aes(y=gene,x=celltype,size=Experssionvalue,color=avg_log2FC)) +
  geom_point() +
  scale_size_continuous(name='Mean Expression',
                        breaks=c(0.5,1,2,3,4,5),
                        labels=c("0.5","1", "2", "3", "4", "5"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-1,1),oob=scales::squish) +
    theme_minimal(base_size=10, base_line_size =0.1) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')

# Reverse Color
GOI_Sig_AB$BA_avg_log2FC <- GOI_Sig_AB$avg_log2FC * -1
range(GOI_Sig_AB$BA_avg_log2FC)

#flip Axis, no grid
p1 <- ggplot(GOI_Sig_AB, aes(x=gene,y=celltype,size=Experssionvalue,color=BA_avg_log2FC)) +
  geom_point() +
  scale_size_continuous(name='Mean Expression',
                        breaks=c(0.5,1,2,3,4,5),
                        labels=c("0.5","1", "2", "3", "4", "5"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-1,1),oob=scales::squish) +
  theme_minimal(base_size=10, base_line_size =0.001) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')

p2 <- p1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=10))
p2