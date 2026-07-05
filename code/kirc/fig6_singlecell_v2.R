## Fig 6 升级版: 单细胞 UMAP + 二硫死亡活性 + 小提琴 + DotPlot,统一样式
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
seu <- readRDS("results/kirc_scRNA.rds")
ct_cols <- c("#3C5488","#E64B35","#00A087","#4DBBD5","#F39B7F","#8491B4","#B09C85")
names(ct_cols) <- levels(factor(seu$cell_type))[seq_len(length(unique(seu$cell_type)))]
pA <- DimPlot(seu, group.by="cell_type", cols=ct_cols, label=TRUE, repel=TRUE, label.size=2, pt.size=.2) +
  theme_nat() + ggtitle("ccRCC cell types") +
  theme(legend.position="none",
        axis.title.y=element_text(margin=margin(r=1)),
        axis.title.x=element_text(margin=margin(t=1)),
        plot.margin=margin(1,2,1,0))
pB <- FeaturePlot(seu, "Disulfidptosis", pt.size=.2) +
  scale_colour_gradientn(colours=c(PAL$low,"grey90",PAL$high), name="Activity") +
  theme_nat() + ggtitle("Disulfidptosis activity") +
  theme(axis.title.y=element_text(margin=margin(r=1)),
        axis.title.x=element_text(margin=margin(t=1)),
        plot.margin=margin(1,2,1,0))
ord <- names(sort(tapply(seu$Disulfidptosis, seu$cell_type, median), decreasing=TRUE))
seu$cell_type <- factor(seu$cell_type, ord)
cell_lab <- function(x) {
  x <- gsub("Fibroblast", "Fibro", x)
  x <- gsub("Endothelial", "Endo", x)
  x <- gsub("Tumor_ccRCC", "Tumor", x)
  x <- gsub("B_Plasma", "B/Plasma", x)
  x <- gsub("T_NK", "T/NK", x)
  x
}
pC <- VlnPlot(seu, "Disulfidptosis", group.by="cell_type", pt.size=0, cols=ct_cols[ord]) +
  scale_x_discrete(labels=cell_lab) +
  theme_nat() + NoLegend() +
  theme(axis.text.x=element_text(angle=0,hjust=.5,vjust=.5,size=5),
        plot.margin=margin(2, 8, 8, 3)) +
  ggtitle("Disulfidptosis by cell type") + xlab(NULL)
key <- intersect(c("SLC7A11","TLN1","ACTB","MYH10","NCKAP1","FLNA","GYS1"), rownames(seu))
pD <- DotPlot(seu, features=key, group.by="cell_type") +
  scale_colour_gradient(low=PAL$low, high=PAL$high) + theme_nat() +
  labs(size="Pct\nexpr.", colour="Avg\nexpr.") +
  guides(size=guide_legend(title.position="top", order=1),
         colour=guide_colorbar(title.position="top", order=2,
                               barheight=unit(12,"mm"), barwidth=unit(2.6,"mm"))) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_text(size=5.4),
        legend.text=element_text(size=5.2),
        legend.key.size=unit(2.6,"mm"),
        legend.spacing.y=unit(.5,"mm")) +
  ggtitle("Key disulfidptosis genes") + xlab(NULL) + ylab(NULL)
fig <- add_tags((pA | pB) / (pC | pD))
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig7_singlecell_v2", w_mm=183, h_mm=125)
cat("Fig6 升级版完成\n")
