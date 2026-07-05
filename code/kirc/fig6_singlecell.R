## Fig 6 单细胞 ccRCC: A UMAP细胞类型 | B 二硫死亡活性(找热点细胞) | C 关键基因
suppressPackageStartupMessages({ library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(patchwork) })
set.seed(1)
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9",
         "MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11",
         "NDUFS1","OXSM","LRPPRC","NUBPL")
h5 <- list.files("data/GSE159115", pattern="\\.h5$", full.names=TRUE)
objs <- lapply(h5, function(f){
  s <- CreateSeuratObject(Read10X_h5(f), project=gsub(".h5","",basename(f)), min.cells=3, min.features=200)
  s$orig.ident <- gsub(".h5","",basename(f)); s[["percent.mt"]] <- PercentageFeatureSet(s,"^MT-")
  subset(s, subset=nFeature_RNA>200 & nFeature_RNA<6000 & percent.mt<25)
})
seu <- Reduce(function(a,b) merge(a,b), objs); seu <- JoinLayers(seu)
cat("细胞:", ncol(seu), "\n")
seu <- NormalizeData(seu)|>FindVariableFeatures(nfeatures=2000)|>ScaleData(verbose=FALSE)|>RunPCA(npcs=30,verbose=FALSE)
seu <- RunHarmony(seu,"orig.ident",verbose=FALSE)
seu <- RunUMAP(seu,reduction="harmony",dims=1:30,verbose=FALSE)
seu <- FindNeighbors(seu,reduction="harmony",dims=1:30,verbose=FALSE)|>FindClusters(resolution=0.4,verbose=FALSE)

mk <- list(Tumor_ccRCC=c("CA9","NDUFA4L2","VEGFA","NNMT","ANGPTL4"),
           ProxTubule=c("LRP2","CUBN","SLC34A1","GATM"),
           Endothelial=c("PECAM1","CLDN5","VWF","KDR"),
           Myeloid=c("LYZ","CD68","CD14","FCGR3A"),
           T_NK=c("CD3D","CD3E","NKG7","GNLY"),
           B_Plasma=c("MS4A1","CD79A","MZB1"),
           Fibroblast=c("PDGFRB","COL1A1","ACTA2","DCN"))
mk <- lapply(mk, intersect, rownames(seu))
seu <- AddModuleScore(seu, features=mk, name="ct_")
cs <- sapply(seq_along(mk), function(k) tapply(seu@meta.data[[paste0("ct_",k)]], seu$seurat_clusters, mean)); colnames(cs)<-names(mk)
cl2t <- setNames(colnames(cs)[apply(cs,1,which.max)], rownames(cs))
seu$cell_type <- unname(cl2t[as.character(seu$seurat_clusters)])
print(table(seu$cell_type))

du <- intersect(dis, rownames(seu))
seu <- AddModuleScore(seu, features=list(du), name="dis"); seu$Disulfidptosis <- seu$dis1

## A UMAP
pA <- DimPlot(seu, group.by="cell_type", label=TRUE, repel=TRUE) + ggtitle("A  ccRCC cell types") + theme(legend.position="none")
## B 二硫死亡活性
pB1 <- FeaturePlot(seu, "Disulfidptosis") + scale_color_gradientn(colors=c("#3C5488","white","#E64B35")) + ggtitle("B  Disulfidptosis activity")
ord <- names(sort(tapply(seu$Disulfidptosis, seu$cell_type, median), decreasing=TRUE))
seu$cell_type <- factor(seu$cell_type, ord)
pB2 <- VlnPlot(seu, "Disulfidptosis", group.by="cell_type", pt.size=0) + NoLegend() +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + ggtitle("")
## C 关键基因
pC <- DotPlot(seu, features=intersect(c("SLC7A11","TLN1","ACTB","MYH10","NCKAP1","FLNA","GYS1"),rownames(seu)), group.by="cell_type") +
  RotatedAxis() + ggtitle("C  Key disulfidptosis genes") + theme(axis.text.x=element_text(size=8))

fig <- (pA | pB1) / (pB2 | pC) + plot_layout(heights=c(1,.9))
ggsave("results/Fig6_singlecell.pdf", fig, width=14, height=11)
ggsave("results/Fig6_singlecell.png", fig, width=14, height=11, dpi=130)
cat("Fig6 完成 | 二硫死亡活性最高细胞:", ord[1], "\n")
saveRDS(seu, "results/kirc_scRNA.rds")
