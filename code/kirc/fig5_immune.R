## Fig 5 免疫微环境(风险高 vs 低组): A 免疫细胞浸润 | B 检查点 | C 免疫/基质评分
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(ggplot2); library(patchwork); library(scales); library(tidyr); library(dplyr) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df          # sample, risk, rgrp
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
tum <- expr[, df$sample]; df$Risk <- factor(df$rgrp, c("Low","High"))

## 免疫细胞 signature(ssGSEA)
imm <- list(CD8_T=c("CD8A","CD8B","GZMK"), CD4_T=c("CD4","IL7R","CD40LG"), Treg=c("FOXP3","IL2RA","IKZF2"),
  Th1=c("TBX21","IFNG","STAT1"), Bcell=c("CD19","MS4A1","CD79A"), Plasma=c("MZB1","IGHG1","XBP1"),
  NK=c("KLRD1","NKG7","GNLY","KLRF1"), DC=c("CD1C","CLEC9A","LILRA4"), Macrophage=c("CD68","CSF1R"),
  M1=c("NOS2","CXCL9","CXCL10"), M2=c("CD163","MRC1","MSR1"), Monocyte=c("CD14","FCN1","S100A8"),
  Neutrophil=c("FCGR3B","CSF3R"), Mast=c("TPSAB1","CPA3","MS4A2"), Cytotoxic=c("GZMB","PRF1","GNLY"))
imm <- lapply(imm, intersect, rownames(tum)); imm <- imm[sapply(imm,length)>=2]
isc <- tryCatch({gp<-ssgseaParam(tum,imm); gsva(gp)}, error=function(e) gsva(tum,imm,method="ssgsea"))
isc <- as.data.frame(t(isc)); isc$Risk <- df$Risk
long <- pivot_longer(isc, -Risk, names_to="cell", values_to="val")
pv <- long %>% group_by(cell) %>% summarise(p=wilcox.test(val~Risk)$p.value, .groups="drop") %>% mutate(sig=p<0.05)
pA <- ggplot(long, aes(reorder(cell, val), val, fill=Risk)) + geom_boxplot(outlier.size=.15, lwd=.3) +
  scale_fill_manual(values=c("#3C5488","#E64B35")) + coord_flip() + theme_bw() +
  labs(title="A  Immune cell infiltration (ssGSEA)", x=NULL, y="score") +
  geom_text(data=pv[pv$sig,], aes(x=cell, y=max(long$val)*.95, label="*"), inherit.aes=FALSE, size=5)

## B 免疫检查点
ckp <- intersect(c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","IDO1","CD276","VTCN1","SIGLEC15"), rownames(tum))
cdf <- do.call(rbind, lapply(ckp, function(gn) data.frame(gene=gn, val=tum[gn,], Risk=df$Risk)))
pvc <- cdf %>% group_by(gene) %>% summarise(p=wilcox.test(val~Risk)$p.value, .groups="drop")
cdf$gene <- factor(cdf$gene, levels=ckp)
pB <- ggplot(cdf, aes(gene, val, fill=Risk)) + geom_boxplot(outlier.size=.15, lwd=.3) +
  scale_fill_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="B  Immune checkpoints", x=NULL, y="log2 expression") + theme(axis.text.x=element_text(angle=45,hjust=1))

## C 免疫/基质评分(ESTIMATE 式 ssGSEA)
es <- list(ImmuneScore=intersect(c("PTPRC","CD3D","CD2","CD8A","CD19","NKG7","CD68","IL2RB","CCL5","GZMA","CXCL9","CXCL10"),rownames(tum)),
           StromalScore=intersect(c("COL1A1","COL1A2","COL3A1","FAP","PDGFRB","ACTA2","DCN","LUM","THBS2","FN1","SPARC"),rownames(tum)))
esc <- tryCatch({gp<-ssgseaParam(tum,es); gsva(gp)}, error=function(e) gsva(tum,es,method="ssgsea"))
edf <- as.data.frame(t(esc)); edf$Risk <- df$Risk
elong <- pivot_longer(edf, -Risk, names_to="set", values_to="val")
pe <- elong %>% group_by(set) %>% summarise(p=wilcox.test(val~Risk)$p.value, .groups="drop")
elong <- left_join(elong, pe, by="set"); elong$set <- sprintf("%s (p=%.1e)", elong$set, elong$p)
pC <- ggplot(elong, aes(Risk, val, fill=Risk)) + geom_violin(alpha=.6) + geom_boxplot(width=.2, outlier.size=.2) +
  facet_wrap(~set, scales="free_y") + scale_fill_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="C  Immune / Stromal score", x=NULL, y="ssGSEA") + theme(legend.position="none")

fig <- (pA | pC) / pB + plot_layout(heights=c(1.1,1))
ggsave("results/Fig6_immune.pdf", fig, width=14, height=10)
ggsave("results/Fig6_immune.png", fig, width=14, height=10, dpi=130)
cat("Fig5 完成 | 显著差异免疫细胞:", paste(pv$cell[pv$sig], collapse=","), "\n")
