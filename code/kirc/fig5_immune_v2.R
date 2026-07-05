## Fig 5 升级版: 免疫浸润 + 检查点 + 免疫/基质评分,统一样式
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(ggplot2); library(tidyr); library(dplyr) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
tum <- expr[, df$sample]; df$Risk <- factor(df$rgrp, c("Low","High"))
imm <- list(CD8_T=c("CD8A","CD8B","GZMK"), CD4_T=c("CD4","IL7R","CD40LG"), Treg=c("FOXP3","IL2RA","IKZF2"),
  Th1=c("TBX21","IFNG","STAT1"), Bcell=c("CD19","MS4A1","CD79A"), Plasma=c("MZB1","IGHG1","XBP1"),
  NK=c("KLRD1","NKG7","GNLY","KLRF1"), DC=c("CD1C","CLEC9A","LILRA4"), Macrophage=c("CD68","CSF1R"),
  M1=c("NOS2","CXCL9","CXCL10"), M2=c("CD163","MRC1","MSR1"), Monocyte=c("CD14","FCN1","S100A8"),
  Neutrophil=c("FCGR3B","CSF3R"), Mast=c("TPSAB1","CPA3","MS4A2"), Cytotoxic=c("GZMB","PRF1","GNLY"))
imm <- lapply(imm, intersect, rownames(tum)); imm <- imm[sapply(imm,length)>=2]
isc <- tryCatch({gp<-ssgseaParam(tum,imm); gsva(gp)}, error=function(e) gsva(tum,imm,method="ssgsea"))
isc <- as.data.frame(t(isc)); isc$Risk <- df$Risk
long <- pivot_longer(isc, -Risk, names_to="cell", values_to="val")
pv <- long %>% group_by(cell) %>% summarise(p=wilcox.test(val~Risk)$p.value, .groups="drop")
pA <- ggplot(long, aes(reorder(cell,val), val, fill=Risk)) + geom_boxplot(outlier.size=.12, linewidth=.22) +
  scale_fill_manual(values=risk_cols) + coord_flip() + labs(title="Immune cell infiltration (ssGSEA)", x=NULL, y=NULL) +
  geom_text(data=pv[pv$p<.05,], aes(x=cell, y=max(long$val)*.99, label="*"), inherit.aes=FALSE, size=2.2) +
  theme(legend.position=c(.015,.985),
        legend.justification=c(0,1),
        legend.background=element_rect(fill=scales::alpha("white", .82), colour=NA),
        legend.key.size=unit(2.4, "mm"),
        legend.title=element_text(size=5.5),
        legend.text=element_text(size=5.2),
        plot.margin=margin(3, 5, 2, 3))
ckp <- intersect(c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","IDO1","CD276"), rownames(tum))
cdf <- do.call(rbind, lapply(ckp, function(gn) data.frame(gene=factor(gn,ckp), val=tum[gn,], Risk=df$Risk)))
pB <- ggplot(cdf, aes(gene, val, fill=Risk)) + geom_boxplot(outlier.size=.12, linewidth=.22) +
  scale_fill_manual(values=risk_cols) + labs(title="Immune checkpoints", x=NULL, y=expression(log[2]*" expression")) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none",
        plot.margin=margin(3, 5, 2, 3))
es <- list(ImmuneScore=intersect(c("PTPRC","CD3D","CD2","CD8A","CD19","NKG7","CD68","IL2RB","CCL5","GZMA","CXCL9","CXCL10"),rownames(tum)),
           StromalScore=intersect(c("COL1A1","COL1A2","COL3A1","FAP","PDGFRB","ACTA2","DCN","LUM","THBS2","FN1","SPARC"),rownames(tum)))
esc <- tryCatch({gp<-ssgseaParam(tum,es); gsva(gp)}, error=function(e) gsva(tum,es,method="ssgsea"))
edf <- as.data.frame(t(esc)); edf$Risk<-df$Risk; elong <- pivot_longer(edf,-Risk,names_to="set",values_to="val")
pe <- elong %>% group_by(set) %>% summarise(p=wilcox.test(val~Risk)$p.value,.groups="drop")
elong <- left_join(elong,pe,by="set"); elong$set <- sub("Score"," score", elong$set)   # 短strip标签,p值移到面板内(避免facet标题截断)
plab <- elong %>% group_by(set) %>% summarise(p=first(p), y=max(val)+0.06*diff(range(val)), .groups="drop")
pC <- ggplot(elong, aes(Risk,val,fill=Risk)) + geom_violin(alpha=.7,linewidth=.3) + geom_boxplot(width=.2,outlier.size=.2,linewidth=.3) +
  geom_text(data=plab, aes(x=1.5,y=y,label=sprintf("p=%.1e",p)), inherit.aes=FALSE, size=2, family="Arial") +
  facet_wrap(~set, scales="free_y") + scale_fill_manual(values=risk_cols, guide="none") + labs(title="Immune / Stromal score", x=NULL, y="ssGSEA")
fig <- add_tags((pA | pC) / pB + plot_layout(heights=c(1.1,1)))
save_pub(fig, "results/Fig6_immune_v2", w_mm=183, h_mm=125)
cat("Fig5 升级版完成\n")
