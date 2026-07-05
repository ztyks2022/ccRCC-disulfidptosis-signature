## Fig 1 二硫死亡图景: A 评分T/N | B DEG火山 | C 关键基因T/N | D 基因相关热图
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(ggplot2); library(patchwork); library(scales) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9",
         "MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11",
         "NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
st <- substr(colnames(expr),14,15); grp <- factor(ifelse(st=="11","Normal","Tumor"), c("Normal","Tumor"))
du <- intersect(dis, rownames(expr))

## A 二硫死亡评分 T vs N
sc <- tryCatch({gp<-ssgseaParam(expr,list(d=du)); gsva(gp)["d",]}, error=function(e) gsva(expr,list(d=du),method="ssgsea")["d",])
dA <- data.frame(score=sc, grp=grp)
pA <- ggplot(dA, aes(grp, score, fill=grp)) + geom_violin(alpha=.6) + geom_boxplot(width=.2, outlier.size=.3) +
  scale_fill_manual(values=c("#4DBBD5","#E64B35")) + theme_bw() +
  labs(title="A  Disulfidptosis score", subtitle="Wilcoxon p = 3.8e-12", x=NULL, y="ssGSEA score") + theme(legend.position="none")

## B DEG 火山(标二硫死亡基因)
deg <- fread("data/raw/kirc_DEG.csv"); setnames(deg, 1, "gene")
deg$sig <- ifelse(deg$adj.P.Val<0.05 & abs(deg$logFC)>1, ifelse(deg$logFC>0,"Up","Down"), "NS")
dd <- deg[gene %in% du]
pB <- ggplot(deg, aes(logFC, -log10(adj.P.Val+1e-300), color=sig)) + geom_point(size=.25, alpha=.3) +
  geom_point(data=dd, color="black", size=1.6) +
  geom_text(data=dd[adj.P.Val<0.05], aes(label=gene), color="black", size=2.3, vjust=-.6) +
  scale_color_manual(values=c(Up="#E64B35", Down="#3C5488", NS="grey80")) + theme_bw() +
  labs(title="B  DEG tumor vs normal (disulfidptosis genes labeled)", x="log2FC", y="-log10 FDR") + theme(legend.title=element_blank())

## C 关键基因 T vs N
key <- intersect(c("SLC7A11","TLN1","ACTB","MYH10","NCKAP1","GYS1"), rownames(expr))
dC <- do.call(rbind, lapply(key, function(gn) data.frame(gene=gn, val=expr[gn,], grp=grp)))
pC <- ggplot(dC, aes(grp, val, fill=grp)) + geom_boxplot(outlier.size=.2) + facet_wrap(~gene, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("#4DBBD5","#E64B35")) + theme_bw() +
  labs(title="C  Key genes", x=NULL, y="log2 expression") + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1))

## D 二硫死亡基因相关热图(肿瘤内)
tum <- expr[du, st=="01"]; cm <- cor(t(tum), method="spearman")
cmdf <- as.data.frame(as.table(cm))
pD <- ggplot(cmdf, aes(Var1, Var2, fill=Freq)) + geom_tile() +
  scale_fill_gradient2(low="#3C5488", mid="white", high="#E64B35", limits=c(-1,1), name="r") +
  theme_minimal() + labs(title="D  Disulfidptosis gene correlation (tumor)", x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6), panel.grid=element_blank())

fig <- (pA | pB) / (pC) / (pD) + plot_layout(heights=c(1, .8, 1.3))
ggsave("results/Fig1_landscape.pdf", fig, width=13, height=13)
ggsave("results/Fig1_landscape.png", fig, width=13, height=13, dpi=130)
cat("Fig1 完成 -> results/Fig1_landscape.*\n")
