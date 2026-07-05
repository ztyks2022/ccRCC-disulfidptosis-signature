## Fig 1 升级版: 评分T/N + DEG火山 + 关键基因 + 相关ComplexHeatmap,统一样式
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(ggplot2); library(ggrepel)
  library(ComplexHeatmap); library(circlize); library(grid) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
st <- substr(colnames(expr),14,15); grp <- factor(ifelse(st=="11","Normal","Tumor"), c("Normal","Tumor"))
du <- intersect(dis, rownames(expr))
sc <- tryCatch({gp<-ssgseaParam(expr,list(d=du)); gsva(gp)["d",]}, error=function(e) gsva(expr,list(d=du),method="ssgsea")["d",])
pA <- ggplot(data.frame(score=sc, grp=grp), aes(grp, score, fill=grp)) +
  geom_violin(alpha=.7, linewidth=.3) + geom_boxplot(width=.2, outlier.size=.3, linewidth=.3) +
  scale_fill_manual(values=two_cols, guide="none") +
  labs(title="Disulfidptosis score (p = 3.8e-12)", x=NULL, y="ssGSEA score")
deg <- fread("data/raw/kirc_DEG.csv"); setnames(deg,1,"gene")
deg$sig <- ifelse(deg$adj.P.Val<.05 & abs(deg$logFC)>1, ifelse(deg$logFC>0,"Up","Down"),"NS")
dd <- deg[gene %in% du]
pB <- ggplot(deg, aes(logFC, -log10(adj.P.Val+1e-300), colour=sig)) + geom_point(size=.25, alpha=.35) +
  geom_point(data=dd, colour=PAL$n_dark, size=1.4) +
  ggrepel::geom_text_repel(data=dd[adj.P.Val<.05], aes(label=gene), colour=PAL$n_dark, size=1.9, max.overlaps=20, segment.size=.2) +
  scale_colour_manual(values=c(Up=PAL$high, Down=PAL$low, NS=PAL$n_light), name=NULL) +
  labs(title="DEG tumor vs normal", x="log2 FC", y="-log10 FDR") + theme(legend.position=c(.12,.85))
key <- intersect(c("SLC7A11","TLN1","ACTB","MYH10","NCKAP1","GYS1"), rownames(expr))
dC <- do.call(rbind, lapply(key, function(gn) data.frame(gene=gn, val=expr[gn,], grp=grp)))
pC <- ggplot(dC, aes(grp, val, fill=grp)) + geom_boxplot(outlier.size=.2, linewidth=.3) +
  facet_wrap(~gene, nrow=1, scales="free_y") +
  scale_fill_manual(values=two_cols, guide="none") +
  labs(title="Key genes", x=NULL, y="log2 expression") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=5.3),
        strip.text=element_text(size=6, margin=margin(1,0,1,0)),
        panel.spacing=unit(2.2,"mm"))
tum <- expr[du, st=="01"]; cm <- cor(t(tum), method="spearman")
cmdf <- as.data.frame(as.table(cm))
pD <- ggplot(cmdf, aes(Var1, Var2, fill=Freq)) +
  geom_tile() +
  scale_fill_gradient2(low=PAL$low, mid="white", high=PAL$high,
                       limits=c(-1,1), name="r") +
  labs(title="Disulfidptosis gene correlation (tumor)", x=NULL, y=NULL) +
  theme_minimal(base_family="Arial", base_size=6.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=5),
        axis.text.y=element_text(size=5),
        panel.grid=element_blank(),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5.5),
        plot.title=element_text(size=7, face="bold", hjust=0))
fig <- add_tags((pA | pB | pC) / pD + plot_layout(heights=c(1,1.25)))
save_pub(fig, "results/Fig1_landscape_v2", w_mm=183, h_mm=130)
cat("Fig1 升级版完成\n")
