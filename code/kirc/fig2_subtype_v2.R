## Fig 2 升级版: ComplexHeatmap 分型热图 + 统一样式
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(data.table)
  library(survival); library(survminer); library(ggplot2); library(tidyr); library(dplyr); library(grid) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
mg <- readRDS("data/raw/kirc_subtype_immune.rds"); mg$sub <- factor(mg$sub)
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
du <- intersect(dis, rownames(expr)); ord <- order(mg$sub)
mat <- t(scale(t(expr[du, mg$sample[ord]])))
af <- gpar(fontsize=6, fontfamily="Arial")
ta <- HeatmapAnnotation(Subtype=mg$sub[ord], col=list(Subtype=c(C1=PAL$low, C2=PAL$high)),
  annotation_name_gp=af, simple_anno_size=unit(3,"mm"),
  annotation_legend_param=list(title_gp=af, labels_gp=af))
ht <- Heatmap(mat, name="z-score", col=div_ramp(), top_annotation=ta,
  cluster_columns=FALSE, show_column_names=FALSE, row_names_gp=gpar(fontsize=5, fontfamily="Arial"),
  column_title="Disulfidptosis genes by subtype", column_title_gp=gpar(fontsize=7, fontface="bold", fontfamily="Arial"),
  heatmap_legend_param=list(title_gp=af, labels_gp=af, legend_height=unit(14,"mm")))
pA <- wrap_elements(grid.grabExpr(draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")))

gg <- ggsurvplot(survfit(Surv(OS.time/365,OS)~sub, data=mg), data=mg, palette=c(PAL$low,PAL$high),
  pval=TRUE, pval.size=2.4, size=.55, censor.size=1.6, legend.title="Subtype")
pB <- gg$plot + theme_nat() + labs(title="Subtype survival", x="Time (years)", y="Overall survival") +
  theme(legend.position=c(.78,.86))

imm <- mg %>% select(sub, CD8_Tcell, Cytolytic, Checkpoint, Treg) %>% pivot_longer(-sub, names_to="set", values_to="val")
pv  <- imm %>% group_by(set) %>% summarise(p=wilcox.test(val~sub)$p.value, .groups="drop")
imm$set <- factor(sub("CD8_Tcell","CD8 T",imm$set), levels=c("CD8 T","Cytolytic","Checkpoint","Treg"))  # 改单面板分组箱线,x轴放全名,消除facet截断
pv$set <- factor(sub("CD8_Tcell","CD8 T",pv$set), levels=c("CD8 T","Cytolytic","Checkpoint","Treg"))
star <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","ns")))
plab <- imm %>% group_by(set) %>% summarise(y=max(val)+0.04*diff(range(val)), .groups="drop") %>% left_join(pv,by="set")
pC <- ggplot(imm, aes(set, val, fill=sub)) +
  geom_boxplot(outlier.size=.2, linewidth=.3, width=.7, position=position_dodge(.8)) +
  geom_text(data=plab, aes(x=set, y=y, label=star(p)), inherit.aes=FALSE, size=2.6) +
  scale_fill_manual(values=c(PAL$low,PAL$high), name="Subtype") +
  labs(title="Immune infiltration by subtype", x=NULL, y="ssGSEA") +
  theme(legend.position="right", legend.key.size=unit(3,"mm"))
pD <- ggplot(mg, aes(sub, dscore, fill=sub)) + geom_violin(alpha=.7, linewidth=.3) +
  geom_boxplot(width=.18, outlier.size=.2, linewidth=.3) + scale_fill_manual(values=c(PAL$low,PAL$high), guide="none") +
  labs(title="Disulfidptosis activity", x=NULL, y="Mean z-score")

fig <- add_tags((pA | pB) / (pC | pD) + plot_layout(heights=c(1.15,1)))
save_pub(fig, "results/Fig3_subtype_v2", w_mm=183, h_mm=122)
cat("Fig2 升级版完成\n")
