## Fig 8 升级版: 候选基因 Venn + KM,统一样式
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(survival); library(survminer); library(VennDiagram); library(grid) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
deg <- fread("data/raw/kirc_DEG.csv"); setnames(deg,1,"gene")
up  <- deg$gene[deg$adj.P.Val<.05 & deg$logFC>1]; wg <- readRDS("results/kirc_wgcna.rds")$genes
sets <- list(`Up-DEG`=up, `WGCNA-ME9`=wg, `Disulfidptosis`=dis)
cand <- Reduce(intersect, list(up, wg))
futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger")
#vg <- venn.diagram(sets, filename=NULL, fill=c(PAL$high, PAL$low, PAL$accent), alpha=.45,
#  col=c(PAL$high,PAL$low,PAL$accent), lwd=1.2, cex=.7, cat.cex=.7, margin=.06,
#  cat.dist=c(.055,.055,.045), cat.pos=c(-25,25,180),
#  fontfamily="Arial", cat.fontfamily="Arial")
vg <- venn.diagram(
  sets,
  filename=NULL,
  fill=c(PAL$high, PAL$low, PAL$accent),
  alpha=.45,
  col=c(PAL$high,PAL$low,PAL$accent),
  lwd=1.2,
  cex=.68,                 # 里面数字大小
  cat.cex=.68,             # 外面标签大小
  margin=.035,
  cat.dist=c(.075,.045,.04),
  cat.pos=c(8,22,172),     # 第一个控制 Up-DEG
  fontfamily="Arial",
  cat.fontfamily="Arial"
)
# 整个 Venn(含外置类别标签)尽量填满面板,同时保留一点安全边距避免长标签被裁掉
pV <- wrap_elements(
  full = grid::grobTree(vg, vp=grid::viewport(width=.96, height=.92)),
  clip = FALSE
) +
  ggtitle("Candidate gene intersection") +
  theme(plot.title=element_text(size=7,face="bold"), plot.margin=margin(1, 1, 1, 1))
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
cg <- intersect(cand, rownames(expr)); X <- t(expr[cg, df$sample]); colnames(X)<-cg
cp <- sapply(cg, function(gn) tryCatch(summary(coxph(Surv(df$OS.time,df$OS)~X[,gn]))$coefficients[1,5], error=function(e)NA))
best <- names(sort(cp))[1]; df$gg <- ifelse(X[,best]>median(X[,best]),"High","Low")
gg <- ggsurvplot(survfit(Surv(OS.time/365,OS)~gg, data=df), data=df, palette=c(PAL$high,PAL$low),
  pval=TRUE, pval.size=2.4, size=.55, censor.size=1.6, legend.title=best, legend.labs=c("High","Low"))
pK <- gg$plot + theme_nat() + labs(title=sprintf("%s (top candidate)", best), x="Time (years)", y="Overall survival") +
  theme(legend.position=c(.78,.86))
fig <- add_tags(pV | pK)
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig10_candidate_v2", w_mm=170, h_mm=82)
cat("Fig8 升级版完成\n")
