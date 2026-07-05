## Fig 8 候选基因: 上调DEG ∩ 二硫死亡 ∩ WGCNA模块 取交集 + 候选基因生存
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(survival); library(survminer); library(patchwork) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
deg <- fread("data/raw/kirc_DEG.csv"); setnames(deg,1,"gene")
up  <- deg$gene[deg$adj.P.Val<0.05 & deg$logFC>1]
wg  <- readRDS("results/kirc_wgcna.rds")$genes
sets <- list(`Up-DEG`=up, `WGCNA-ME9`=wg, `Disulfidptosis`=dis)
cand <- Reduce(intersect, list(up, wg))                 # 上调 + 二硫死亡共表达模块
cat("候选基因(上调DEG ∩ ME9模块):", length(cand), "\n")

## Venn via VennDiagram -> grob
suppressPackageStartupMessages({ library(VennDiagram); library(grid) })
futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger")
vg <- venn.diagram(sets, filename=NULL, fill=c("#E64B35","#4DBBD5","#00A087"), alpha=.5,
                   cex=1.1, cat.cex=1, margin=.1, main="A  Candidate gene intersection")
pV <- wrap_elements(grid::grobTree(vg))

## 候选基因 univariate Cox -> 取最显著做 KM
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
cg <- intersect(cand, rownames(expr)); X <- t(expr[cg, df$sample]); colnames(X)<-cg
cp <- sapply(cg, function(gn) tryCatch(summary(coxph(Surv(df$OS.time,df$OS)~X[,gn]))$coefficients[1,5], error=function(e)NA))
best <- names(sort(cp))[1]
df$gg <- ifelse(X[,best] > median(X[,best]), "High", "Low")
pK <- ggsurvplot(survfit(Surv(OS.time/365,OS)~gg,df), data=df, pval=TRUE, palette=c("#E64B35","#3C5488"),
  legend.title=best, xlab="Years", title=sprintf("B  %s (top candidate)", best))$plot

fig <- pV | pK
ggsave("results/Fig10_candidate.pdf", fig, width=13, height=6)
ggsave("results/Fig10_candidate.png", fig, width=13, height=6, dpi=130)
write.csv(data.frame(gene=cg, cox_p=cp[cg]), "results/kirc_candidates.csv", row.names=FALSE)
cat(sprintf("Fig8 完成 | 候选 %d 个, 最显著 %s (Cox p=%.1e)\n", length(cg), best, min(cp,na.rm=TRUE)))
