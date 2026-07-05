## Fig 3 预后核心图: A 单因素Cox森林 | D 风险KM | E time-ROC | C 风险三联图
suppressPackageStartupMessages({ library(data.table); library(survival); library(survminer)
  library(timeROC); library(ggplot2); library(patchwork); library(scales) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df; lg <- P$lasso
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9",
         "MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11",
         "NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
tum <- expr[, substr(colnames(expr),14,15)=="01"]; du <- intersect(dis, rownames(tum))
X <- t(tum[du, df$sample]); colnames(X) <- du
df$Risk <- factor(df$rgrp, c("Low","High")); yr <- df$OS.time/365

## A 单因素 Cox 森林
fr <- do.call(rbind, lapply(du, function(gn){ s <- summary(coxph(Surv(df$OS.time, df$OS) ~ X[,gn]))
  data.frame(gene=gn, HR=s$coefficients[1,2], lo=s$conf.int[1,3], hi=s$conf.int[1,4], p=s$coefficients[1,5]) }))
fr <- fr[fr$p < 0.05, ]; fr <- fr[order(fr$HR), ]; fr$gene <- factor(fr$gene, fr$gene)
pA <- ggplot(fr, aes(HR, gene)) + geom_errorbarh(aes(xmin=lo, xmax=hi), height=.25, color="grey40") +
  geom_point(aes(color=HR>1), size=2) + geom_vline(xintercept=1, lty=2) + scale_x_log10() +
  scale_color_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="A  Univariate Cox (disulfidptosis genes)", x="Hazard ratio", y=NULL) + theme(legend.position="none")

## D 风险 KM
pD <- ggsurvplot(survfit(Surv(yr, df$OS) ~ Risk, data=df), data=df, pval=TRUE,
  palette=c("#3C5488","#E64B35"), legend.title="Risk", legend.labs=c("Low","High"),
  xlab="Years", ylab="Overall survival", title="D  Risk model KM")$plot

## E time-ROC 1/3/5 年
tr <- timeROC(T=yr, delta=df$OS, marker=df$risk, cause=1, times=c(1,3,5), iid=FALSE)
rocdf <- do.call(rbind, lapply(1:3, function(i)
  data.frame(FPR=tr$FP[,i], TPR=tr$TP[,i], yr=sprintf("%d-yr  AUC=%.2f", c(1,3,5)[i], tr$AUC[i]))))
pE <- ggplot(rocdf, aes(FPR, TPR, color=yr)) + geom_line(linewidth=.8) +
  geom_abline(lty=2, color="grey") + scale_color_manual(values=c("#E64B35","#4DBBD5","#00A087")) +
  theme_bw() + labs(title="E  Time-dependent ROC", x="1-Specificity", y="Sensitivity") +
  theme(legend.position=c(.65,.22), legend.title=element_blank())

## C 风险三联图
ord <- order(df$risk); d2 <- df[ord, ]; d2$rank <- seq_len(nrow(d2))
pC1 <- ggplot(d2, aes(rank, risk, color=Risk)) + geom_point(size=.5) +
  scale_color_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="C  Risk score / status / expression", y="Risk score", x=NULL) + theme(legend.position="none")
pC2 <- ggplot(d2, aes(rank, yr[ord], color=factor(OS))) + geom_point(size=.5) +
  scale_color_manual(values=c("grey70","#E64B35"), labels=c("Alive","Dead")) + theme_bw() +
  labs(y="Years", x=NULL) + theme(legend.title=element_blank(), legend.position="right")
hm <- t(scale(X[ord, lg, drop=FALSE])); colnames(hm) <- seq_len(ncol(hm))
hmdf <- as.data.frame(as.table(as.matrix(hm))); hmdf$Var2 <- as.integer(as.character(hmdf$Var2))
pC3 <- ggplot(hmdf, aes(Var2, Var1, fill=Freq)) + geom_tile() +
  scale_fill_gradient2(low="#3C5488", mid="white", high="#E64B35", limits=c(-2,2), oob=squish, name="z") +
  theme_minimal() + labs(x="Patients (ranked by risk)", y=NULL) +
  theme(axis.text.x=element_blank(), panel.grid=element_blank())

fig <- (pA | pD | pE) / (pC1 / pC2 / pC3) + plot_layout(heights=c(1, 1.1))
ggsave("results/Fig2_prognosis.pdf", fig, width=15, height=10)
ggsave("results/Fig2_prognosis.png", fig, width=15, height=10, dpi=140)
cat("Fig3 完成 -> results/Fig3_prognosis.pdf/png | time-ROC AUC:",
    sprintf("1yr=%.2f 3yr=%.2f 5yr=%.2f", tr$AUC[1], tr$AUC[2], tr$AUC[3]), "\n")
