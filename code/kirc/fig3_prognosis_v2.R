## Fig 3 升级版(Nature 审美): 统一配色/theme/角标/矢量导出
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(survival); library(survminer); library(timeROC) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df; lg <- P$lasso
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
tum <- expr[, substr(colnames(expr),14,15)=="01"]; du<-intersect(dis,rownames(tum))
X <- t(tum[du, df$sample]); colnames(X)<-du
df$Risk <- factor(df$rgrp, c("Low","High")); df$yr <- df$OS.time/365; yr <- df$yr

## A 单因素 Cox 森林(中性轴线 + 红高蓝低,直接标 HR)
fr <- do.call(rbind, lapply(du, function(gn){ s<-summary(coxph(Surv(df$OS.time,df$OS)~X[,gn]))
  data.frame(gene=gn, HR=s$coefficients[1,2], lo=s$conf.int[1,3], hi=s$conf.int[1,4], p=s$coefficients[1,5])}))
fr <- fr[fr$p<0.05,]; fr<-fr[order(fr$HR),]; fr$gene<-factor(fr$gene, fr$gene)
pA <- ggplot(fr, aes(HR, gene)) +
  geom_vline(xintercept=1, linewidth=.3, colour=PAL$n_mid, linetype=2) +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=0, linewidth=.45, colour=PAL$n_light) +
  geom_point(aes(colour=HR>1), size=1.5) +
  scale_colour_manual(values=c(`TRUE`=PAL$high,`FALSE`=PAL$low), guide="none") +
  scale_x_log10() + labs(title="Univariate Cox (disulfidptosis genes)", x="Hazard ratio (95% CI)", y=NULL)

## B 风险 KM(统一红蓝)
gg <- ggsurvplot(survfit(Surv(yr, OS)~Risk, data=df), data=df, palette=unname(risk_cols),
  pval=TRUE, pval.size=2.4, size=.55, censor.size=1.6, legend.title="Risk", legend.labs=c("Low","High"))
pB <- gg$plot + theme_nat() + labs(title="Risk-model survival", x="Time (years)", y="Overall survival") +
  theme(legend.position=c(.78,.86)) + scale_y_continuous(expand=expansion(mult=c(.02,.04)))

## C time-ROC
tr <- timeROC(T=yr, delta=df$OS, marker=df$risk, cause=1, times=c(1,3,5), iid=FALSE)
rocdf <- do.call(rbind, lapply(1:3, function(i) data.frame(FPR=tr$FP[,i], TPR=tr$TP[,i],
  yr=sprintf("%d-yr  AUC=%.2f", c(1,3,5)[i], tr$AUC[i]))))
pC <- ggplot(rocdf, aes(FPR,TPR,colour=yr)) +
  geom_abline(linewidth=.3, linetype=2, colour=PAL$n_light) + geom_line(linewidth=.7) +
  scale_colour_manual(values=c(PAL$high, PAL$accent3, PAL$low), name=NULL) +
  labs(title="Time-dependent ROC", x="1 - specificity", y="Sensitivity") +
  coord_equal() + theme(legend.position=c(.66,.18))

## D 风险三联图(评分 / 状态 / 基因热图,统一红蓝发散)
ord <- order(df$risk); d2 <- df[ord,]; d2$rank<-seq_len(nrow(d2)); d2$Risk<-factor(d2$rgrp,c("Low","High"))
pD1 <- ggplot(d2, aes(rank, risk, colour=Risk)) + geom_point(size=.35) +
  scale_colour_manual(values=risk_cols, guide="none") +
  labs(title="Risk score / survival status / signature expression", y="Risk score", x=NULL) +
  theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank())
pD2 <- ggplot(d2, aes(rank, yr[ord], colour=factor(OS))) + geom_point(size=.35) +
  scale_colour_manual(values=c(`0`=PAL$n_mid, `1`=PAL$high), labels=c("Alive","Dead"), name=NULL) +
  labs(y="Years", x=NULL) + theme(axis.text.x=element_blank(), axis.line.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position=c(.06,.8), legend.direction="horizontal")
hm <- t(scale(X[ord, lg, drop=FALSE])); colnames(hm)<-seq_len(ncol(hm))
hmdf <- as.data.frame(as.table(as.matrix(hm))); hmdf$Var2<-as.integer(as.character(hmdf$Var2))
pD3 <- ggplot(hmdf, aes(Var2, Var1, fill=Freq)) + geom_raster() +
  scale_fill_gradient2(low=PAL$low, mid="white", high=PAL$high, limits=c(-2,2), oob=squish, name="z-score") +
  labs(x="Patients (ranked by risk score)", y=NULL) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line=element_blank(),
        axis.text.y=element_text(size=5.2), legend.key.height=unit(3,"mm"))

fig <- (pA | pB | pC) / (pD1 / pD2 / pD3 + plot_layout(heights=c(1,1,1.5))) +
  plot_layout(heights=c(1, 1.25))
fig <- add_tags(fig)
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig2_prognosis_v2", w_mm=183, h_mm=125)
cat("Fig3 升级版完成 -> results/Fig3_prognosis_v2.{svg,pdf,png}\n")
