## Fig 11 验证图(3联): A 内部留出测试 / B 外部 CPTAC RNA-seq(平台匹配,阳性) / C 外部 GSE29609 双色(平台不匹配,阴性)
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(survival); library(ggplot2); library(patchwork) })

km_panel <- function(time, event, grp, title, sub){
  grp <- factor(grp, levels=c("Low","High"))
  fit <- survfit(Surv(time, event) ~ grp); st <- summary(fit)
  df  <- data.frame(time=st$time, surv=st$surv, grp=sub("grp=","",as.character(st$strata)))
  df  <- rbind(data.frame(time=0, surv=1, grp=levels(grp)), df)
  sd  <- survdiff(Surv(time,event)~grp); p <- 1-pchisq(sd$chisq, length(sd$n)-1)
  plab <- ifelse(p<1e-4, sprintf("log-rank p = %.1e", p), sprintf("log-rank p = %.3f", p))
  ggplot(df, aes(time, surv, colour=grp)) + geom_step(linewidth=.7) +
    scale_colour_manual(values=risk_cols, name="Risk") +
    scale_y_continuous(limits=c(0,1), expand=expansion(mult=c(0,.02))) +
    labs(title=title, subtitle=sub, x="Time (months)", y="Overall survival") +
    annotate("text", x=0, y=.06, hjust=0, label=plab, size=2.4, family="Arial") +
    theme_nat()
}

## A. TCGA 留出测试集(seed 42)
set.seed(42)
DIS <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); gg<-m[[1]]; X<-as.matrix(m[,-1]); rownames(X)<-gg
tum <- X[, substr(colnames(X),14,15)=="01"]; colnames(tum)<-substr(colnames(tum),1,15)
P <- readRDS("data/raw/kirc_prognosis.rds")$df; P$sample<-substr(P$sample,1,15)
sur <- P[!duplicated(P$sample), c("sample","OS","OS.time")]
gi <- intersect(DIS, rownames(tum)); E <- t(tum[gi,,drop=FALSE])
d <- merge(data.frame(sample=rownames(E), E, check.names=FALSE), sur, by="sample"); d <- d[complete.cases(d)&d$OS.time>0,]
idx <- sample(nrow(d), round(.7*nrow(d))); tr<-d[idx,]; te<-d[-idx,]
cv <- cv.glmnet(as.matrix(tr[,gi]), Surv(tr$OS.time,tr$OS), family="cox", alpha=1, nfolds=10)
co <- as.matrix(coef(cv, s="lambda.min")); sel <- rownames(co)[co[,1]!=0]
rte <- as.numeric(as.matrix(te[,sel,drop=FALSE]) %*% co[sel,1]); ci<-summary(coxph(Surv(te$OS.time,te$OS)~rte))$concordance[1]
pA <- km_panel(te$OS.time/30.44, te$OS, ifelse(rte>median(rte),"High","Low"),
               "TCGA-KIRC held-out test", sprintf("internal | n = %d | C-index = %.3f", nrow(te), ci))

## B. CPTAC RNA-seq 外部(阳性)
cp <- readRDS("results/kirc_extval_CPTAC.rds")
pB <- km_panel(cp$d$time, cp$d$ev, as.character(cp$d$grp),
               "External CPTAC (RNA-seq)", sprintf("platform-matched | n = %d | C-index = %.3f", cp$n, cp$cindex))

## C. GSE29609 双色芯片外部(阴性,平台不匹配)
ev <- readRDS("results/kirc_extval_GSE29609.rds")$d
pC <- km_panel(ev$time, ev$event, as.character(ev$grp),
               "External GSE29609 (two-colour array)", sprintf("platform-mismatched | n = %d | C-index = 0.627", nrow(ev)))

fig <- add_tags(pA + pB + pC) + plot_layout(guides="collect") & theme(legend.position="bottom")
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig5_validation_v2", w_mm=183, h_mm=72)
cat("Fig11 三联验证图完成\n")
