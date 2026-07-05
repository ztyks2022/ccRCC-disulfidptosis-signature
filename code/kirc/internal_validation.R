## 内部验证 (同平台,稳健): TCGA-KIRC 70/30 拆分,训练集重训 LASSO-Cox,测试集评估
## 目的: 给出诚实的、非过拟合乐观值的判别力(对比"全队列表观 C-index 0.674")
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(survival) })
set.seed(42)
DIS <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4",
         "MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1",
         "GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; X <- as.matrix(m[,-1]); rownames(X) <- g
tum <- X[, substr(colnames(X),14,15)=="01"]; colnames(tum) <- substr(colnames(tum),1,15)
P <- readRDS("data/raw/kirc_prognosis.rds")$df; P$sample <- substr(P$sample,1,15)
sur <- P[!duplicated(P$sample), c("sample","OS","OS.time")]
gg  <- intersect(DIS, rownames(tum))
E   <- t(tum[gg, , drop=FALSE])                          # sample x gene
d   <- merge(data.frame(sample=rownames(E), E, check.names=FALSE), sur, by="sample")
d   <- d[complete.cases(d) & d$OS.time>0, ]
cat("可用肿瘤样本:", nrow(d), "| 基因:", length(gg), "\n")

n <- nrow(d); idx <- sample(n, round(.7*n)); tr <- d[idx,]; te <- d[-idx,]
Xtr <- as.matrix(tr[,gg]); Xte <- as.matrix(te[,gg])
ytr <- Surv(tr$OS.time, tr$OS)
cv  <- cv.glmnet(Xtr, ytr, family="cox", alpha=1, nfolds=10)
co  <- as.matrix(coef(cv, s="lambda.min")); sel <- rownames(co)[co[,1]!=0]
cat("训练集 LASSO 选出", length(sel), "基因:", paste(sel, collapse=","), "\n")

rsk_tr <- as.numeric(Xtr[,sel,drop=FALSE] %*% co[sel,1])
rsk_te <- as.numeric(Xte[,sel,drop=FALSE] %*% co[sel,1])
c_tr <- summary(coxph(Surv(tr$OS.time,tr$OS)~rsk_tr))$concordance[1]
cxte <- summary(coxph(Surv(te$OS.time,te$OS)~rsk_te))
c_te <- cxte$concordance[1]; hr_te <- cxte$coef[1,"exp(coef)"]
te$grp <- ifelse(rsk_te>median(rsk_te),"High","Low")
sd <- survdiff(Surv(OS.time,OS)~grp, te); p_te <- 1-pchisq(sd$chisq,1)
cat(sprintf(">>> 内部验证: 训练C-index=%.3f | 测试C-index=%.3f | 测试KM p=%.3g | 测试HR=%.2f (n_test=%d)\n",
            c_tr, c_te, p_te, hr_te, nrow(te)))
saveRDS(list(sel=sel, c_train=c_tr, c_test=c_te, p_test=p_te, hr_test=hr_te,
             n_train=nrow(tr), n_test=nrow(te)), "results/kirc_internal_val.rds")
