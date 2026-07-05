## KIRC 二硫死亡 预后主线:评分KM + 单因素Cox + LASSO-Cox风险模型 + C-index
## 数据: data/raw/KIRC_HiSeqV2.gz(表达) + KIRC_survival.txt(生存)
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(survival); library(glmnet) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4",
         "MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1",
         "GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")

m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[, -1]); rownames(expr) <- g
tum <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]   # 原发肿瘤
du  <- intersect(dis, rownames(tum))

## 二硫死亡 ssGSEA 评分
sc <- tryCatch({ gp <- ssgseaParam(tum, list(d = du)); gsva(gp)["d", ] },
               error = function(e) gsva(tum, list(d = du), method = "ssgsea")["d", ])

## 配生存
sv <- fread("data/raw/KIRC_survival.txt")
df <- merge(data.frame(sample = colnames(tum), score = as.numeric(sc)),
            sv[, .(sample, OS, OS.time)], by = "sample")
df <- df[!is.na(df$OS) & !is.na(df$OS.time) & df$OS.time > 0, ]
cat(sprintf("肿瘤+生存样本: %d\n", nrow(df)))

## (1) 评分高/低 KM
df$grp <- ifelse(df$score > median(df$score), "High", "Low")
p_km <- 1 - pchisq(survdiff(Surv(OS.time, OS) ~ grp, df)$chisq, 1)

## (2) 单因素 Cox(每个二硫死亡基因)
X <- t(tum[du, df$sample]); colnames(X) <- du
cox_p <- sapply(du, function(gn) tryCatch(
  summary(coxph(Surv(df$OS.time, df$OS) ~ X[, gn]))$coefficients[1, 5], error = function(e) NA))
sig <- names(cox_p)[which(cox_p < 0.05)]
cat(sprintf("单因素 Cox 显著基因: %d/%d (%s)\n", length(sig), length(du), paste(head(sig,8),collapse=",")))

## (3) LASSO-Cox 风险模型
set.seed(1); xs <- X[, sig, drop = FALSE]; ysurv <- Surv(df$OS.time, df$OS)
cv <- cv.glmnet(xs, ysurv, family = "cox", nfolds = 10)
cf <- as.numeric(coef(cv, s = "lambda.min")); names(cf) <- colnames(xs)
lg <- names(cf)[cf != 0]
df$risk <- as.numeric(xs[, lg, drop = FALSE] %*% cf[lg])
df$rgrp <- ifelse(df$risk > median(df$risk), "High", "Low")
p_risk <- 1 - pchisq(survdiff(Surv(OS.time, OS) ~ rgrp, df)$chisq, 1)
cidx <- summary(coxph(Surv(OS.time, OS) ~ risk, df))$concordance[1]

cat("\n================ KIRC 二硫死亡 预后结果 ================\n")
cat(sprintf("二硫死亡评分 高vs低  KM log-rank p = %.2e\n", p_km))
cat(sprintf("LASSO-Cox 风险模型 (%d 基因: %s)\n", length(lg), paste(lg, collapse = ",")))
cat(sprintf("  风险高vs低  KM log-rank p = %.2e   |   C-index = %.3f\n", p_risk, cidx))
saveRDS(list(df = df, lasso = lg, coef = cf[lg]), "data/raw/kirc_prognosis.rds")
cat(ifelse(p_risk < 0.01, ">>> 真实 KIRC 上二硫死亡预后模型显著 —— 肿瘤题强信号兑现 ✓\n", ">>> 信号一般\n"))
