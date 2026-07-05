## 外部验证(平台匹配 RNA-seq): CPTAC ccRCC 表达(cBioPortal)+ 完整OS(GDC)→ 12基因签名
suppressPackageStartupMessages({ library(data.table); library(survival) })
COEF <- c(ACTB=0.508, CD2AP=0.252, SLC7A11=0.224, GYS1=0.167, NCKAP1=0.017,
          LRPPRC=-0.418, CYFIP1=-0.386, TLN1=-0.099, MYH10=-0.082,
          NUBPL=-0.043, NDUFS1=-0.012, WASF2=-0.009)
d <- fread("data/raw/extval/cptac_ccrcc_os.csv")
g <- names(COEF)
X <- log2(as.matrix(d[, ..g]) + 1)
Z <- scale(X)
d$risk <- as.numeric(Z[, g] %*% COEF[g])
d$time <- suppressWarnings(as.numeric(d$time_days)) / 30.44   # 天->月
d$ev   <- suppressWarnings(as.integer(d$event))
dd <- d[complete.cases(d[, .(risk, time, ev)]) & time > 0]
cat("可分析 n =", nrow(dd), "| 死亡:", sum(dd$ev), "| 删失:", sum(dd$ev==0), "\n")

dd$grp <- factor(ifelse(dd$risk > median(dd$risk), "High", "Low"), c("Low","High"))
sd <- survdiff(Surv(time, ev) ~ grp, dd); p <- 1 - pchisq(sd$chisq, 1)
cx <- summary(coxph(Surv(time, ev) ~ risk, dd))
ci <- cx$concordance[1]; hr <- cx$coef[1,"exp(coef)"]; hp <- cx$coef[1,"Pr(>|z|)"]
cat(sprintf(">>> CPTAC 外部验证(RNA-seq, n=%d): KM p=%.4g | C-index=%.3f | 连续HR=%.2f (p=%.3g)\n",
            nrow(dd), p, ci, hr, hp))
saveRDS(list(d=dd, p=p, cindex=ci, hr=hr, hp=hp, n=nrow(dd)), "results/kirc_extval_CPTAC.rds")
fwrite(dd[, .(sample, risk, time, ev, grp)], "results/kirc_extval_CPTAC.csv")
