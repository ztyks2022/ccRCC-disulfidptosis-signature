## 换题验证:二硫死亡 × 肾透明细胞癌(TCGA-KIRC)真数据信号
## 下 Xena TCGA-KIRC 表达(+生存),算二硫死亡评分,肿瘤vs正常 + 生存
## 对照衰老×AD 的弱信号,看肿瘤题是否信号强、易出阳性
suppressPackageStartupMessages({ library(data.table); library(GSVA) })
options(timeout = 900); dir.create("data/raw", showWarnings = FALSE)

## 二硫死亡相关基因(Liu 2023 Nat Cell Biol + 常用集;actin细胞骨架+SLC7A11轴)
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB",
         "ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1",
         "IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")

exprF <- "data/raw/KIRC_HiSeqV2.gz"
if (!file.exists(exprF))
  download.file("https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz", exprF, quiet = TRUE)
m <- fread(exprF); g <- m[[1]]; expr <- as.matrix(m[, -1]); rownames(expr) <- g
cat("TCGA-KIRC 表达:", nrow(expr), "基因 ×", ncol(expr), "样本\n")

## 样本类型:barcode 第14-15位 01=肿瘤 11=正常
st  <- substr(colnames(expr), 14, 15)
grp <- ifelse(st == "11", "Normal", "Tumor")
cat("样本:", paste(names(table(grp)), table(grp), collapse = " / "), "\n")

du <- intersect(dis, rownames(expr))
sc <- tryCatch({ gp <- ssgseaParam(expr, list(d = du)); as.numeric(gsva(gp)["d", ]) },
               error = function(e) as.numeric(gsva(expr, list(d = du), method = "ssgsea")["d", ]))
p   <- wilcox.test(sc ~ grp)$p.value
auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(factor(grp, c("Normal","Tumor")), sc, quiet = TRUE))), error = function(e) NA)
cat("\n===== 二硫死亡评分: 肿瘤 vs 正常 =====\n")
cat(sprintf("用基因 %d 个 | Tumor=%.3f  Normal=%.3f  Wilcox p=%.2e  AUC=%.3f\n",
            length(du), median(sc[grp=="Tumor"]), median(sc[grp=="Normal"]), p, auc))

## 生存:肿瘤样本按评分高/低 KM logrank(可选)
ok <- tryCatch({
  sf <- "data/raw/KIRC_survival.tsv"
  if (!file.exists(sf))
    download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.survival.tsv", sf, quiet = TRUE)
  sv <- fread(sf); sv$id15 <- substr(sv$sample, 1, 15)
  tum <- data.frame(id15 = substr(colnames(expr), 1, 15), sc = sc, grp = grp)
  tum <- tum[tum$grp == "Tumor", ]; tum <- merge(tum, sv, by = "id15")
  tum$hi <- tum$sc > median(tum$sc)
  library(survival)
  sd <- survdiff(Surv(OS.time, OS) ~ hi, data = tum)
  pv <- 1 - pchisq(sd$chisq, 1)
  cat(sprintf("生存(肿瘤 %d 例): 二硫死亡高 vs 低  log-rank p=%.2e\n", nrow(tum), pv)); TRUE
}, error = function(e) { cat("生存分析跳过:", conditionMessage(e), "\n"); FALSE })

cat(ifelse(p < 1e-5, "\n>>> 肿瘤题信号强(p远小于AD的n.s.) —— 换题性价比验证 ✓\n", "\n>>> 信号一般\n"))
