## 外部验证: 把 12 基因 LASSO 签名搬到独立 ccRCC 队列 (GSE29609, Affy, n~39, 有生存)
## 跨平台迁移标准做法: 验证集内每基因 z-score 标准化后乘以 TCGA 训练得到的 LASSO 系数
suppressPackageStartupMessages({ library(GEOquery); library(survival) })
dir.create("results", showWarnings = FALSE)
COEF <- c(ACTB=0.508, CD2AP=0.252, SLC7A11=0.224, GYS1=0.167, NCKAP1=0.017,
          LRPPRC=-0.418, CYFIP1=-0.386, TLN1=-0.099, MYH10=-0.082,
          NUBPL=-0.043, NDUFS1=-0.012, WASF2=-0.009)

ok <- tryCatch({
  gse <- getGEO("GSE29609", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
  TRUE
}, error = function(e) { message("GEO 下载失败: ", conditionMessage(e)); FALSE })
if (!ok) quit(save = "no")

ex <- exprs(gse); pd <- pData(gse); fd <- fData(gse)
cat("队列维度: ", nrow(ex), "探针 x", ncol(ex), "样本\n")
cat("=== pData 列 ===\n"); print(colnames(pd))
chr <- grep("characteristics", colnames(pd), value = TRUE)
cat("=== characteristics 示例(首样本) ===\n"); print(unlist(pd[1, chr]))

## ---- 探针 -> 基因符号 (取平均表达最高的探针代表该基因) ----
sym_col <- grep("symbol|Symbol|GENE_SYMBOL|Gene Symbol", colnames(fd), value = TRUE)[1]
if (is.na(sym_col)) { message("无基因符号列, 终止"); quit(save="no") }
sym <- fd[[sym_col]]; sym <- sub(" ///.*", "", sym)        # 多映射取第一个
rs <- rowMeans(ex, na.rm = TRUE)
ord <- order(rs, decreasing = TRUE)
ex2 <- ex[ord, ]; sym2 <- sym[ord]
keep <- !duplicated(sym2) & sym2 != "" & !is.na(sym2)
mat <- ex2[keep, ]; rownames(mat) <- sym2[keep]

g <- intersect(names(COEF), rownames(mat))
cat("命中签名基因: ", length(g), "/12:", paste(g, collapse=","), "\n")
if (length(g) < 6) { message("命中基因过少, 验证不可靠, 仍输出"); }

## ---- z-score + 风险分 ----
z <- t(scale(t(mat[g, , drop = FALSE])))
risk <- as.numeric(COEF[g] %*% z)
names(risk) <- colnames(mat)

## ---- 解析生存 (OS 时间 + 状态), 尽量稳健地从 characteristics 抓 ----
allc <- apply(pd[, chr, drop = FALSE], 1, paste, collapse = " | ")
get_num <- function(pat) suppressWarnings(as.numeric(sub(".*?([0-9.]+).*", "\\1",
                          regmatches(allc, regexpr(pat, allc, ignore.case = TRUE)))))
time <- get_num("(survival|follow[- ]?up|os)[^0-9]*[0-9.]+")
ev_raw <- tolower(regmatches(allc, regexpr("status[^|]*", allc, ignore.case = TRUE)))
event <- ifelse(grepl("dead|decease|1|yes", ev_raw), 1, ifelse(grepl("alive|0|no", ev_raw), 0, NA))
cat("生存解析: 非缺失 time =", sum(!is.na(time)), " event =", sum(!is.na(event)), "\n")

res <- data.frame(sample = names(risk), risk = risk, time = time, event = event)
write.csv(res, "results/kirc_extval_GSE29609.csv", row.names = FALSE)

ev_ok <- res[complete.cases(res[, c("risk","time","event")]) & res$time > 0, ]
if (nrow(ev_ok) >= 10) {
  ev_ok$grp <- ifelse(ev_ok$risk > median(ev_ok$risk), "High", "Low")
  sd <- survdiff(Surv(time, event) ~ grp, ev_ok)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  cc <- summary(coxph(Surv(time, event) ~ risk, ev_ok))$concordance[1]
  cat(sprintf("外部验证 GSE29609: n=%d, KM log-rank p=%.4g, C-index=%.3f\n", nrow(ev_ok), p, cc))
  saveRDS(list(res = ev_ok, p = p, cindex = cc), "results/kirc_extval_GSE29609.rds")
} else {
  cat("可用样本不足(", nrow(ev_ok), "), 仅保存风险分; 生存字段需人工核对 characteristics\n")
}
