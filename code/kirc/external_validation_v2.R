## 外部验证 v2: GSE29609 (n=39, GPL1708 Agilent, ccRCC, 有 OS) 直接解析,不靠 GEOquery
## 把 TCGA 训练得到的 12 基因 LASSO 系数搬过去: 验证集内每基因 z-score -> 风险分 -> KM/C-index
suppressPackageStartupMessages({ library(survival) })
COEF <- c(ACTB=0.508, CD2AP=0.252, SLC7A11=0.224, GYS1=0.167, NCKAP1=0.017,
          LRPPRC=-0.418, CYFIP1=-0.386, TLN1=-0.099, MYH10=-0.082,
          NUBPL=-0.043, NDUFS1=-0.012, WASF2=-0.009)
F  <- "data/raw/geo/GSE29609_series_matrix.txt.gz"
AN <- "data/raw/geo/GPL1708.annot.gz"

## ---- 1. 读表达矩阵 (probe x sample) ----
L   <- readLines(gzfile(F));
b   <- grep("series_matrix_table_begin", L); e <- grep("series_matrix_table_end", L)
expr <- read.table(text = L[(b+1):(e-1)], header = TRUE, sep = "\t",
                   row.names = 1, quote = "\"", check.names = FALSE)
expr <- as.matrix(expr); cat("表达矩阵:", nrow(expr), "probe x", ncol(expr), "样本\n")

## ---- 2. 解析每样本生存 (characteristics 行是 ragged,逐样本扫) ----
acc   <- gsub('"', '', strsplit(L[grep("!Sample_geo_accession", L)[1]], "\t")[[1]][-1])
chl   <- L[grep("!Sample_characteristics_ch1", L)]
cells <- do.call(rbind, lapply(chl, function(x){ v <- strsplit(x, "\t")[[1]][-1]; gsub('^"|"$', '', v) }))
colnames(cells) <- acc
pull <- function(col, pat){ h <- grep(pat, col, ignore.case = TRUE, value = TRUE)
  if(length(h)==0) NA else suppressWarnings(as.numeric(sub(".*:\\s*", "", h[1]))) }
time  <- apply(cells, 2, pull, pat = "survival time")
death <- apply(cells, 2, pull, pat = "^death \\(1=yes")     # 全因死亡 = OS event
cat("生存: time 非缺失", sum(!is.na(time)), "/ death 非缺失", sum(!is.na(death)), "\n")

## ---- 3. GPL1708 probe -> gene symbol ----
A  <- readLines(gzfile(AN))
pb <- grep("platform_table_begin", A); pe <- grep("platform_table_end", A)
A2 <- if(length(pb) && length(pe)) A[(pb+1):(pe-1)] else A[grep("^ID\t", A)[1]:length(A)]
ann <- read.table(text = A2, header = TRUE, sep = "\t", quote = "",
                  comment.char = "", fill = TRUE, check.names = FALSE)
sym_col <- grep("^Gene symbol$", colnames(ann))[1]
map <- setNames(ann[[sym_col]], as.character(ann$ID))
sym <- map[rownames(expr)]; sym <- sub(" ?///.*", "", sym)   # 多映射取首个

## 每基因取平均表达最高的 probe 代表
rs  <- rowMeans(expr, na.rm = TRUE); ord <- order(rs, decreasing = TRUE)
expr2 <- expr[ord, ]; sym2 <- sym[ord]
keep  <- !is.na(sym2) & sym2 != "" & !duplicated(sym2)
mat   <- expr2[keep, ]; rownames(mat) <- sym2[keep]
cat("基因级矩阵:", nrow(mat), "基因\n")

g <- intersect(names(COEF), rownames(mat))
cat("命中签名基因", length(g), "/12:", paste(g, collapse = ","), "\n")
miss <- setdiff(names(COEF), g); if(length(miss)) cat("缺失:", paste(miss, collapse = ","), "\n")

## ---- 4. z-score + 风险分 ----
z    <- t(scale(t(mat[g, , drop = FALSE])))           # 每基因标准化
risk <- as.numeric(COEF[g] %*% z); names(risk) <- colnames(mat)

## ---- 5. 合表 + 生存分析 ----
d <- data.frame(sample = colnames(mat), risk = risk,
                time = time[colnames(mat)], event = death[colnames(mat)])
write.csv(d, "results/kirc_extval_GSE29609.csv", row.names = FALSE)
d <- d[complete.cases(d[, c("risk","time","event")]) & d$time > 0, ]
cat("可分析样本 n =", nrow(d), "\n")

if(nrow(d) >= 10){
  d$grp <- factor(ifelse(d$risk > median(d$risk), "High", "Low"), levels = c("Low","High"))
  sd  <- survdiff(Surv(time, event) ~ grp, d); p <- 1 - pchisq(sd$chisq, length(sd$n)-1)
  cox <- summary(coxph(Surv(time, event) ~ risk, d))
  ci  <- cox$concordance[1]; hr <- cox$coefficients[1, "exp(coef)"]; hp <- cox$coefficients[1, "Pr(>|z|)"]
  cat(sprintf(">>> 外部验证 GSE29609: n=%d | KM log-rank p=%.4g | C-index=%.3f | 连续风险 HR=%.2f (p=%.3g)\n",
              nrow(d), p, ci, hr, hp))
  saveRDS(list(d = d, p = p, cindex = ci, hr = hr, hp = hp, n = nrow(d), genes = g),
          "results/kirc_extval_GSE29609.rds")
} else cat("可用样本不足,仅存风险分\n")
