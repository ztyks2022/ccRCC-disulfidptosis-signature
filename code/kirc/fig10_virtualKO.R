## Fig 10 虚拟敲除(scTenifoldKnk): 敲 ACTB(签名最大权重基因),看扰动网络
## ★关键: 基因数控制在 ~2000(HVG ∪ 二硫死亡基因),否则 pcNet 在 1.2万基因上要跑数小时
suppressPackageStartupMessages({ library(Seurat); library(scTenifoldKnk); library(dplyr) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
seu <- readRDS("results/kirc_scRNA.rds")
m <- as.matrix(GetAssayData(seu, layer = "counts"))
m <- m[rowSums(m > 0) >= 0.02 * ncol(m), ]
## 取高变基因 top2000 ∪ 二硫死亡基因,确保 ACTB 在内 → 网络可算
ln <- log1p(sweep(m, 2, pmax(colSums(m),1), "/") * 1e4)
v  <- apply(ln, 1, var)
hvg  <- names(sort(v, decreasing = TRUE))[1:2000]
keep <- intersect(unique(c(hvg, dis, "ACTB")), rownames(m))
m2 <- m[keep, ]
gko <- "ACTB"
cat(sprintf("虚拟敲除: %d 基因(原%d)× %d 细胞; 敲除 %s\n", nrow(m2), nrow(m), ncol(m2), gko))

ko <- scTenifoldKnk(countMatrix = m2, gKO = gko)
dr <- ko$diffRegulation %>% arrange(p.value)
write.csv(dr, sprintf("results/kirc_scTenifoldKnk_%s.csv", gko), row.names = FALSE)
cat("完成 | 敲", gko, "最扰动:", paste(head(dr$gene, 8), collapse = ","), "\n")

## GO 富集 top 扰动基因
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  eg <- bitr(head(dr$gene, 30), "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  go <- tryCatch(enrichGO(eg, org.Hs.eg.db, ont = "BP", readable = TRUE), error = function(e) NULL)
  if (!is.null(go) && nrow(go) > 0) saveRDS(go, "results/kirc_KO_GO.rds")
}
