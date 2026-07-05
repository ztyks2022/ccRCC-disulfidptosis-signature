## 补充图: 高 vs 低风险组 GSEA 通路差异(离线 gseGO)
suppressPackageStartupMessages({ library(data.table); library(limma); library(clusterProfiler)
  library(org.Hs.eg.db); library(ggplot2) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
tum <- expr[, df$sample]; grp <- factor(df$rgrp, c("Low","High"))
fit <- eBayes(lmFit(tum, model.matrix(~ grp)))
res <- topTable(fit, coef = 2, number = Inf)
rl  <- sort(setNames(res$logFC, rownames(res)), decreasing = TRUE)

gse <- gseGO(rl, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", pvalueCutoff = 0.05, verbose = FALSE)
gd  <- as.data.frame(gse)
cat("GSEA 显著通路:", nrow(gd), "\n")
top <- rbind(head(gd[order(-gd$NES),], 10), head(gd[order(gd$NES),], 10))
top$Description <- factor(top$Description, levels = top$Description[order(top$NES)])
top$dir <- ifelse(top$NES > 0, "High-risk up", "High-risk down")
p <- ggplot(top, aes(NES, Description, fill = dir)) + geom_col() +
  scale_fill_manual(values = c("High-risk up" = "#E64B35", "High-risk down" = "#3C5488")) +
  theme_bw() + labs(title = "GSEA: high vs low risk (GO-BP)", x = "Normalized enrichment score", y = NULL) +
  theme(axis.text.y = element_text(size = 8), legend.title = element_blank())
ggsave("results/SuppFig_GSEA.pdf", p, width = 10, height = 7)
ggsave("results/SuppFig_GSEA.png", p, width = 10, height = 7, dpi = 130)
saveRDS(gse, "results/kirc_gsea.rds")
cat("SuppFig GSEA 完成 | 高危上调 top:", paste(head(top$Description[top$NES>0],3), collapse=" | "), "\n")
