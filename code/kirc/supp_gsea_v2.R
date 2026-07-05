## SuppFig GSEA 升级版: 高 vs 低风险 GO-BP,统一样式(棒棒糖,语义红蓝)
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(clusterProfiler); library(ggplot2) })
gse <- readRDS("results/kirc_gsea.rds")
gd  <- as.data.frame(gse)
top <- rbind(head(gd[order(-gd$NES), ], 12), head(gd[order(gd$NES), ], 12))
top$Description <- factor(top$Description, levels = top$Description[order(top$NES)])
top$dir <- ifelse(top$NES > 0, "High-risk up", "High-risk down")
p <- ggplot(top, aes(NES, Description, colour = dir)) +
  geom_segment(aes(xend = 0, yend = Description), linewidth = .5, colour = PAL$n_light) +
  geom_point(aes(size = -log10(p.adjust))) +
  scale_colour_manual(values = c("High-risk up" = PAL$high, "High-risk down" = PAL$low), name = NULL) +
  scale_size(range = c(1.5, 4), name = "-log10 p.adj") +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey60") +
  labs(title = "GSEA: high- vs low-risk (GO biological process)",
       x = "Normalized enrichment score", y = NULL) +
  theme_nat() + theme(axis.text.y = element_text(size = 6), legend.position = "right")
save_pub(p, "results/SuppFig_GSEA_v2", w_mm = 160, h_mm = 120)
cat("SuppFig GSEA 升级版完成 | 高危上调 top:", paste(head(top$Description[top$NES > 0], 3), collapse = " | "), "\n")
