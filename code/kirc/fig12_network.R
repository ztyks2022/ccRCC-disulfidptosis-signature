## Fig 12 网络合成图(信息密度升级,学 Cell Metab Fig2C / Nat Commun Fig2B):
## A. 上调 DEG 的 GO-BP 富集"团块网络"(emapplot,替代柱状图)
## B. WGCNA 模块 ME9 共表达 hub 力导向网络(ggraph)
source("code/kirc/_style.R")
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(igraph); library(ggraph); library(ggplot2); library(patchwork)
})

## ---------- A. 富集网络 ----------
deg <- fread("data/raw/kirc_DEG.csv"); setnames(deg, 1, "gene")
deg$gene <- sub("\\|.*", "", deg$gene)
up <- deg[logFC > 1 & adj.P.Val < 0.05 & gene != "?" , gene]
eg <- suppressWarnings(bitr(up, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID)
cat("上调 DEG:", length(up), "→ entrez:", length(eg), "\n")
ego <- enrichGO(eg, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
ego <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego <- pairwise_termsim(ego)
cat("富集条目:", nrow(as.data.frame(ego)), "\n")
## 手搓 term-相似网络(完全控制标签大小,和 B 面板统一为 ggraph 风格)
df <- as.data.frame(ego); df <- head(df[order(df$p.adjust), ], 22)
S  <- ego@termsim; common <- intersect(df$Description, rownames(S))
S  <- S[common, common, drop = FALSE]; S[is.na(S)] <- 0; S <- pmax(S, t(S)); diag(S) <- 0
gA <- graph_from_adjacency_matrix(S * (S >= 0.2), mode = "undirected", weighted = TRUE)
idx <- match(V(gA)$name, df$Description)
V(gA)$count <- df$Count[idx]; V(gA)$padj <- -log10(df$p.adjust)[idx]
V(gA)$lab <- sapply(V(gA)$name, function(s) paste(strwrap(s, 20), collapse = "\n"))
set.seed(2)
pA <- ggraph(gA, layout = "fr") +
  geom_edge_link(aes(edge_width = weight), colour = PAL$n_light, alpha = .55) +
  geom_node_point(aes(size = count, fill = padj), shape = 21, colour = "grey30", stroke = .25) +
  geom_node_text(aes(label = lab), size = 1.7, family = "Arial", lineheight = .8,
                 repel = TRUE, max.overlaps = Inf, segment.size = .15) +
  scale_edge_width(range = c(.15, .8), guide = "none") +
  scale_fill_gradient(low = PAL$high_soft, high = PAL$high, name = "-log10 p.adj") +
  scale_size(range = c(2, 8), name = "Gene count") +
  labs(title = "Up-regulated DEG pathway network (GO-BP)") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.4),
        legend.box = "horizontal")

## ---------- B. WGCNA ME9 共表达 hub 网络 ----------
w  <- readRDS("results/kirc_wgcna.rds"); g9 <- sub("\\|.*", "", w$genes)
m  <- fread("data/raw/KIRC_HiSeqV2.gz"); gg <- m[[1]]; X <- as.matrix(m[,-1]); rownames(X) <- gg
tum <- X[, substr(colnames(X),14,15) == "01"]
g9 <- intersect(g9, rownames(tum)); cat("ME9 命中表达:", length(g9), "基因\n")
C  <- cor(t(tum[g9, ])); diag(C) <- 0
ig <- graph_from_adjacency_matrix(C * (abs(C) >= 0.45), mode = "undirected", weighted = TRUE, diag = FALSE)
ig <- delete_vertices(ig, degree(ig) == 0)
E(ig)$cor    <- E(ig)$weight            # 保留正负号给边着色
E(ig)$weight <- abs(E(ig)$weight)       # FR 布局要求正权重
V(ig)$deg <- degree(ig)
hub <- names(sort(degree(ig), decreasing = TRUE))[1:12]
V(ig)$lab <- ifelse(V(ig)$name %in% hub, V(ig)$name, "")
cat("网络:", vcount(ig), "节点,", ecount(ig), "边\n")
set.seed(1)
pB <- ggraph(ig, layout = "fr") +
  geom_edge_link(aes(edge_width = weight, edge_colour = cor > 0), alpha = .35) +
  geom_node_point(aes(size = deg, fill = deg), shape = 21, colour = "grey30", stroke = .25) +
  geom_node_text(aes(label = lab), size = 2, family = "Arial", repel = TRUE,
                 max.overlaps = Inf, segment.size = .2) +
  scale_edge_width(range = c(.15, .9), guide = "none") +
  scale_edge_colour_manual(values = c(`TRUE` = PAL$high_soft, `FALSE` = PAL$low_soft), guide = "none") +
  scale_fill_gradient(low = PAL$low_soft, high = PAL$high, name = "Degree") +
  scale_size(range = c(1, 6), guide = "none") +
  labs(title = "WGCNA module ME9 co-expression hubs") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold"),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.6))

fig <- add_tags(pA + pB)
save_pub(fig, "results/Fig9_network_v2", w_mm = 183, h_mm = 95)
cat("Fig12 网络图完成\n")
