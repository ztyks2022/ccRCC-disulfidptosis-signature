## KIRC 二硫死亡 分子分型 + 免疫微环境
suppressPackageStartupMessages({ library(data.table); library(survival) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4",
         "MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1",
         "GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[, -1]); rownames(expr) <- g
tum <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
du  <- intersect(dis, rownames(tum)); sv <- fread("data/raw/KIRC_survival.txt")

## (1) 分子分型:基于二硫死亡基因的一致性聚类(k=2)
d <- tum[du, ]; d <- t(scale(t(d)))
if (requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  cc <- ConsensusClusterPlus::ConsensusClusterPlus(d, maxK = 4, reps = 50, pItem = 0.8,
        clusterAlg = "km", distance = "euclidean", seed = 1, plot = NULL)
  sub <- cc[[2]]$consensusClass
} else { set.seed(1); sub <- kmeans(t(d), 2)$cluster }
sub <- factor(paste0("C", sub))

## 分型 -> 生存
df <- merge(data.frame(sample = colnames(tum), sub = sub, dscore = colMeans(d)),
            sv[, .(sample, OS, OS.time)], by = "sample")
df <- df[!is.na(df$OS) & df$OS.time > 0, ]
p_sub <- 1 - pchisq(survdiff(Surv(OS.time, OS) ~ sub, df)$chisq, length(unique(df$sub)) - 1)

## (2) 免疫微环境:免疫/检查点 signature ssGSEA,比分型
suppressPackageStartupMessages(library(GSVA))
imm_sets <- list(
  CD8_Tcell = c("CD8A","CD8B","GZMA","GZMB","PRF1","IFNG"),
  Treg      = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Checkpoint= c("CD274","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","IDO1"),
  Cytolytic = c("GZMA","PRF1","GZMK","NKG7"))
imm_sets <- lapply(imm_sets, intersect, rownames(tum))
isc <- tryCatch({ gp <- ssgseaParam(tum, imm_sets); gsva(gp) },
                error = function(e) gsva(tum, imm_sets, method = "ssgsea"))
isc <- as.data.frame(t(isc)); isc$sample <- rownames(isc)
mg  <- merge(df, isc, by = "sample")

cat("\n================ KIRC 二硫死亡 分型 + 免疫 ================\n")
cat(sprintf("一致性聚类得 %d 型: %s\n", nlevels(sub), paste(levels(sub), table(df$sub), collapse=" ")))
cat(sprintf("分型 -> 生存  KM log-rank p = %.2e\n", p_sub))
cat("各亚型免疫评分(均值) + 亚型间差异(wilcox p):\n")
for (s in names(imm_sets)) {
  pv <- tryCatch(wilcox.test(mg[[s]] ~ mg$sub)$p.value, error = function(e) NA)
  cat(sprintf("  %-11s C1=%.3f C2=%.3f  p=%.1e\n", s,
      mean(mg[[s]][mg$sub=="C1"]), mean(mg[[s]][mg$sub=="C2"]), pv))
}
saveRDS(mg, "data/raw/kirc_subtype_immune.rds")
