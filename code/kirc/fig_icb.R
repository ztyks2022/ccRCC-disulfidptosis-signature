## Immunotherapy-response-relevant signatures by risk group (TCGA-KIRC)
## Validated ICB-response predictors scored by ssGSEA (local, reproducible; not the TIDE web tool).
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(GSVA); library(ggplot2); library(ggpubr) })

m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- as.data.table(P$df)
df$grp <- factor(ifelse(df$risk > median(df$risk), "High", "Low"), c("Low","High"))
cs <- intersect(df$sample, colnames(expr)); expr <- expr[, cs]
grp <- df$grp[match(cs, df$sample)]

## --- validated ICB-response signatures ---
sigs <- list(
  `T-cell-inflamed GEP` = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
                            "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
                            "PSMB10","STAT1","TIGIT"),
  `IFN-gamma`           = c("IFNG","STAT1","IDO1","CXCL9","CXCL10","CXCL11","HLA-DRA"),
  `T-cell exhaustion`   = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","VSIR"),
  `Cytolytic activity`  = c("GZMA","GZMB","PRF1","GNLY","NKG7")
)
sigs <- lapply(sigs, function(x) intersect(x, rownames(expr)))
par_gsva <- ssgseaParam(expr, sigs)
sc <- gsva(par_gsva)                                   # signatures x samples

long <- rbindlist(lapply(rownames(sc), function(s)
  data.table(sig=s, score=sc[s,], grp=grp)))
long$sig <- factor(long$sig, levels=names(sigs))

pv <- long[, .(p=wilcox.test(score~grp)$p.value), by=sig]
cat("ICB signatures High-vs-Low Wilcoxon p:\n"); print(pv)

pICB <- ggplot(long, aes(grp, score, fill=grp)) +
  geom_violin(alpha=.45, colour=NA) + geom_boxplot(width=.3, outlier.size=.4, alpha=.9) +
  facet_wrap(~sig, scales="free_y", nrow=1) +
  stat_compare_means(comparisons=list(c("Low","High")), method="wilcox.test",
                     label="p.signif", size=2.6, tip.length=.01) +
  scale_fill_manual(values=c(Low=PAL$low, High=PAL$high), guide="none") +
  labs(title="Immunotherapy-response signatures by risk group", x=NULL, y="ssGSEA score") +
  theme_nat() + theme(strip.text=element_text(size=6.5))

saveRDS(list(scores=sc, grp=grp, pv=pv), "results/kirc_icb.rds")
if (!exists("ASSEMBLE_MODE")) {
  save_pub(pICB, "results/Fig_ICB", w_mm=170, h_mm=68)
  cat("ICB signatures done\n")
}
