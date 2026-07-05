## Fig 9 升级版: GDSC 预测药敏 高 vs 低风险。A=差异药物方向棒棒糖  B=top6 箱线
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(tidyr); library(patchwork) })
dp <- fread("results/gdsc/calcPhenotype_Output/DrugPredictions.csv"); setnames(dp, 1, "sample")
P  <- readRDS("data/raw/kirc_prognosis.rds")$df
dp <- dp[sample %in% P$sample]
dp$Risk <- factor(P$rgrp[match(dp$sample, P$sample)], c("Low","High"))
drugs <- setdiff(colnames(dp), c("sample","Risk"))

## 每药: wilcox p + 方向(High 中位 - Low 中位; <0 = 高危更敏感)
st <- rbindlist(lapply(drugs, function(d){
  x <- dp[[d]]; p <- tryCatch(wilcox.test(x ~ dp$Risk)$p.value, error=function(e) NA)
  data.table(drug=d, p=p, dHL=median(x[dp$Risk=="High"],na.rm=T)-median(x[dp$Risk=="Low"],na.rm=T))
}))
st <- st[!is.na(p)]; st$fdr <- p.adjust(st$p,"BH")
st$dir <- ifelse(st$dHL<0,"High-risk more sensitive","Low-risk more sensitive")
st$lab <- gsub("_[0-9]+$","",st$drug)
top <- st[order(p)][1:15]; top$lab <- factor(top$lab, levels=rev(top$lab))

pA <- ggplot(top, aes(-log10(p), lab, colour=dir)) +
  geom_segment(aes(xend=0, yend=lab), linewidth=.5, colour=PAL$n_light) +
  geom_point(aes(size=abs(dHL))) +
  scale_colour_manual(values=c("High-risk more sensitive"=PAL$high,"Low-risk more sensitive"=PAL$low), name=NULL) +
  scale_size(range=c(1.5,4), name="|ΔIC50|") +
  labs(title="Differential predicted drug sensitivity (GDSC)", x="-log10 Wilcoxon p", y=NULL) +
  theme_nat() + theme(legend.position="right", legend.box="vertical")

top6 <- as.character(top$drug[1:6])
long <- melt(dp[,c("Risk",top6),with=FALSE], id.vars="Risk", variable.name="drug", value.name="IC50")
long$drug <- factor(gsub("_[0-9]+$","",long$drug), levels=gsub("_[0-9]+$","",top6))
pB <- ggplot(long, aes(Risk, IC50, fill=Risk)) +
  geom_boxplot(outlier.size=.15, linewidth=.3, width=.6) +
  facet_wrap(~drug, scales="free_y", nrow=2) +
  scale_fill_manual(values=risk_cols, guide="none") +
  labs(title="Top differential drugs", x=NULL, y="Predicted IC50 (log)") +
  theme_nat() + theme(strip.text=element_text(size=6.2))

fig <- add_tags(pA | pB) + plot_layout(widths=c(1,1.15))
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig11_drug_v2", w_mm=183, h_mm=95)
fwrite(st[order(p)], "results/kirc_drug_pvals.csv")
cat("Fig9 升级版完成 | 高危更敏感 top:", paste(head(st[dHL<0][order(p)]$lab,5),collapse=", "), "\n")
