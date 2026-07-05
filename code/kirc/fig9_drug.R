## Fig 9 药敏: GDSC 预测 IC50,高 vs 低风险组差异药物
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(tidyr); library(dplyr) })
dp <- fread("results/gdsc/calcPhenotype_Output/DrugPredictions.csv")
setnames(dp, 1, "sample")
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
dp <- dp[sample %in% df$sample]
dp$Risk <- factor(df$rgrp[match(dp$sample, df$sample)], c("Low","High"))
drugs <- setdiff(colnames(dp), c("sample","Risk"))

## 每药 高vs低 wilcox,取差异最大的 12 个
pv <- sapply(drugs, function(d) tryCatch(wilcox.test(dp[[d]] ~ dp$Risk)$p.value, error=function(e) NA))
top <- names(sort(pv))[1:12]
long <- pivot_longer(dp[, c("Risk", top), with=FALSE], -Risk, names_to="drug", values_to="IC50")
long$drug <- factor(long$drug, levels=top)
pl <- ggplot(long, aes(Risk, IC50, fill=Risk)) + geom_boxplot(outlier.size=.2) +
  facet_wrap(~drug, scales="free_y", nrow=3) + scale_fill_manual(values=c("#3C5488","#E64B35")) +
  theme_bw() + labs(title="GDSC predicted drug sensitivity: high vs low risk (top differential)",
                    x=NULL, y="Predicted IC50") + theme(legend.position="none", strip.text=element_text(size=7))
ggsave("results/Fig9_drug.pdf", pl, width=12, height=8)
ggsave("results/Fig9_drug.png", pl, width=12, height=8, dpi=130)
write.csv(data.frame(drug=names(pv), p=pv), "results/kirc_drug_pvals.csv", row.names=FALSE)
cat("Fig9 药敏完成 | 差异最大药物:", paste(head(top,5), collapse=", "), "\n")
