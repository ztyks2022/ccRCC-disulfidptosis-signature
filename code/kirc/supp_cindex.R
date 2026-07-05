## SuppFig S2: 预后判别力头对头 我们的签名 vs ClearCode34(TCGA 发现集 + CPTAC 外部)
source("code/kirc/_style.R")
suppressPackageStartupMessages(library(ggplot2))
df <- data.frame(
  cohort = factor(rep(c("TCGA-KIRC (discovery)","CPTAC (external)"), each=2),
                  levels=c("TCGA-KIRC (discovery)","CPTAC (external)")),
  method = rep(c("Disulfidptosis 12-gene signature","ClearCode34"), 2),
  cindex = c(0.678, 0.653, 0.633, 0.593))
df$method <- factor(df$method, levels=c("Disulfidptosis 12-gene signature","ClearCode34"))
p <- ggplot(df, aes(cohort, cindex, fill=method)) +
  geom_col(position=position_dodge(.72), width=.62, colour="grey30", linewidth=.25) +
  geom_text(aes(label=sprintf("%.3f", cindex)), position=position_dodge(.72), vjust=-.5, size=2.4, family="Arial") +
  geom_hline(yintercept=.5, linetype=2, colour=PAL$n_mid, linewidth=.3) +
  scale_fill_manual(values=c("Disulfidptosis 12-gene signature"=PAL$high, "ClearCode34"=PAL$n_mid), name=NULL) +
  scale_y_continuous(limits=c(0,.74), expand=expansion(mult=c(0,.04))) +
  labs(title="Prognostic concordance vs an established ccRCC classifier", x=NULL, y="C-index (overall survival)") +
  theme_nat() + theme(legend.position="top")
if (!exists("ASSEMBLE_MODE")) save_pub(p, "results/SuppFig2_cindex_compare_v2", w_mm=120, h_mm=85)
cat("SuppFig S2 完成\n")
