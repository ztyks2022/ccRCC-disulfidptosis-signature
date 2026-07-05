## Fig 10 升级版: ACTB 虚拟敲除。A=最扰动基因  B=扰动基因 GO 富集(统一样式)
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
dr <- fread("results/kirc_scTenifoldKnk_ACTB.csv")
dr <- dr[gene != "ACTB"]                       # 去掉自身(必然最高,无信息)
top <- head(dr[order(p.value)], 15); top$gene <- factor(top$gene, levels=rev(top$gene))
pA <- ggplot(top, aes(-log10(p.value+1e-300), gene)) +
  geom_segment(aes(xend=0, yend=gene), colour=PAL$n_light, linewidth=.5) +
  geom_point(aes(size=distance), colour=PAL$high) +
  scale_size(range=c(1.5,4.5), name="Network\ndistance") +
  labs(title="Most perturbed genes after in-silico ACTB knockout", x="-log10 p", y=NULL) +
  theme_nat() + theme(legend.position="right")

go <- as.data.frame(readRDS("results/kirc_KO_GO.rds"))
gt <- head(go[order(go$p.adjust), ], 8)
gt$lab <- sapply(as.character(gt$Description), function(s) paste(strwrap(s, 34), collapse="\n"))
gt$lab <- factor(gt$lab, levels=rev(gt$lab))
pB <- ggplot(gt, aes(-log10(p.adjust), lab)) +
  geom_col(fill=PAL$accent2, width=.66) +
  labs(title="Enriched processes (perturbed genes)", x="-log10 p.adjust", y=NULL) +
  theme_nat() + theme(axis.text.y=element_text(size=5.4))

fig <- add_tags(pA | pB) + plot_layout(widths=c(1, 1.1))
if (!exists("ASSEMBLE_MODE")) save_pub(fig, "results/Fig12_virtualKO_v2", w_mm=183, h_mm=80)
cat("Fig10 升级版完成 | top 扰动:", paste(head(as.character(top$gene[order(-(-log10(top[order(p.value)]$p.value)))]),6), collapse=","), "\n")
