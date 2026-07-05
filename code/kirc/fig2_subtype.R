## Fig 2 分子分型: A 分型基因热图 | B 分型KM | C 分型免疫 | D 分型二硫死亡评分
suppressPackageStartupMessages({ library(data.table); library(survival); library(survminer)
  library(ggplot2); library(patchwork); library(scales); library(tidyr); library(dplyr) })
mg  <- readRDS("data/raw/kirc_subtype_immune.rds")   # sample,sub,dscore,OS,OS.time,CD8_Tcell,Treg,Checkpoint,Cytolytic
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9",
         "MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11",
         "NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
du <- intersect(dis, rownames(expr)); mg$sub <- factor(mg$sub)

## A 分型基因热图(+顶部分型注释条)
ord <- order(mg$sub); hm <- t(scale(t(expr[du, mg$sample[ord]]))); colnames(hm) <- seq_len(ncol(hm))
hmdf <- as.data.frame(as.table(as.matrix(hm))); hmdf$Var2 <- as.integer(as.character(hmdf$Var2))
anno <- data.frame(x=seq_len(nrow(mg)), sub=mg$sub[ord])
pAa <- ggplot(anno, aes(x,1,fill=sub)) + geom_tile() + scale_fill_manual(values=c("#3C5488","#E64B35"), name="Subtype") +
  theme_void() + theme(legend.position="top")
pAh <- ggplot(hmdf, aes(Var2, Var1, fill=Freq)) + geom_tile() +
  scale_fill_gradient2(low="#3C5488", mid="white", high="#E64B35", limits=c(-2,2), oob=squish, name="z") +
  theme_minimal() + labs(x="Patients (by subtype)", y=NULL) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=6), panel.grid=element_blank())
pA <- (pAa / pAh + plot_layout(heights=c(.06,1))) + plot_annotation(title="A  Disulfidptosis genes by subtype")

## B 分型 KM
pB <- ggsurvplot(survfit(Surv(OS.time/365, OS) ~ sub, data=mg), data=mg, pval=TRUE,
  palette=c("#3C5488","#E64B35"), legend.title="Subtype", xlab="Years", ylab="Overall survival",
  title="B  Subtype survival")$plot

## C 分型免疫
imm <- mg %>% select(sub, CD8_Tcell, Cytolytic, Checkpoint, Treg) %>% pivot_longer(-sub, names_to="set", values_to="val")
pv <- imm %>% group_by(set) %>% summarise(p=wilcox.test(val~sub)$p.value, .groups="drop")
imm <- left_join(imm, pv, by="set"); imm$set <- sprintf("%s (p=%.1e)", imm$set, imm$p)
pC <- ggplot(imm, aes(sub, val, fill=sub)) + geom_boxplot(outlier.size=.2) + facet_wrap(~set, nrow=1, scales="free_y") +
  scale_fill_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="C  Immune infiltration by subtype", x=NULL, y="ssGSEA") + theme(legend.position="none")

## D 分型二硫死亡评分
pD <- ggplot(mg, aes(sub, dscore, fill=sub)) + geom_violin(alpha=.6) + geom_boxplot(width=.2, outlier.size=.2) +
  scale_fill_manual(values=c("#3C5488","#E64B35")) + theme_bw() +
  labs(title="D  Disulfidptosis activity", x=NULL, y="Mean z-score") + theme(legend.position="none")

fig <- (wrap_elements(pA) | pB) / (pC | pD) + plot_layout(heights=c(1.1,1))
ggsave("results/Fig2_subtype.pdf", fig, width=14, height=11)
ggsave("results/Fig2_subtype.png", fig, width=14, height=11, dpi=130)
cat("Fig2 完成 -> results/Fig2_subtype.*\n")
