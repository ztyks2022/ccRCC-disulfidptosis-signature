## TMB + mutation landscape by risk group (TCGA-KIRC)
## Downloads masked somatic MAF via TCGAbiolinks, computes TMB, compares high/low risk, oncoplot.
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(TCGAbiolinks); library(maftools); library(data.table)
  library(ggplot2); library(survival) })

## --- risk groups from the locked signature ---
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- as.data.table(P$df)
df$grp <- ifelse(df$risk > median(df$risk), "High", "Low")
df$bcr <- substr(df$sample, 1, 12)                      # patient barcode

## --- download / cache MAF ---
maf_rds <- "data/raw/kirc_maf.rds"
if (!file.exists(maf_rds)) {
  q <- GDCquery(project="TCGA-KIRC", data.category="Simple Nucleotide Variation",
                data.type="Masked Somatic Mutation",
                workflow.type="Aliquot Ensemble Somatic Variant Merging and Masking")
  GDCdownload(q, directory="data/raw/GDCdata")
  maf <- GDCprepare(q, directory="data/raw/GDCdata")
  saveRDS(maf, maf_rds)
} else maf <- readRDS(maf_rds)

maf <- as.data.table(maf)
maf$bcr <- substr(maf$Tumor_Sample_Barcode, 1, 12)
maf <- maf[bcr %in% df$bcr]
grp_map <- setNames(df$grp, df$bcr)
maf$grp <- grp_map[maf$bcr]

## --- TMB per patient (mutations / 38 Mb exome) ---
tmb <- maf[, .(n=.N), by=bcr]
tmb$TMB <- tmb$n / 38
tmb$grp <- grp_map[tmb$bcr]
tmb <- tmb[!is.na(grp)]
wp <- wilcox.test(TMB ~ grp, tmb)$p.value
cat(sprintf("TMB: High median %.2f vs Low median %.2f | Wilcoxon p=%.3g\n",
            median(tmb[grp=="High"]$TMB), median(tmb[grp=="Low"]$TMB), wp))

pTMB <- ggplot(tmb, aes(grp, TMB, fill=grp)) +
  geom_violin(alpha=.5, colour=NA) + geom_boxplot(width=.28, outlier.size=.5, alpha=.9) +
  scale_fill_manual(values=c(High=PAL$high, Low=PAL$low), guide="none") +
  scale_y_log10() +
  labs(title=sprintf("TMB by risk (p=%.2g)", wp), x=NULL, y="TMB (mut/Mb)") +
  theme_nat()

## --- top mutated genes by group (for oncoplot via maftools) ---
maf_obj <- read.maf(maf[!is.na(grp)], clinicalData=data.frame(
  Tumor_Sample_Barcode=unique(maf[!is.na(grp)]$Tumor_Sample_Barcode),
  Risk=grp_map[substr(unique(maf[!is.na(grp)]$Tumor_Sample_Barcode),1,12)]))

saveRDS(list(tmb=tmb, wp=wp), "results/kirc_tmb.rds")

if (!exists("ASSEMBLE_MODE")) {
  save_pub(pTMB, "results/Fig_TMB", w_mm=80, h_mm=80)
  png("results/Fig_oncoplot.png", width=1400, height=1500, res=200)
  oncoplot(maf_obj, top=20, clinicalFeatures="Risk", sortByAnnotation=TRUE)
  dev.off()
  cat("TMB + oncoplot done\n")
}
