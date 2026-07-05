## KIRC 肿瘤 vs 正常 DEG + 二硫死亡 DEG + 富集
suppressPackageStartupMessages({ library(data.table); library(limma) })
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4",
         "MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1",
         "GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[, -1]); rownames(expr) <- g
st  <- substr(colnames(expr), 14, 15); keep <- st %in% c("01","11")
expr <- expr[, keep]; grp <- factor(ifelse(st[keep]=="11","Normal","Tumor"), c("Normal","Tumor"))

fit <- eBayes(lmFit(expr, model.matrix(~ grp)))
res <- topTable(fit, coef = 2, number = Inf)
up   <- rownames(res)[res$adj.P.Val < 0.05 & res$logFC >  1]
down <- rownames(res)[res$adj.P.Val < 0.05 & res$logFC < -1]
cat(sprintf("DEG(肿瘤vs正常, FDR<0.05,|logFC|>1): 上调 %d, 下调 %d\n", length(up), length(down)))

dis_de <- res[intersect(dis, rownames(res)), c("logFC","adj.P.Val")]
dis_de <- dis_de[order(dis_de$logFC), ]
cat(sprintf("二硫死亡基因中差异显著(FDR<0.05): %d/%d\n",
            sum(dis_de$adj.P.Val < 0.05), nrow(dis_de)))
cat("  关键: SLC7A11 logFC=%.2f(FDR=%.1e)\n" |>
    sprintf(res["SLC7A11","logFC"], res["SLC7A11","adj.P.Val"]))

## 上调 DEG KEGG 富集
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  eg <- bitr(up, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID
  kk <- enrichKEGG(eg, organism = "hsa")
  if (!is.null(kk) && nrow(kk) > 0) {
    cat("\n上调 DEG KEGG Top5:\n")
    print(head(as.data.frame(kk)[, c("Description","p.adjust","Count")], 5))
  }
}
write.csv(res, "data/raw/kirc_DEG.csv")
cat("\n03 完成 -> data/raw/kirc_DEG.csv\n")
