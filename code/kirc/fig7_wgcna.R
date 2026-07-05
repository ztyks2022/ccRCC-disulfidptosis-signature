## Fig 7 WGCNA: 共表达模块 - 性状(二硫死亡/风险/分期/分级)相关 + hub 基因
suppressPackageStartupMessages({ library(WGCNA); library(data.table); library(GSVA); library(ggplot2) })
options(stringsAsFactors=FALSE); enableWGCNAThreads(4)
dis <- c("SLC7A11","SLC3A2","RPN1","NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g<-m[[1]]; expr<-as.matrix(m[,-1]); rownames(expr)<-g
tum <- expr[, substr(colnames(expr),14,15)=="01"]
## top 5000 高变基因
mad <- apply(tum,1,mad); top <- names(sort(mad,decreasing=TRUE))[1:5000]
datExpr <- t(tum[top,])

## 性状:二硫死亡评分 / 风险评分 / 分期 / 分级
du <- intersect(dis,rownames(tum))
dscore <- tryCatch({gp<-ssgseaParam(tum,list(d=du)); gsva(gp)["d",]}, error=function(e) gsva(tum,list(d=du),method="ssgsea")["d",])
P<-readRDS("data/raw/kirc_prognosis.rds"); cl<-fread("data/raw/KIRC_clinical.tsv")
traits <- data.frame(row.names=rownames(datExpr))
traits$Disulfidptosis <- dscore[rownames(datExpr)]
traits$RiskScore <- P$df$risk[match(rownames(datExpr), P$df$sample)]
traits$Stage <- as.numeric(factor(gsub("Stage ","",cl$pathologic_stage[match(rownames(datExpr),cl$sampleID)]),levels=c("I","II","III","IV")))
traits$Grade <- suppressWarnings(as.numeric(gsub("G","",cl$neoplasm_histologic_grade[match(rownames(datExpr),cl$sampleID)])))

sft <- pickSoftThreshold(datExpr, powerVector=1:20, verbose=0)
pw <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
net <- blockwiseModules(datExpr, power=pw, TOMType="unsigned", minModuleSize=30,
                        mergeCutHeight=0.25, numericLabels=TRUE, maxBlockSize=6000, verbose=0)
MEs <- orderMEs(net$MEs)
mtc <- cor(MEs, traits, use="p"); mtp <- corPvalueStudent(mtc, nrow(datExpr))

## 模块-性状热图
md <- reshape2::melt(mtc); md$p <- reshape2::melt(mtp)$value
md$lab <- sprintf("%.2f\n(%.0e)", md$value, md$p)
pdf("results/Fig8_WGCNA.pdf", width=9, height=8); par(mar=c(6,9,3,3))
labeledHeatmap(Matrix=mtc, xLabels=colnames(traits), yLabels=rownames(mtc), ySymbols=rownames(mtc),
  colors=blueWhiteRed(50), textMatrix=sprintf("%.2f\n(%.1e)",mtc,mtp), setStdMargins=FALSE,
  cex.text=0.6, zlim=c(-1,1), main="Fig 7  WGCNA module-trait relationships")
dev.off()
png("results/Fig8_WGCNA.png", width=900, height=800, res=110); par(mar=c(6,9,3,3))
labeledHeatmap(Matrix=mtc, xLabels=colnames(traits), yLabels=rownames(mtc), ySymbols=rownames(mtc),
  colors=blueWhiteRed(50), textMatrix=sprintf("%.2f\n(%.1e)",mtc,mtp), setStdMargins=FALSE,
  cex.text=0.6, zlim=c(-1,1), main="Fig 7  WGCNA module-trait relationships"); dev.off()

## 与二硫死亡最相关模块的 hub 基因
key_mod <- rownames(mtc)[which.max(abs(mtc[,"Disulfidptosis"]))]
mod_col <- gsub("ME","",key_mod); genes_in <- top[net$colors==as.integer(mod_col)]
cat(sprintf("Fig7 完成 | 与二硫死亡最相关模块: %s (r=%.2f), 含 %d 基因\n", key_mod, max(abs(mtc[,"Disulfidptosis"])), length(genes_in)))
saveRDS(list(net=net, traits=traits, mtc=mtc, key_mod=key_mod, genes=genes_in), "results/kirc_wgcna.rds")
