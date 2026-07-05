## MR: 二硫死亡签名基因表达(blood cis-eQTL, GTEx QTD000356)→ ccRCC 风险(FinnGen R12)
## 全程远程 tabix,不下大文件。单 lead cis-SNP Wald ratio(cis-MR/SMR 式)。
suppressPackageStartupMessages({ library(Rsamtools); library(data.table) })
source("code/kirc/_style.R")
EQTL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz"
FIN  <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_KIDNEY_CLEAR_CELL_CARCINOMA_EXALLC.gz"
PCUT <- 5e-6
## 远程 tabix 易掉线 → 每次开新连接 + 重试
q_tabix <- function(url, gr, tries=5){
  for(t in 1:tries){
    r <- tryCatch(scanTabix(TabixFile(url, index=paste0(url,".tbi")), param=gr)[[1]],
                  error=function(e) NULL)
    if(!is.null(r)) return(r)
    Sys.sleep(2)
  }
  character(0)
}

## 12 基因 hg38 坐标 + Ensembl(硬编码,稳)
GI <- data.table(
  gene = c("ACTB","CD2AP","SLC7A11","GYS1","NCKAP1","LRPPRC","CYFIP1","TLN1","MYH10","NUBPL","NDUFS1","WASF2"),
  ens  = c("ENSG00000075624","ENSG00000198087","ENSG00000151012","ENSG00000104812","ENSG00000061676","ENSG00000138095","ENSG00000273749","ENSG00000137076","ENSG00000133026","ENSG00000151413","ENSG00000023228","ENSG00000158195"),
  chr  = c("7","6","4","19","2","2","15","9","17","14","2","1"),
  start= c(5527147,47477745,138178486,48964861,182750621,43861766,22931684,35696604,8373418,31789030,206128479,27318990),
  end  = c(5530601,47627263,138233440,49002442,182890977,43968277,23003630,35732236,8520536,32134837,206159696,27372544))
ALT  <- list(CYFIP1 = c("ENSG00000273749","ENSG00000068793"))
strip <- function(x) sub("\\..*", "", x)

res <- list()
for(i in 1:nrow(GI)){
  gi <- GI[i]; ens_set <- unique(c(gi$ens, ALT[[gi$gene]]))
  reg <- GRanges(gi$chr, IRanges(max(1, gi$start-1e6), gi$end+1e6))
  ex <- q_tabix(EQTL, reg)
  if(!length(ex)){ cat(sprintf("%-8s eQTL 区无数据\n", gi$gene)); next }
  M <- fread(text=paste(ex, collapse="\n"), header=FALSE, sep="\t")
  sub <- M[strip(V1) %in% ens_set | strip(V17) %in% ens_set]
  sub <- sub[as.numeric(V9) < PCUT]
  if(!nrow(sub)){ cat(sprintf("%-8s 无显著 cis 工具(p<%.0e)\n", gi$gene, PCUT)); next }
  ld <- sub[which.min(as.numeric(V9))]
  chr<-as.character(ld$V2); pos<-as.integer(ld$V3); ref_e<-ld$V4; alt_e<-ld$V5
  be<-as.numeric(ld$V10); se_e<-as.numeric(ld$V11); rs<-ld$V19; pe<-as.numeric(ld$V9)
  pal <- (toupper(ref_e)=="A"&toupper(alt_e)=="T")|(toupper(ref_e)=="T"&toupper(alt_e)=="A")|(toupper(ref_e)=="C"&toupper(alt_e)=="G")|(toupper(ref_e)=="G"&toupper(alt_e)=="C")
  if(pal){ cat(sprintf("%-8s lead %s 回文SNP,跳过\n", gi$gene, rs)); next }
  fx <- q_tabix(FIN, GRanges(chr, IRanges(pos,pos)))
  if(!length(fx)){ cat(sprintf("%-8s FinnGen 无该位点 (%s:%d)\n", gi$gene, chr, pos)); next }
  F <- fread(text=paste(fx, collapse="\n"), header=FALSE, sep="\t")
  hit <- F[(V3==ref_e & V4==alt_e) | (V3==alt_e & V4==ref_e)]
  if(!nrow(hit)){ cat(sprintf("%-8s FinnGen 等位不匹配\n", gi$gene)); next }
  hit <- hit[1]; bo <- as.numeric(hit$V9); se_o <- as.numeric(hit$V10)
  if(hit$V4!=alt_e) bo <- -bo
  bmr <- bo/be; smr <- sqrt(se_o^2/be^2 + bo^2*se_e^2/be^4); z<-bmr/smr; pmr<-2*pnorm(-abs(z))
  res[[gi$gene]] <- data.table(gene=gi$gene, rsid=rs, p_eqtl=pe, beta_exp=round(be,3), beta_out=round(bo,3),
                         OR=exp(bmr), L=exp(bmr-1.96*smr), U=exp(bmr+1.96*smr), p_MR=pmr)
  cat(sprintf("%-8s ✓ %s | OR=%.2f (%.2f-%.2f) p_MR=%.3g\n", gi$gene, rs, exp(bmr), exp(bmr-1.96*smr), exp(bmr+1.96*smr), pmr))
}
R <- rbindlist(res)
if(nrow(R)){
  fwrite(R, "results/kirc_MR_results.csv"); saveRDS(R, "results/kirc_MR_results.rds")
  cat("\n==== MR 汇总 ====\n"); print(R)
} else cat("无可用工具变量,MR 未产出\n")
