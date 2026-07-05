## KIRC disulfidptosis: does the prognostic signal live in the CYTOSKELETAL subset or the METABOLIC subset?
## Fair head-to-head: out-of-fold (10-fold CV) C-index for cyto-only vs metab-only LASSO-Cox,
## a gene-count-matched control (cyto sub-sampled to 9), and per-gene univariable Cox.
## Same data prep / cohort as 01_disulfidptosis_prognosis.R. Honest: report whatever comes out.
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet) })
set.seed(1)

## a priori classification of the 26-gene disulfidptosis set
cyto  <- c("NCKAP1","WASF2","BRK1","ABI2","CYFIP1","ACTB","ACTN4","MYH9","MYH10","MYL6",
           "FLNA","FLNB","TLN1","INF2","CD2AP","PDLIM1","IQGAP1")            # actin / mechanotransduction
metab <- c("SLC7A11","SLC3A2","RPN1","GYS1","NDUFA11","NDUFS1","OXSM","LRPPRC","NUBPL") # transport/glycogen/mito

## data prep (mirror 01_*)
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[, -1]); rownames(expr) <- g
tum <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
sv  <- fread("data/raw/KIRC_survival.txt")
samp <- intersect(colnames(tum), sv$sample)
sv2 <- sv[match(samp, sv$sample), ]
keep <- !is.na(sv2$OS) & !is.na(sv2$OS.time) & sv2$OS.time > 0
samp <- samp[keep]; y <- Surv(sv2$OS.time[keep], sv2$OS[keep])
cyto  <- intersect(cyto,  rownames(tum)); metab <- intersect(metab, rownames(tum))
cat(sprintf("samples=%d | cyto genes=%d | metab genes=%d\n", length(samp), length(cyto), length(metab)))

Xall <- t(tum[, samp])           # samples x genes (log2)

## out-of-fold CV C-index: fit LASSO-Cox per fold, score the held-out fold, take C WITHIN each fold, average (event-weighted)
oofC <- function(genes, folds=10, seed=1){
  genes <- intersect(genes, colnames(Xall)); if(length(genes) < 2) return(NA)
  X <- Xall[, genes, drop=FALSE]; n <- nrow(X)
  set.seed(seed); fold <- sample(rep(1:folds, length.out=n))
  cs <- w <- numeric(0)
  for(k in 1:folds){
    tr <- fold!=k; te <- fold==k
    cv <- tryCatch(cv.glmnet(X[tr,,drop=FALSE], y[tr], family="cox", nfolds=10), error=function(e) NULL)
    if(is.null(cv)) next
    lp <- as.numeric(predict(cv, newx=X[te,,drop=FALSE], s="lambda.min", type="link"))
    if(sd(lp)==0) next
    cc <- tryCatch(as.numeric(summary(coxph(y[te] ~ lp))$concordance[1]), error=function(e) NA)
    ev <- sum(y[te][,2])                       # events in this fold
    if(is.finite(cc) && ev>0){ cs <- c(cs, cc); w <- c(w, ev) }
  }
  if(!length(cs)) return(NA)
  sum(cs*w)/sum(w)
}
## apparent (in-sample) C-index, same reporting style as the manuscript's 0.68
appC <- function(genes){
  genes <- intersect(genes, colnames(Xall)); if(length(genes) < 2) return(NA)
  X <- Xall[, genes, drop=FALSE]
  cv <- cv.glmnet(X, y, family="cox", nfolds=10)
  lp <- as.numeric(predict(cv, newx=X, s="lambda.min", type="link")); if(sd(lp)==0) return(NA)
  as.numeric(summary(coxph(y ~ lp))$concordance[1])
}

## (1) head-to-head (full subsets) — both apparent (paper style) and out-of-fold CV (honest)
a_cyto<-appC(cyto); a_metab<-appC(metab); a_full<-appC(c(cyto,metab))
c_cyto  <- oofC(cyto)
c_metab <- oofC(metab)
c_full  <- oofC(c(cyto, metab))
cat(sprintf("\n[0] Apparent C-index (anchor; full should ~match paper 0.68)\n  cyto %.3f | metab %.3f | full %.3f\n", a_cyto,a_metab,a_full))
cat(sprintf("\n[1] Out-of-fold CV C-index (event-weighted mean of per-fold C)\n  cytoskeletal (%d genes): %.3f\n  metabolic   (%d genes): %.3f\n  full set    (%d genes): %.3f\n",
            length(cyto), c_cyto, length(metab), c_metab, length(cyto)+length(metab), c_full))

## (2) gene-count-matched control: cyto randomly down-sampled to size of metab
k <- length(metab); B <- 30
set.seed(2); cs <- replicate(B, oofC(sample(cyto, k), seed=sample(1e6,1)))
cat(sprintf("\n[2] Size-matched (cyto down-sampled to %d genes, %d draws)\n  cyto-%d mean C = %.3f (sd %.3f)\n  metab-%d   C = %.3f\n",
            k, B, k, mean(cs,na.rm=TRUE), sd(cs,na.rm=TRUE), k, c_metab))

## (3) per-gene univariable Cox C-index by class
uniC <- function(gn){ x <- Xall[,gn]; if(sd(x)==0) return(c(NA,NA))
  fit <- coxph(y ~ x); c(summary(fit)$concordance[1], summary(fit)$coefficients[1,5]) }
tab <- t(sapply(c(cyto,metab), uniC)); colnames(tab) <- c("C","p")
cls <- c(rep("cyto",length(cyto)), rep("metab",length(metab)))
cat(sprintf("\n[3] Per-gene univariable Cox\n  cyto : mean C = %.3f | sig(p<0.05) = %d/%d\n  metab: mean C = %.3f | sig(p<0.05) = %d/%d\n",
            mean(tab[cls=="cyto","C"],na.rm=TRUE), sum(tab[cls=="cyto","p"]<0.05,na.rm=TRUE), length(cyto),
            mean(tab[cls=="metab","C"],na.rm=TRUE), sum(tab[cls=="metab","p"]<0.05,na.rm=TRUE), length(metab)))

## save source data
out <- data.frame(gene=rownames(tab), class=cls, uni_C=tab[,"C"], uni_p=tab[,"p"])
write.csv(out, "results/subset_cyto_vs_metab_pergene.csv", row.names=FALSE)
saveRDS(list(a_cyto=a_cyto,a_metab=a_metab,a_full=a_full,c_cyto=c_cyto,c_metab=c_metab,c_full=c_full,
             cyto_matched=cs,pergene=out), "results/subset_cyto_vs_metab.rds")
cat("\nsaved: results/subset_cyto_vs_metab_pergene.csv + .rds\n")

## ---- Supplementary Fig. S3 ----
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (!exists("theme_nat")) source("code/kirc/_style.R")
th <- theme_nat() + theme(legend.position="top", legend.title=element_blank())
dA <- data.frame(set=factor(rep(c("Cytoskeletal","Metabolic","Full"),2),levels=c("Cytoskeletal","Metabolic","Full")),
                 metric=rep(c("Apparent","Out-of-fold CV"),each=3), C=c(a_cyto,a_metab,a_full,c_cyto,c_metab,c_full))
pA <- ggplot(dA,aes(set,C,fill=metric))+geom_col(position=position_dodge(.7),width=.65)+
  geom_hline(yintercept=.5,linetype=2,colour="grey40")+
  geom_text(aes(label=sprintf("%.2f",C)),position=position_dodge(.7),vjust=-.45,size=2.25)+
  scale_fill_manual(values=c(Apparent="#bdc3c7","Out-of-fold CV"="#34495e"))+
  coord_cartesian(ylim=c(.4,.78))+labs(title="Gene-class subset C-index",x=NULL,y="C-index")+th
dB <- data.frame(set=factor(c("Cytoskeletal\n(9, matched)","Metabolic\n(9)"),levels=c("Cytoskeletal\n(9, matched)","Metabolic\n(9)")),
                 C=c(mean(cs,na.rm=TRUE),c_metab), sd=c(sd(cs,na.rm=TRUE),NA))
pB <- ggplot(dB,aes(set,C,fill=set))+geom_col(width=.6)+
  geom_hline(yintercept=.5,linetype=2,colour="grey40")+geom_text(aes(label=sprintf("%.2f",C)),vjust=-.8,size=2.9)+
  scale_fill_manual(values=c("#C0392B","#2C7FB8"),guide="none")+
  coord_cartesian(ylim=c(.4,.78))+labs(title="Gene-count-matched (9 vs 9)",x=NULL,y="Out-of-fold C-index")+th
dC <- out; dC$class <- factor(ifelse(dC$class=="cyto","Cytoskeletal","Metabolic"))
pC <- ggplot(dC,aes(class,uni_C,fill=class))+geom_boxplot(width=.5,outlier.shape=NA,alpha=.5)+
  geom_jitter(width=.12,size=1.6)+geom_hline(yintercept=.5,linetype=2,colour="grey40")+
  scale_fill_manual(values=c(Cytoskeletal="#C0392B",Metabolic="#2C7FB8"),guide="none")+
  labs(title="Per-gene univariable C",x=NULL,y="Univariable C-index")+th
fig <- pA | pB | pC
if (!exists("ASSEMBLE_MODE")) {
  ggsave("results/SuppFigS3_subset_v2.png", fig, width=11, height=3.4, dpi=300)
  ggsave("results/SuppFigS3_subset_v2.svg", fig, width=11, height=3.4)
  cat("saved figure: results/SuppFigS3_subset_v2.{png,svg}\n")
}
