## fig_ml_shap.R — Multi-algorithm ML + SHAP explainability for the 12-gene signature
## Mirrors the house-style explainable-ML block (GBM relative influence, RF importance,
## XGBoost SHAP importance + beeswarm). Real TCGA-KIRC data only.
suppressPackageStartupMessages({
  library(data.table); library(gbm); library(ranger); library(xgboost)
  library(SHAPforxgboost); library(pROC); library(ggplot2); library(patchwork)
})
set.seed(42)
source("code/kirc/_style.R")

## --- 1. data: 12 signature genes x 531 prognosis samples ---
prog <- readRDS("data/raw/kirc_prognosis.rds")
genes <- prog$lasso
expr  <- fread("data/raw/KIRC_HiSeqV2.gz")            # gene x sample, log2(norm+1)
gm    <- as.data.frame(expr); rownames(gm) <- gm[[1]]; gm[[1]] <- NULL
stopifnot(all(genes %in% rownames(gm)))
X <- t(gm[genes, prog$df$sample, drop=FALSE]); X <- as.data.frame(X)
d <- cbind(prog$df[,c("sample","OS","OS.time")], X)

## --- 2. clean binary outcome: died (event) vs confirmed long-term survivor (>3y alive) ---
d$label <- NA
d$label[d$OS==1] <- 1L                                # died
d$label[d$OS==0 & d$OS.time > 1095] <- 0L             # alive & followed >3y
dd <- d[!is.na(d$label), ]
cat(sprintf("ML cohort: %d (died=%d, long-term survivor=%d) of %d\n",
            nrow(dd), sum(dd$label==1), sum(dd$label==0), nrow(d)))
Xm <- as.matrix(dd[, genes]); y <- dd$label

## --- 3. 5-fold CV AUC for each algorithm (show these aren't overfit toys) ---
folds <- sample(rep(1:5, length.out=nrow(dd)))
cv_auc <- function(fitpred){
  pr <- numeric(nrow(dd))
  for(k in 1:5){ tr <- folds!=k; te <- !tr; pr[te] <- fitpred(Xm[tr,],y[tr],Xm[te,]) }
  as.numeric(pROC::auc(y, pr, quiet=TRUE))
}
auc_gbm <- cv_auc(function(xtr,ytr,xte){
  m <- gbm::gbm.fit(xtr, ytr, distribution="bernoulli", n.trees=500,
                    interaction.depth=3, shrinkage=0.01, verbose=FALSE)
  gbm::predict.gbm(m, as.data.frame(xte), n.trees=500, type="response") })
auc_rf <- cv_auc(function(xtr,ytr,xte){
  m <- ranger::ranger(x=as.data.frame(xtr), y=factor(ytr), probability=TRUE, num.trees=800)
  predict(m, as.data.frame(xte))$predictions[,"1"] })
xgb_params <- list(objective="binary:logistic", max_depth=3, eta=0.05)
auc_xgb <- cv_auc(function(xtr,ytr,xte){
  m <- xgboost::xgb.train(params=xgb_params, data=xgboost::xgb.DMatrix(xtr,label=ytr),
                          nrounds=120, verbose=0)
  predict(m, xgboost::xgb.DMatrix(xte)) })
cat(sprintf("5-fold CV AUC | GBM %.3f | RF %.3f | XGBoost %.3f\n", auc_gbm, auc_rf, auc_xgb))

## --- 4. full-data fits for importance / SHAP ---
gbm_full <- gbm::gbm.fit(Xm, y, distribution="bernoulli", n.trees=800,
                         interaction.depth=3, shrinkage=0.01, verbose=FALSE)
gbm_imp <- summary(gbm_full, plotit=FALSE)                       # rel.inf
rf_full <- ranger::ranger(x=as.data.frame(Xm), y=factor(y), probability=TRUE,
                          num.trees=1000, importance="permutation")
rf_imp  <- ranger::importance(rf_full)
xgb_full <- xgboost::xgb.train(params=xgb_params, data=xgboost::xgb.DMatrix(Xm,label=y),
                               nrounds=200, verbose=0)
shap <- SHAPforxgboost::shap.values(xgb_model=xgb_full, X_train=Xm)
shap_imp <- sort(shap$mean_shap_score, decreasing=TRUE)

cat("\nGBM relative influence (top):\n"); print(round(head(gbm_imp$rel.inf[order(-gbm_imp$rel.inf)],6),2))
cat("mean|SHAP| (top):\n"); print(round(head(shap_imp,6),4))

saveRDS(list(genes=genes, dd=dd, gbm_imp=gbm_imp, rf_imp=rf_imp,
             shap=shap, xgb=xgb_full, Xm=Xm,
             auc=c(GBM=auc_gbm, RF=auc_rf, XGBoost=auc_xgb),
             coef=prog$coef), "results/kirc_ml_shap.rds")
cat("\nsaved results/kirc_ml_shap.rds\n")

## ================= FIGURE =================
res <- readRDS("results/kirc_ml_shap.rds")
lvl <- names(sort(res$shap$mean_shap_score))          # gene order by SHAP (low->high)

## Panel A: cross-validated AUC across algorithms
dfa <- data.frame(alg=factor(names(res$auc), levels=names(res$auc)), auc=as.numeric(res$auc))
pA <- ggplot(dfa, aes(alg, auc, fill=alg)) +
  geom_col(width=.62, colour=NA) +
  geom_text(aes(label=sprintf("%.3f",auc)), vjust=-0.4, size=2.4) +
  scale_fill_manual(values=c(GBM=PAL$high, RF="#5AAE61", XGBoost=PAL$low), guide="none") +
  scale_y_continuous(limits=c(0,0.8), expand=expansion(c(0,.08))) +
  labs(x=NULL, y="5-fold CV AUC", title="Multi-algorithm cross-validation") + theme_nat()

## Panel B: GBM relative influence
gi <- res$gbm_imp; gi <- gi[order(gi$rel.inf),]; gi$var <- factor(gi$var, levels=gi$var)
pB <- ggplot(gi, aes(rel.inf, var)) +
  geom_col(fill=PAL$high, width=.7) +
  geom_text(aes(label=sprintf("%.1f",rel.inf)), hjust=-0.15, size=2.2) +
  scale_x_continuous(expand=expansion(c(0,.12))) +
  labs(x="Relative influence", y=NULL, title="GBM gene importance") + theme_nat()

## Panel C: mean |SHAP| importance
sc <- data.frame(gene=names(res$shap$mean_shap_score), val=as.numeric(res$shap$mean_shap_score))
sc <- sc[order(sc$val),]; sc$gene <- factor(sc$gene, levels=sc$gene)
pC <- ggplot(sc, aes(val, gene)) +
  geom_col(fill=PAL$low, width=.7) +
  geom_text(aes(label=sprintf("%.3f",val)), hjust=-0.15, size=2.2) +
  scale_x_continuous(expand=expansion(c(0,.15))) +
  labs(x="mean(|SHAP value|)", y=NULL, title="XGBoost SHAP importance") + theme_nat()

## Panel D: SHAP beeswarm
shap_long <- SHAPforxgboost::shap.prep(xgb_model=res$xgb, X_train=res$Xm)
pD <- SHAPforxgboost::shap.plot.summary(shap_long, scientific=FALSE) +
  scale_colour_gradient(low=PAL$low, high=PAL$high, breaks=c(0,1),
                        labels=c("Low","High"), guide=guide_colourbar(barwidth=.5,barheight=3)) +
  labs(title="SHAP summary (directionality)") + theme_nat() +
  theme(legend.position="right", legend.title=element_text(size=6))

fig <- (pA | pB) / (pC | pD) + plot_layout(heights=c(1,1.15)) +
  plot_annotation(tag_levels="A") & theme(plot.tag=element_text(face="bold", size=9))
save_pub(fig, "results/Fig_ML_SHAP", w_mm=183, h_mm=150)
cat("saved results/Fig_ML_SHAP\n")
