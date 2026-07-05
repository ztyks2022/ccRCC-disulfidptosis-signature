## GDSC 药敏预测(oncoPredict calcPhenotype)-> results/gdsc/calcPhenotype_Output/
suppressPackageStartupMessages({ library(oncoPredict); library(data.table) })
GE <- readRDS("data/raw/gdsc/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GR <- readRDS("data/raw/gdsc/GDSC2_Res.rds")
m <- fread("data/raw/KIRC_HiSeqV2.gz"); g <- m[[1]]; expr <- as.matrix(m[,-1]); rownames(expr) <- g
tum <- expr[, substr(colnames(expr),14,15) == "01"]
dir.create("results/gdsc", showWarnings = FALSE, recursive = TRUE)
setwd("results/gdsc")                         # calcPhenotype 输出写到当前目录
calcPhenotype(trainingExprData = GE, trainingPtype = GR, testExprData = tum,
              batchCorrect = "eb", powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2, minNumSamples = 10,
              printOutput = FALSE, removeLowVaringGenesFrom = "rawData")
cat("GDSC calcPhenotype 完成 -> results/gdsc/calcPhenotype_Output/DrugPredictions.csv\n")
