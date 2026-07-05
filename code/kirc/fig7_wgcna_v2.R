## Fig 7 升级版: WGCNA 模块-性状 ComplexHeatmap(替换简陋 labeledHeatmap)
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(WGCNA); library(ComplexHeatmap); library(circlize); library(grid) })
w <- readRDS("results/kirc_wgcna.rds")
MEs <- orderMEs(w$net$MEs)
tr  <- w$traits[rownames(MEs), , drop=FALSE]
mtc <- cor(MEs, tr, use="p"); mtp <- corPvalueStudent(mtc, nrow(MEs))
rownames(mtc) <- gsub("ME","M", rownames(mtc))
af <- gpar(fontsize=6, fontfamily="Arial")
cor_ramp <- circlize::colorRamp2(c(-.6,0,.6), c(PAL$low,"white",PAL$high))   # 相关性专用色阶,不偏淡
sig <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))  # 单行:相关+显著性星号(避免两行数字上下叠印)
ht <- Heatmap(mtc, name="Correlation", col=cor_ramp, cluster_rows=FALSE, cluster_columns=FALSE,
  width = unit(15,"mm")*ncol(mtc), height = unit(4.2,"mm")*nrow(mtc),  # 固定格子尺寸,避免拼装时被宽扁槽拉扁
  rect_gp=gpar(col="white", lwd=.6),
  cell_fun=function(j,i,x,y,wd,ht2,fill) grid.text(sprintf("%.2f%s", mtc[i,j], sig(mtp[i,j])), x, y,
    gp=gpar(fontsize=5.5, fontfamily="Arial", col=ifelse(abs(mtc[i,j])>.55,"white","black"))),
  row_names_gp=gpar(fontsize=6,fontfamily="Arial"), column_names_gp=gpar(fontsize=6.5,fontfamily="Arial"),
  column_names_rot=0, column_names_centered=TRUE,
  column_title="WGCNA module-trait relationships", column_title_gp=gpar(fontsize=7,fontface="bold",fontfamily="Arial"),
  heatmap_legend_param=list(title_gp=af, labels_gp=af, legend_height=unit(18,"mm")))
pHt <- wrap_elements(
  full = grid.grabExpr(draw(ht, padding=unit(c(4, 7, 7, 10), "mm"))),
  clip = FALSE
)     # 转 grob,显式 padding 避免拼装后标题/边缘标签被裁掉
if (!exists("ASSEMBLE_MODE")) save_pub(pHt, "results/Fig8_WGCNA_v2", w_mm = 120, h_mm = 130)
cat("Fig7 升级版完成\n")
