## fig_cellchat.R — cell-cell communication in ccRCC single-cell (GSE159115), CellChat.
## Deepens the single-cell story: which compartments signal, and how the tumour
## compartment relates to the disulfidptosis programme. Real data only.
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(patchwork); library(ggplot2)
})
set.seed(1)
source("code/kirc/_style.R")

s <- readRDS("results/kirc_scRNA.rds")
## log-normalised data matrix (robust to Seurat v4/v5)
data.input <- tryCatch(GetAssayData(s, assay="RNA", slot="data"),
                       error=function(e) SeuratObject::LayerData(s, assay="RNA", layer="data"))
meta <- data.frame(labels=as.character(s$cell_type), row.names=colnames(s))

cc <- createCellChat(object=data.input, meta=meta, group.by="labels")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type="triMean")
cc <- filterCommunication(cc, min.cells=10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
cc <- netAnalysis_computeCentrality(cc, slot.name="netP")

cat("groups:", paste(levels(cc@idents), collapse=", "), "\n")
cat("significant pathways:\n"); print(cc@netP$pathways)
saveRDS(cc, "results/kirc_cellchat.rds")
cat("saved results/kirc_cellchat.rds\n")

## ================= FIGURE =================
suppressPackageStartupMessages(library(ggplotify))
cc <- readRDS("results/kirc_cellchat.rds")
gcol <- c(B_Plasma="#8C6BB1", Endothelial="#4EB3D3", Fibroblast="#FDBB84",
          Myeloid="#5AAE61", T_NK="#2166AC", Tumor_ccRCC="#B2182B")
gcol <- gcol[levels(cc@idents)]
gg <- cc@net$count

## A: number of interactions ; B: interaction strength (circle nets)
pA <- as.ggplot(function() netVisual_circle(cc@net$count, vertex.weight=as.numeric(table(cc@idents)),
        weight.scale=TRUE, label.edge=FALSE, color.use=gcol, title.name="Number of interactions",
        vertex.label.cex=.7))
pB <- as.ggplot(function() netVisual_circle(cc@net$weight, vertex.weight=as.numeric(table(cc@idents)),
        weight.scale=TRUE, label.edge=FALSE, color.use=gcol, title.name="Interaction strength",
        vertex.label.cex=.7))

## C: dominant senders vs receivers
pC <- netAnalysis_signalingRole_scatter(cc, color.use=gcol) +
      labs(title="Dominant senders / receivers") + theme_nat() +
      theme(legend.position="none")

## D: MHC-II signalling heatmap (sender x receiver; ties to ACTB-KO / antigen presentation)
pD <- as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(
        netVisual_heatmap(cc, signaling="MHC-II", color.heatmap="Reds",
        color.use=gcol, title.name="MHC-II signalling"))))

## E: tumour -> immune ligand-receptor bubble (top pathways)
imm <- intersect(c("Myeloid","T_NK","B_Plasma"), levels(cc@idents))
pE <- netVisual_bubble(cc, sources.use="Tumor_ccRCC", targets.use=imm,
        signaling=c("MHC-II","MHC-I","SPP1","VEGF","MIF","APP","GALECTIN"),
        remove.isolate=TRUE, return.data=FALSE) +
      labs(title="Tumour-to-immune signalling") + theme_nat() +
      theme(axis.text.x=element_text(angle=45,hjust=1,size=5), axis.text.y=element_text(size=5),
            legend.position="right", legend.key.size=unit(3,"mm"))

fig <- (pA | pB | pC) / (pD | pE) + plot_layout(heights=c(1,1.1)) +
  plot_annotation(tag_levels="A") & theme(plot.tag=element_text(face="bold",size=9))
save_pub(fig, "results/Fig_CellChat", w_mm=183, h_mm=150)
cat("saved results/Fig_CellChat\n")
