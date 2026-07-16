## fig_cellchat.R — cell-cell communication in ccRCC single-cell (GSE159115), CellChat.
## Deepens the single-cell story: which compartments signal, and how the tumour
## compartment relates to the disulfidptosis programme. Real data only.
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(patchwork); library(ggplot2)
})
set.seed(1)
source("code/kirc/_style.R")

if (!file.exists("results/kirc_cellchat.rds")) {
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
}

## ================= FIGURE =================
suppressPackageStartupMessages(library(ggplotify))
cc <- readRDS("results/kirc_cellchat.rds")
gcol <- c(B_Plasma="#8C6BB1", Endothelial="#4EB3D3", Fibroblast="#FDBB84",
          Myeloid="#5AAE61", T_NK="#2166AC", Tumor_ccRCC="#B2182B")
gcol <- gcol[levels(cc@idents)]
gg <- cc@net$count

## A: number of interactions ; B: interaction strength (circle nets)
## Rebuild the circular networks with fixed node and label positions. CellChat's
## base plot places labels on top of nodes for this dense six-group network.
circle_layout <- function(groups) {
  data.frame(
    group = groups,
    x = c(B_Plasma=0.96, Endothelial=0.40, Fibroblast=-0.45,
          Myeloid=-0.98, T_NK=-0.45, Tumor_ccRCC=0.38)[groups],
    y = c(B_Plasma=0.03, Endothelial=0.82, Fibroblast=0.82,
          Myeloid=0.02, T_NK=-0.82, Tumor_ccRCC=-0.82)[groups],
    stringsAsFactors = FALSE
  )
}
label_layout <- function(pos) {
  offs <- data.frame(
    group = pos$group,
    dx = c(B_Plasma=0.37, Endothelial=0.28, Fibroblast=-0.28,
           Myeloid=-0.32, T_NK=-0.14, Tumor_ccRCC=0.10)[pos$group],
    dy = c(B_Plasma=0.02, Endothelial=0.30, Fibroblast=0.30,
           Myeloid=-0.07, T_NK=-0.34, Tumor_ccRCC=-0.46)[pos$group],
    hjust = c(B_Plasma=0, Endothelial=0, Fibroblast=1,
              Myeloid=1, T_NK=1, Tumor_ccRCC=0)[pos$group],
    stringsAsFactors = FALSE
  )
  merge(pos, offs, by="group")
}
make_circle_net <- function(mat, title, x_shift=0, right_space=1.86, title_hjust=.52) {
  groups <- names(gcol)
  mat <- as.matrix(mat)[groups, groups, drop=FALSE]
  diag(mat) <- 0
  idx <- which(mat > 0, arr.ind=TRUE)
  edges <- data.frame(source=rownames(mat)[idx[,1]], target=colnames(mat)[idx[,2]],
                      value=mat[idx], stringsAsFactors=FALSE)
  edges$a <- ifelse(match(edges$source, groups) < match(edges$target, groups),
                    edges$source, edges$target)
  edges$b <- ifelse(match(edges$source, groups) < match(edges$target, groups),
                    edges$target, edges$source)
  edges <- aggregate(value ~ a + b, edges, sum)
  pos <- circle_layout(groups)
  pos$x <- pos$x + x_shift
  edges <- merge(edges, pos, by.x="a", by.y="group")
  edges <- merge(edges, pos, by.x="b", by.y="group", suffixes=c(".a",".b"))

  node_n <- as.numeric(table(cc@idents)[groups])
  nodes <- transform(pos, n=node_n,
                     node_size=4.7 + 6.1 * (sqrt(node_n) - min(sqrt(node_n))) /
                       diff(range(sqrt(node_n))))
  labs <- label_layout(pos)

  ggplot() +
    geom_curve(data=edges, aes(x=x.a, y=y.a, xend=x.b, yend=y.b, linewidth=value),
               curvature=.16, colour="#808080", alpha=.38, lineend="round") +
    geom_point(data=nodes, aes(x=x, y=y, size=node_size, fill=group),
               shape=21, colour="white", stroke=.35, alpha=.97) +
    geom_text(data=labs, aes(x=x + dx, y=y + dy, label=group, hjust=hjust),
              family="Arial", size=2.35, colour="black", lineheight=.92) +
    scale_fill_manual(values=gcol, guide="none") +
    scale_size_identity() +
    scale_linewidth(range=c(.18, 1.55), guide="none") +
    coord_equal(xlim=c(-1.70, right_space), ylim=c(-1.42, 1.32), clip="off") +
    labs(title=title) +
    theme_void(base_family="Arial") +
    theme(plot.title=element_text(size=7, face="bold", hjust=title_hjust),
          plot.margin=margin(4, 8, 4, 0))
}
pA <- make_circle_net(cc@net$count, "Number of interactions",
                      x_shift=-.78, right_space=2.05, title_hjust=.18)
pB <- make_circle_net(cc@net$weight, "Interaction strength")

## C: dominant senders vs receivers
role_df <- data.frame(group=names(gcol),
                      outgoing=rowSums(cc@net$weight[names(gcol), names(gcol)]),
                      incoming=colSums(cc@net$weight[names(gcol), names(gcol)]),
                      n=as.numeric(table(cc@idents)[names(gcol)]),
                      stringsAsFactors=FALSE)
role_lab <- transform(role_df,
  dx = c(B_Plasma=.40, Endothelial=-.28, Fibroblast=.34,
         Myeloid=.52, T_NK=.30, Tumor_ccRCC=-.42)[group],
  dy = c(B_Plasma=.12, Endothelial=.02, Fibroblast=-.08,
         Myeloid=.18, T_NK=-.02, Tumor_ccRCC=.18)[group],
  hjust = c(B_Plasma=0, Endothelial=1, Fibroblast=0,
            Myeloid=0, T_NK=0, Tumor_ccRCC=1)[group])
pC <- ggplot(role_df, aes(outgoing, incoming)) +
  geom_point(aes(size=n, fill=group), shape=21, colour="white", stroke=.35, alpha=.98) +
  geom_text(data=role_lab, aes(x=outgoing + dx, y=incoming + dy,
                               label=group, hjust=hjust, colour=group),
            family="Arial", size=2.35, show.legend=FALSE) +
  scale_fill_manual(values=gcol, guide="none") +
  scale_colour_manual(values=gcol, guide="none") +
  scale_size_continuous(range=c(3.8, 7.2), guide="none") +
  coord_cartesian(xlim=c(1.75, 7.85), ylim=c(1.32, 7.55), clip="off") +
  labs(title="Dominant senders / receivers",
       x="Outgoing interaction strength", y="Incoming interaction strength") +
  theme_nat() +
  theme(plot.margin=margin(4, 9, 4, 5))

## D: MHC-II signalling as a clean sender x receiver tile heatmap (ties to ACTB-KO / antigen presentation)
mhc <- cc@netP$prob[,, "MHC-II"]
dfh <- expand.grid(Sender=rownames(mhc), Receiver=colnames(mhc), KEEP.OUT.ATTRS=FALSE)
dfh$prob <- as.vector(mhc)
pD <- ggplot(dfh, aes(Receiver, Sender, fill=prob)) +
  geom_tile(colour="white", linewidth=.5) +
  scale_fill_gradient(low="white", high=PAL$high, name="Comm.\nprob.",
                      guide=guide_colourbar(barwidth=.5, barheight=3)) +
  labs(x="Receiver", y="Sender", title="MHC-II signalling") +
  coord_equal() + theme_nat() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_text(size=6))

## E: tumour -> immune ligand-receptor bubble (top pathways)
imm <- intersect(c("Myeloid","T_NK","B_Plasma"), levels(cc@idents))
pE_raw <- netVisual_bubble(cc, sources.use="Tumor_ccRCC", targets.use=imm,
        signaling=c("MHC-II","MHC-I","SPP1","VEGF","MIF","APP","GALECTIN"),
        remove.isolate=TRUE, return.data=FALSE)
bub <- pE_raw$data
bub$target <- factor(bub$target, levels=imm)
bub$interaction_name_2 <- factor(bub$interaction_name_2,
                                 levels=rev(unique(as.character(bub$interaction_name_2))))
pE <- ggplot(bub, aes(target, interaction_name_2)) +
  geom_point(aes(size=prob, fill=prob), shape=21, colour="grey25",
             stroke=.12, alpha=.96) +
  scale_fill_gradient(low="white", high=PAL$high, name="Comm.\nprob.",
                      guide=guide_colourbar(barwidth=.55, barheight=3.2)) +
  scale_size_continuous(range=c(.65, 2.35), guide="none") +
  scale_x_discrete(position="top", expand=expansion(add=.62)) +
  scale_y_discrete(expand=expansion(add=.45)) +
  labs(x=NULL, y=NULL, title="Tumour-to-immune signalling") +
  theme_nat(base=6.1) +
  theme(axis.text.x=element_text(angle=45, hjust=0, vjust=.5, size=5.2),
        axis.text.y=element_text(size=4.7),
        axis.line=element_blank(), axis.ticks=element_blank(),
        panel.grid.major=element_line(linewidth=.18, colour="#E5E5E5"),
        legend.position="right", legend.key.size=unit(2.8,"mm"),
        plot.margin=margin(4, 2, 3, 3))

fig <- (pA | pB | pC) / (pD | pE) + plot_layout(heights=c(1,1.1)) +
  plot_annotation(tag_levels="A") & theme(plot.tag=element_text(face="bold",size=9))
save_pub(fig, "results/Fig_CellChat", w_mm=183, h_mm=150)
tryCatch({
  ragg::agg_tiff("results/Fig_CellChat.tiff", width=183/25.4, height=150/25.4,
                 units="in", res=300, compression="lzw")
  print(fig)
  grDevices::dev.off()
}, error=function(e) message("  (TIFF skipped: ", conditionMessage(e), ")"))
cat("saved results/Fig_CellChat\n")
