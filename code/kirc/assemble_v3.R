## assemble_v3.R — consolidate the 12 fragmented v2 figures into 7 themed,
## continuously-lettered submission figures. Contributing fig*_v2.R scripts are sourced
## with ASSEMBLE_MODE=TRUE so they build panels but skip their own save_pub().
## ggplot panels compose directly; base-R panels (nomogram/calibration/forest) are
## wrapped with ggplotify::as.ggplot so everything lives in one patchwork layout.
ASSEMBLE_MODE <- TRUE
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(patchwork); library(ggplotify) })
TAG <- function(g) g + plot_annotation(tag_levels="A") &
  theme(plot.tag = element_text(size=11, face="bold", family="Arial"))

grab <- function(script, ...) {            # isolated source -> named panel objects
  env <- new.env(); env$ASSEMBLE_MODE <- TRUE
  sys.source(script, envir = env)
  vars <- c(...); setNames(lapply(vars, function(v) get(v, envir = env)), vars)
}

## ===== Fig 3 — Validation & benchmarking (old Fig5 validation + ClearCode34) =====
v  <- grab("code/kirc/fig11_validation.R", "pA", "pB", "pC")
cc <- grab("code/kirc/supp_cindex.R", "p")
fig3 <- TAG( (v$pA | v$pB | v$pC) / (cc$p) + plot_layout(heights=c(1, 1.05)) )
save_pub(fig3, "results/Fig3_validation_v3", w_mm=183, h_mm=150)
cat("Fig3 (validation + ClearCode34) done\n")

## ===== Fig 6 — Cellular & network context (single-cell + WGCNA + candidates) =====
sc <- grab("code/kirc/fig6_singlecell_v2.R", "pA", "pB", "pC", "pD")
wg <- grab("code/kirc/fig7_wgcna_v2.R", "pHt")
cd <- grab("code/kirc/fig8_candidate_v2.R", "pV", "pK")
scA <- wrap_elements(
  full=grid::grobTree(ggplot2::ggplotGrob(sc$pA),
                      vp=grid::viewport(x=.52, y=.49, width=.90, height=.90)),
  clip=FALSE
)
fig6_top <- scA | sc$pB | sc$pC
fig6_mid <- (sc$pD | wg$pHt) + plot_layout(widths=c(1.20, 1.55))
fig6_bottom <- (cd$pV | cd$pK) + plot_layout(widths=c(1.28, 1))
fig6 <- TAG( fig6_top / fig6_mid / fig6_bottom +
             plot_layout(heights=c(1.13, 1.05, 1)) )
save_pub(fig6, "results/Fig6_context_v3", w_mm=183, h_mm=215)
cat("Fig6 (single-cell + WGCNA + candidates) done\n")

## ===== Fig 4 — Independent prognosis & molecular subtypes =====
## nomogram/calibration/forest are base-R; source fig4 to global so its objects
## (nom, d, hr, rn, us, cols, labs, datadist) are available to the as.ggplot closures.
source("code/kirc/fig4_nomogram_v2.R")              # computes nom,d,hr,rn,us,cols,labs (+saves Fig4A/BC, harmless)
build_nomogram_panel <- function(nom) {
  row_y <- c("Points"=8, "Risk"=7, "Stage"=6, "Grade"=5, "Age"=4,
             "Total Points"=3, "1-year OS"=2, "3-year OS"=1, "5-year OS"=0)
  total_max <- max(nom$total.points$x, na.rm=TRUE)
  to_panel_x <- function(x, xmax=100) x / xmax * 100
  tick_len <- .12
  axis_df <- data.frame(
    label=names(row_y), y=unname(row_y),
    x0=c(0, 0, min(nom$Stage$points), min(nom$Grade$points), min(nom$Age$points),
         0, min(nom[["1-year OS"]]$x), min(nom[["3-year OS"]]$x), min(nom[["5-year OS"]]$x)),
    x1=c(100, max(nom$Risk$points), max(nom$Stage$points), max(nom$Grade$points),
         max(nom$Age$points), total_max, max(nom[["1-year OS"]]$x),
         max(nom[["3-year OS"]]$x), max(nom[["5-year OS"]]$x)),
    xmax=c(100, 100, 100, 100, 100, total_max, total_max, total_max, total_max)
  )
  axis_df$x0p <- to_panel_x(axis_df$x0, axis_df$xmax)
  axis_df$x1p <- to_panel_x(axis_df$x1, axis_df$xmax)
  point_ticks <- data.frame(row="Points", value=seq(0,100,20),
                            label=as.character(seq(0,100,20)), y=row_y["Points"], xmax=100)
  point_minor <- data.frame(value=setdiff(seq(0,100,2), seq(0,100,10)), y=row_y["Points"], xmax=100)
  predictor_ticks <- rbind(
    data.frame(row="Risk", value=nom$Risk$points, label=nom$Risk$Risk, y=row_y["Risk"], xmax=100),
    data.frame(row="Stage", value=nom$Stage$points, label=nom$Stage$Stage, y=row_y["Stage"], xmax=100),
    data.frame(row="Grade", value=nom$Grade$points, label=nom$Grade$Grade, y=row_y["Grade"], xmax=100),
    data.frame(row="Age", value=nom$Age$points, label=nom$Age$Age, y=row_y["Age"], xmax=100)
  )
  predictor_ticks$label <- as.character(predictor_ticks$label)
  predictor_ticks <- subset(
    predictor_ticks,
    row %in% c("Stage", "Grade") |
      (row == "Risk" & label %in% c("0.5", "1", "2", "3", "4", "4.5")) |
      (row == "Age" & label %in% c("25", "35", "45", "55", "65", "75", "85"))
  )
  total_ticks <- data.frame(row="Total Points", value=seq(0, total_max, 40),
                            label=seq(0, total_max, 40), y=row_y["Total Points"], xmax=total_max)
  total_minor <- data.frame(value=setdiff(seq(0,total_max,5), nom$total.points$x),
                            y=row_y["Total Points"], xmax=total_max)
  surv_ticks <- do.call(rbind, lapply(c("1-year OS", "3-year OS", "5-year OS"), function(row) {
    data.frame(row=row, value=nom[[row]]$x, label=nom[[row]]$fat, y=row_y[row], xmax=total_max)
  }))
  surv_ticks$label <- as.character(surv_ticks$label)
  surv_ticks <- subset(surv_ticks, label %in% c("0.9", "0.7", "0.5", "0.3", "0.1"))
  major_ticks <- rbind(point_ticks, predictor_ticks, total_ticks, surv_ticks)
  major_ticks$x <- to_panel_x(major_ticks$value, major_ticks$xmax)
  point_minor$x <- to_panel_x(point_minor$value, point_minor$xmax)
  total_minor$x <- to_panel_x(total_minor$value, total_minor$xmax)

  ggplot() +
    geom_segment(data=data.frame(x=seq(0,100,10)),
                 aes(x=x, xend=x, y=row_y["Points"], yend=row_y["Age"]-.15),
                 colour=PAL$n_light, linewidth=.25) +
    geom_segment(data=axis_df, aes(x=x0p, xend=x1p, y=y, yend=y), linewidth=.45, colour="black") +
    geom_segment(data=major_ticks, aes(x=x, xend=x, y=y, yend=y-tick_len), linewidth=.35) +
    geom_segment(data=point_minor, aes(x=x, xend=x, y=y, yend=y-tick_len*.55), linewidth=.25) +
    geom_segment(data=total_minor, aes(x=x, xend=x, y=y, yend=y-tick_len*.55), linewidth=.25) +
    geom_text(data=data.frame(label=ifelse(names(row_y)=="Total Points", "Total\nPoints", names(row_y)),
                              y=unname(row_y), x=-5, size=2.35),
              aes(x=x, y=y, label=label, size=size), hjust=1, lineheight=.82,
              family="Arial", colour="black",
              show.legend=FALSE) +
    scale_size_identity() +
    geom_text(data=major_ticks, aes(x=x, y=y-.28, label=label),
              family="Arial", size=1.55, colour="black") +
  annotate("text", x=50, y=8.72, label="Nomogram",
           family="Arial", fontface="bold", size=2.75) +
    coord_cartesian(xlim=c(-22,114), ylim=c(-.45,8.95), clip="off") +
    theme_void(base_family="Arial") +
    theme(plot.margin=margin(3,4,3,4),
          panel.background=element_rect(fill="white", colour=NA),
          plot.background=element_rect(fill="white", colour=NA))
}
nomP <- build_nomogram_panel(nom)
nomP_shift <- wrap_elements(
  full = grid::grobTree(
    ggplot2::ggplotGrob(nomP),
    vp = grid::viewport(x=.54, y=.5, width=.88, height=1)
  ),
  clip = FALSE
)
## calibration + forest rebuilt as clean ggplot (avoids base-R label clipping in small cells)
caldf <- do.call(rbind, lapply(seq_along(us), function(i){
  fi <- rms::cph(Surv(OS.time,OS)~Risk+Stage+Grade+Age, data=d, x=TRUE,y=TRUE,surv=TRUE, time.inc=us[i])
  ca <- tryCatch(rms::calibrate(fi,u=us[i],cmethod="KM",m=60,B=120), error=function(e)NULL)
  if(is.null(ca)) return(NULL)
  data.frame(pred=ca[,"mean.predicted"], obs=ca[,"KM"], horizon=labs[i]) }))
caldf$horizon <- factor(caldf$horizon, levels=labs)
calP <- ggplot(caldf, aes(pred, obs, colour=horizon)) +
  geom_abline(linetype=2, colour=PAL$n_mid) + geom_line(linewidth=.5) + geom_point(size=1.3) +
  scale_colour_manual(values=setNames(cols, labs), name=NULL) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  labs(title="Calibration", x="Nomogram-predicted OS", y="Observed OS") +
  theme_nat() + theme(legend.position=c(.82,.2))
fdf <- data.frame(var=factor(rn, levels=rev(rn)), hr=hr[,1], lo=hr[,2], hi=hr[,3])
forP <- ggplot(fdf, aes(hr, var)) +
  geom_vline(xintercept=1, linetype=2, colour=PAL$n_mid) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=.18, colour=PAL$low, linewidth=.6) +
  geom_point(colour=PAL$high, size=2.2) + scale_x_log10() +
  labs(title="Multivariate Cox", x="Hazard ratio (95% CI)", y=NULL) + theme_nat()
sb <- grab("code/kirc/fig2_subtype_v2.R", "pA", "pB", "pC", "pD")
dca <- grab("code/kirc/fig_dca.R", "pDCA")

dca_compact <- dca$pDCA +
  scale_colour_manual(
    values = c("Nomogram (all)" = PAL$high, "Risk score" = PAL$accent3,
               "Clinical (stage+grade+age)" = PAL$low, "Treat all" = "grey55",
               "Treat none" = "grey20"),
    labels = c("Nomogram", "Risk", "Clinical", "Treat all", "Treat none"),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Nomogram (all)" = 1, "Risk score" = 1,
               "Clinical (stage+grade+age)" = 1, "Treat all" = 2,
               "Treat none" = 2),
    labels = c("Nomogram", "Risk", "Clinical", "Treat all", "Treat none"),
    name = NULL
  ) +
  labs(title = "Decision curve analysis") +
  theme(
    legend.position = c(.67, .80),
    legend.text = element_text(size = 4.8),
    legend.key.height = unit(2.2, "mm"),
    plot.title = element_text(size = 6.2, face = "bold")
  )

sb$pC <- sb$pC +
  theme(
    axis.text.x = element_text(size = 4.8, angle = 25, hjust = 1),
    legend.position = "none"
  )

row1 <- (nomP_shift | calP | forP | dca_compact) +
  plot_layout(widths = c(1.85, .78, .82, .88))

row2 <- (sb$pA | sb$pB | sb$pC | sb$pD) +
  plot_layout(widths = c(2.05, .78, .90, .80))

fig4 <- TAG(
  row1 / row2 +
    plot_layout(heights = c(1, 1))
) & theme(plot.margin = margin(8, 4, 4, 4))

save_pub(fig4, "results/Fig4_prognosis_subtype_v3", w_mm=220, h_mm=145)
cat("Fig4 (independent prognosis + subtypes) done\n")

## ===== Fig 7 — Therapeutic implications (drug + ACTB knockout) =====
d7 <- grab("code/kirc/fig9_drug_v2.R", "pA", "pB")
k7 <- grab("code/kirc/fig10_virtualKO_v2.R", "pA", "pB")
fig7 <- TAG( (d7$pA | d7$pB) / (k7$pA | k7$pB) + plot_layout(heights=c(1,1)) )
save_pub(fig7, "results/Fig7_therapy_v3", w_mm=183, h_mm=165)
cat("Fig7 (therapy) done\n")

## ===== Fig 1 — Disulfidptosis landscape (score / DEG / key genes / correlation) =====
l <- grab("code/kirc/fig1_landscape_v2.R", "pA", "pB", "pC", "pD")
fig1 <- TAG( (l$pA | l$pB) / l$pC / l$pD +
             plot_layout(heights=c(1, .8, 1.3)) )
save_pub(fig1, "results/Fig1_landscape_v3", w_mm=183, h_mm=183)
cat("Fig1 (landscape) done\n")

## ===== Fig 2 — Prognostic signature: HR forest / risk KM / time-ROC + gene-class robustness =====
pr <- grab("code/kirc/fig3_prognosis_v2.R", "pA", "pB", "pC")     # forest, risk-KM, time-ROC
ss <- grab("code/kirc/subset_cyto_vs_metab.R", "pA", "pB", "pC")  # cyto-vs-metab subset robustness
fig2 <- TAG( (pr$pA | pr$pB | pr$pC) / (ss$pA | ss$pB | ss$pC) + plot_layout(heights=c(1, 1)) )
save_pub(fig2, "results/Fig2_signature_v3", w_mm=183, h_mm=145)
cat("Fig2 (signature) done\n")

## ===== Fig 5 — Tumour immune microenvironment by risk =====
im  <- grab("code/kirc/fig5_immune_v2.R", "pA", "pB", "pC")
icb <- grab("code/kirc/fig_icb.R", "pICB")           # ICB-response signatures by risk
tmbg<- grab("code/kirc/fig_tmb.R", "pTMB")           # tumour mutational burden by risk
row_b <- (icb$pICB | tmbg$pTMB) + plot_layout(widths=c(3.2, 1))
fig5 <- TAG( (im$pA | im$pB | im$pC) / row_b + plot_layout(heights=c(1.08, 1)) )
save_pub(fig5, "results/Fig5_immune_v3", w_mm=183, h_mm=155)
cat("Fig5 (immune + ICB + TMB) done\n")

## ===== publish to submission package (guarantees code -> package one-to-one) =====
## Map each generated v3 figure to its submission name; copy PNG/SVG/PDF and
## regenerate the 300-dpi TIFF, so a single run of this script refreshes everything.
PUBMAP <- c(Fig1_landscape_v3="Figure1", Fig2_signature_v3="Figure2",
            Fig3_validation_v3="Figure3", Fig4_prognosis_subtype_v3="Figure4",
            Fig5_immune_v3="Figure5", Fig6_context_v3="Figure6", Fig7_therapy_v3="Figure7")
figdir <- "submission_BMC/figures"; tifdir <- "submission_BMC/figures_TIFF_300dpi"
dir.create(figdir, showWarnings=FALSE, recursive=TRUE); dir.create(tifdir, showWarnings=FALSE, recursive=TRUE)
for (src in names(PUBMAP)) {
  dst <- PUBMAP[[src]]
  for (ext in c("png","svg","pdf")) {
    f <- sprintf("results/%s.%s", src, ext)
    if (file.exists(f)) file.copy(f, sprintf("%s/%s.%s", figdir, dst, ext), overwrite=TRUE)
  }
  png <- sprintf("results/%s.png", src)
  if (file.exists(png) && requireNamespace("magick", quietly=TRUE)) {
    magick::image_write(magick::image_convert(magick::image_read(png), "tiff"),
                        path=sprintf("%s/%s.tiff", tifdir, dst), format="tiff", density="300x300")
  }
}
cat("published 7 figures -> submission_BMC/{figures, figures_TIFF_300dpi}\n")

## ===== supplementary figures (each = one standalone script) =====
SUPP <- list(
  list(script="code/kirc/supp_gsea_v2.R", out="results/SuppFig_GSEA_v2",  dst="Supplementary_Figure_S1"),
  list(script="code/kirc/fig12_network.R", out="results/Fig9_network_v2", dst="Supplementary_Figure_S2"))
sdir <- "submission_BMC/supplementary"; dir.create(sdir, showWarnings=FALSE, recursive=TRUE)
for (s in SUPP) {
  e <- new.env()
  tryCatch(sys.source(s$script, envir=e),
           error=function(err) message("  (", s$script, " skipped: ", conditionMessage(err), ")"))
  for (ext in c("png","svg","pdf")) {
    f <- sprintf("%s.%s", s$out, ext)
    if (file.exists(f)) file.copy(f, sprintf("%s/%s.%s", sdir, s$dst, ext), overwrite=TRUE)
  }
  png <- sprintf("%s.png", s$out)
  if (file.exists(png) && requireNamespace("magick", quietly=TRUE))
    magick::image_write(magick::image_convert(magick::image_read(png), "tiff"),
                        path=sprintf("%s/%s.tiff", tifdir, s$dst), format="tiff", density="300x300")
}
## S3 = somatic mutation oncoplot by risk group (pre-rendered by fig_tmb.R, png only)
if (file.exists("results/Fig_oncoplot.png")) {
  file.copy("results/Fig_oncoplot.png", sprintf("%s/Supplementary_Figure_S3.png", sdir), overwrite=TRUE)
  if (requireNamespace("magick", quietly=TRUE))
    magick::image_write(magick::image_convert(magick::image_read("results/Fig_oncoplot.png"), "tiff"),
                        path=sprintf("%s/Supplementary_Figure_S3.tiff", tifdir), format="tiff", density="300x300")
}
cat("published 3 supplementary figures -> submission_BMC/supplementary\n")
