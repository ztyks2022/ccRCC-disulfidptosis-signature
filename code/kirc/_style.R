## _style.R — KIRC 全套图统一 Nature 审美(配色/主题/导出)。各图 source() 此文件。
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid); library(scales) })

## 语义调色板:中性灰骨架 + 红(高危/上调)+ 蓝(低危/下调)+ 少量强调色
PAL <- list(high="#B2182B", high_soft="#E6A6A1", low="#2166AC", low_soft="#A6C5E0",
            n_dark="#272727", n_mid="#767676", n_light="#D9D9D9",
            accent="#1B7837", accent2="#762A83", accent3="#E08214")
risk_cols <- c(Low = PAL$low, High = PAL$high)             # 风险高/低统一:低蓝高红
two_cols  <- c(PAL$low, PAL$high)                          # 两组对比(正常/肿瘤、C1/C2…)
div_ramp  <- function() circlize::colorRamp2(c(-2, 0, 2), c(PAL$low, "white", PAL$high))
seq_reds  <- function() circlize::colorRamp2(c(0, .5, 1), c("white", PAL$high_soft, PAL$high))

## 统一主题:theme_classic + Arial 7pt + 细轴线 + 无网格 + 无框图例
theme_nat <- function(base = 7) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(axis.line   = element_line(linewidth = .35, colour = "black"),
          axis.ticks  = element_line(linewidth = .35, colour = "black"),
          axis.text   = element_text(colour = "black", size = base - 1),
          axis.title  = element_text(size = base),
          legend.key.size = unit(3.2, "mm"),
          legend.title = element_text(size = base - 1),
          legend.text  = element_text(size = base - 1.3),
          legend.background = element_blank(), legend.key = element_blank(),
          plot.title  = element_text(size = base, face = "bold", hjust = 0),
          plot.tag    = element_text(size = base + 2, face = "bold"),
          strip.background = element_blank(),
          strip.text  = element_text(size = base - .5, face = "bold"),
          panel.grid  = element_blank())
}
theme_set(theme_nat())

## 三格式导出(矢量优先):SVG(可编辑) + PDF + 高分 PNG。单设备失败不影响其他。
.dev_save <- function(open_fn, p, tag) {
  open_fn(); tryCatch(print(p), error = function(e) message("  (", tag, " 跳过: ", conditionMessage(e), ")"))
  if (!is.null(grDevices::dev.list())) grDevices::dev.off()
}
save_pub <- function(p, file, w_mm = 183, h_mm = 120, dpi = 300) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  .dev_save(function() svglite::svglite(paste0(file, ".svg"), width = w, height = h), p, "SVG")
  .dev_save(function() grDevices::cairo_pdf(paste0(file, ".pdf"), width = w, height = h, family = "Arial"), p, "PDF")
  .dev_save(function() ragg::agg_png(paste0(file, ".png"), width = w, height = h, units = "in", res = dpi), p, "PNG")
}
add_tags <- function(fig) fig + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))
