## DCA (decision curve analysis) — 5-year OS net benefit
## Nomogram (Risk+Stage+Grade+Age) vs Risk-score-alone vs Clinical (Stage+Grade+Age)
## Time-to-event net benefit (Vickers), KM-based. Produces ggplot object pDCA.
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(rms); library(survival); library(ggplot2) })

P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
cl <- fread("data/raw/KIRC_clinical.tsv"); cl$sample <- cl$sampleID
d <- merge(df[,c("sample","risk","OS","OS.time")],
           cl[,.(sample, age=age_at_initial_pathologic_diagnosis, grade=neoplasm_histologic_grade, stage=pathologic_stage)], by="sample")
d$Stage <- as.numeric(factor(gsub("Stage ","",d$stage), levels=c("I","II","III","IV")))
d$Grade <- suppressWarnings(as.numeric(gsub("G","",d$grade)))
d$Age <- suppressWarnings(as.numeric(d$age)); d$Risk <- d$risk
d <- as.data.frame(d[complete.cases(d[,c("OS","OS.time","Risk","Stage","Grade","Age")]) & d$OS.time>0,])
dd <- datadist(d); options(datadist="dd")

t_horizon <- 1825  # 5-year
predict_event <- function(form){
  f <- cph(as.formula(form), data=d, x=TRUE, y=TRUE, surv=TRUE)
  s <- survest(f, newdata=d, times=t_horizon)$surv
  1 - as.numeric(s)                       # predicted event prob by 5y
}
risk_models <- list(
  `Nomogram (all)` = predict_event("Surv(OS.time,OS)~Risk+Stage+Grade+Age"),
  `Risk score`     = predict_event("Surv(OS.time,OS)~Risk"),
  `Clinical (stage+grade+age)` = predict_event("Surv(OS.time,OS)~Stage+Grade+Age")
)

km_event <- function(idx){                # KM event prob by t among a subset
  if(sum(idx) < 5) return(NA_real_)
  sf <- survfit(Surv(OS.time,OS)~1, data=d[idx,,drop=FALSE])
  ss <- summary(sf, times=t_horizon, extend=TRUE)$surv
  1 - ss[length(ss)]
}
n <- nrow(d); pt_grid <- seq(0.02, 0.60, by=0.01)
ev_all <- km_event(rep(TRUE,n))

nb_model <- function(risk){
  sapply(pt_grid, function(pt){
    flag <- risk >= pt
    nf <- sum(flag); if(nf==0) return(0)
    pev <- km_event(flag); if(is.na(pev)) return(0)
    (nf/n) * (pev - (1-pev)*(pt/(1-pt)))
  })
}
nb_all <- sapply(pt_grid, function(pt) ev_all - (1-ev_all)*(pt/(1-pt)))

dca <- rbind(
  data.frame(pt=pt_grid, nb=nb_all, model="Treat all"),
  data.frame(pt=pt_grid, nb=0,       model="Treat none"),
  do.call(rbind, lapply(names(risk_models), function(m)
    data.frame(pt=pt_grid, nb=nb_model(risk_models[[m]]), model=m)))
)
dca$model <- factor(dca$model, levels=c("Nomogram (all)","Risk score","Clinical (stage+grade+age)","Treat all","Treat none"))
pal_dca <- c("Nomogram (all)"=PAL$high, "Risk score"=PAL$accent3,
             "Clinical (stage+grade+age)"=PAL$low, "Treat all"="grey55", "Treat none"="grey20")
lty_dca <- c("Nomogram (all)"=1,"Risk score"=1,"Clinical (stage+grade+age)"=1,"Treat all"=2,"Treat none"=2)

pDCA <- ggplot(dca, aes(pt, nb, colour=model, linetype=model)) +
  geom_line(linewidth=.6) +
  scale_colour_manual(values=pal_dca, name=NULL) +
  scale_linetype_manual(values=lty_dca, name=NULL) +
  coord_cartesian(ylim=c(-0.02, max(nb_all)*1.05), xlim=c(0,0.6)) +
  labs(title="Decision curve analysis (5-year OS)",
       x="Threshold probability", y="Net benefit") +
  theme_nat() + theme(legend.position=c(.72,.8), legend.key.height=unit(3.2,"mm"),
                      legend.text=element_text(size=6))

if (!exists("ASSEMBLE_MODE")) {
  save_pub(pDCA, "results/Fig4_DCA", w_mm=95, h_mm=80)
  cat("DCA done | nomogram NB@pt=0.2:",
      round(nb_model(risk_models[["Nomogram (all)"]])[which(abs(pt_grid-0.2)<1e-6)],4), "\n")
}
