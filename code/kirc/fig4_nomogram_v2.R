## Fig 4 升级版: 列线图 + 校准 + 独立预后森林(红蓝配色 + Arial)
source("code/kirc/_style.R")
suppressPackageStartupMessages({ library(data.table); library(rms); library(survival) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
cl <- fread("data/raw/KIRC_clinical.tsv"); cl$sample <- cl$sampleID
d <- merge(df[,c("sample","risk","OS","OS.time")],
           cl[,.(sample, age=age_at_initial_pathologic_diagnosis, grade=neoplasm_histologic_grade, stage=pathologic_stage)], by="sample")
d$Stage <- as.numeric(factor(gsub("Stage ","",d$stage), levels=c("I","II","III","IV")))
d$Grade <- suppressWarnings(as.numeric(gsub("G","",d$grade)))
d$Age <- suppressWarnings(as.numeric(d$age)); d$Risk <- d$risk
d <- as.data.frame(d[complete.cases(d[,c("OS","OS.time","Risk","Stage","Grade","Age")]) & d$OS.time>0,])
dd <- datadist(d); options(datadist="dd")
f <- cph(Surv(OS.time,OS)~Risk+Stage+Grade+Age, data=d, x=TRUE,y=TRUE,surv=TRUE)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365,x), function(x) surv(1095,x), function(x) surv(1825,x)),
                funlabel=c("1-year OS","3-year OS","5-year OS"), lp=FALSE)
sm <- summary(coxph(Surv(OS.time,OS)~Risk+Stage+Grade+Age, d))
hr <- sm$conf.int[,c(1,3,4)]; rn <- rownames(hr)
us <- c(365,1095,1825); cols <- c(PAL$high, PAL$accent3, PAL$low); labs <- c("1-yr","3-yr","5-yr")

dev_base <- function(file, w_mm, h_mm, draw_fn){
  w<-w_mm/25.4; h<-h_mm/25.4
  for(d3 in list(function() svglite::svglite(paste0(file,".svg"),width=w,height=h),
                 function() ragg::agg_png(paste0(file,".png"),width=w,height=h,units="in",res=300)))
    tryCatch({ d3(); par(family="Arial"); draw_fn(); dev.off() }, error=function(e){ if(!is.null(dev.list())) dev.off() })
}
dev_base("results/Fig4A_nomogram_v2", 180, 95, function(){ par(mar=c(4,4,3,2)); plot(nom, xfrac=.32, col.grid=PAL$n_light); title("Nomogram (risk + stage + grade + age)", cex.main=.9) })
dev_base("results/Fig4BC_v2", 180, 80, function(){
  par(mfrow=c(1,2), mar=c(4,4,3,2))
  plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Nomogram-predicted OS",ylab="Observed OS",main="Calibration"); abline(0,1,lty=2,col=PAL$n_mid)
  for(i in seq_along(us)){ fi<-cph(Surv(OS.time,OS)~Risk+Stage+Grade+Age,data=d,x=TRUE,y=TRUE,surv=TRUE,time.inc=us[i])
    ca<-tryCatch(calibrate(fi,u=us[i],cmethod="KM",m=60,B=120),error=function(e)NULL)
    if(!is.null(ca)) lines(ca[,"mean.predicted"],ca[,"KM"],type="b",col=cols[i],pch=19,lwd=2) }
  legend("bottomright",labs,col=cols,lwd=2,pch=19,bty="n")
  par(mar=c(4,6,3,2)); plot(0,0,type="n",xlim=range(c(.5,hr),na.rm=TRUE),ylim=c(.5,length(rn)+.5),yaxt="n",
    xlab="Hazard ratio (95% CI)",ylab="",main="Multivariate Cox",log="x"); abline(v=1,lty=2,col=PAL$n_mid)
  axis(2,at=seq_along(rn),labels=rev(rn),las=1)
  for(i in seq_along(rn)){ j<-length(rn)-i+1; segments(hr[i,2],j,hr[i,3],j,lwd=2.5,col=PAL$low); points(hr[i,1],j,pch=19,col=PAL$high,cex=1.1) }
})
cat("Fig4 升级版完成 (Fig4A_nomogram_v2 + Fig4BC_v2)\n")
