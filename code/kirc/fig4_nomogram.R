## Fig 4 列线图 + 校准 + 独立预后(风险评分 vs 临床)
suppressPackageStartupMessages({ library(data.table); library(rms); library(survival) })
P <- readRDS("data/raw/kirc_prognosis.rds"); df <- P$df
cl <- fread("data/raw/KIRC_clinical.tsv"); cl$sample <- cl$sampleID
d <- merge(df[, c("sample","risk","OS","OS.time")],
           cl[, .(sample, age=age_at_initial_pathologic_diagnosis, gender,
                  grade=neoplasm_histologic_grade, stage=pathologic_stage)], by="sample")
## 清洗
d$Stage <- as.numeric(factor(gsub("Stage ","",d$stage), levels=c("I","II","III","IV")))
d$Grade <- suppressWarnings(as.numeric(gsub("G","", d$grade)))
d$Age   <- suppressWarnings(as.numeric(d$age)); d$Risk <- d$risk
d <- d[complete.cases(d[, c("OS","OS.time","Risk","Stage","Grade","Age")]) & d$OS.time > 0, ]
d <- as.data.frame(d); cat("纳入:", nrow(d), "例\n")

## 独立预后:多因素 Cox(风险评分是否独立于分期/分级/年龄)
mcox <- coxph(Surv(OS.time, OS) ~ Risk + Stage + Grade + Age, data = d); sm <- summary(mcox)
cat("多因素 Cox(独立预后):\n"); print(round(sm$coefficients[, c(2,5)], 4))

## rms 列线图
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(OS.time, OS) ~ Risk + Stage + Grade + Age, data = d, x=TRUE, y=TRUE, surv=TRUE)
surv <- Survival(f)
nom <- nomogram(f, fun = list(function(x) surv(365,x), function(x) surv(1095,x), function(x) surv(1825,x)),
                funlabel = c("1-year OS","3-year OS","5-year OS"), lp = FALSE)

pdf("results/Fig4_nomogram.pdf", width = 12, height = 11)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), heights = c(1.2, 1))
par(mar = c(4,4,3,2))
plot(nom, xfrac = .35); title("A  Nomogram (risk score + stage + grade + age)")

## 校准 1/3/5 年
cols <- c("#E64B35","#4DBBD5","#00A087"); us <- c(365,1095,1825); labs <- c("1-yr","3-yr","5-yr")
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Nomogram-predicted OS",ylab="Observed OS",main="B  Calibration")
abline(0,1,lty=2,col="grey")
for (i in seq_along(us)) {
  fi <- cph(Surv(OS.time,OS)~Risk+Stage+Grade+Age, data=d, x=TRUE,y=TRUE,surv=TRUE,time.inc=us[i])
  ca <- tryCatch(calibrate(fi, u=us[i], cmethod="KM", m=60, B=150), error=function(e) NULL)
  if(!is.null(ca)) lines(ca[,"mean.predicted"], ca[,"KM"], type="b", col=cols[i], pch=19, lwd=2)
}
legend("bottomright", labs, col=cols, lwd=2, pch=19, bty="n")

## 独立预后森林
hr <- sm$conf.int[, c(1,3,4)]; pv <- sm$coefficients[,5]; rn <- rownames(hr)
par(mar=c(4,6,3,2)); plot(0,0,type="n",xlim=range(c(0.5,hr),na.rm=TRUE),ylim=c(.5,length(rn)+.5),
     yaxt="n",xlab="Hazard ratio (95% CI)",ylab="",main="C  Multivariate Cox",log="x")
abline(v=1,lty=2,col="grey"); axis(2,at=seq_along(rn),labels=rev(rn),las=1)
for(i in seq_along(rn)){ j<-length(rn)-i+1
  segments(hr[i,2],j,hr[i,3],j,lwd=2,col="#3C5488"); points(hr[i,1],j,pch=19,col="#E64B35")
  text(max(hr,na.rm=TRUE),j,sprintf("p=%.1e",pv[i]),pos=4,xpd=NA,cex=.8) }
dev.off()
## 列线图单独全图(layout 里 plot.nomogram 易空白)
png("results/Fig4A_nomogram.png", width=1400, height=760, res=120)
par(mar=c(4,4,3,2)); plot(nom, xfrac=.32); title("A  Nomogram (risk score + stage + grade + age)"); dev.off()
## B 校准 + C 独立预后森林
png("results/Fig4BC.png", width=1400, height=620, res=120); par(mfrow=c(1,2), mar=c(4,4,3,2))
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Nomogram-predicted OS",ylab="Observed OS",main="B  Calibration"); abline(0,1,lty=2,col="grey")
for(i in seq_along(us)){fi<-cph(Surv(OS.time,OS)~Risk+Stage+Grade+Age,data=d,x=TRUE,y=TRUE,surv=TRUE,time.inc=us[i]);ca<-tryCatch(calibrate(fi,u=us[i],cmethod="KM",m=60,B=150),error=function(e)NULL);if(!is.null(ca))lines(ca[,"mean.predicted"],ca[,"KM"],type="b",col=cols[i],pch=19,lwd=2)}
legend("bottomright",labs,col=cols,lwd=2,pch=19,bty="n")
par(mar=c(4,6,3,2)); plot(0,0,type="n",xlim=range(c(0.5,hr),na.rm=TRUE),ylim=c(.5,length(rn)+.5),yaxt="n",xlab="HR (95% CI)",ylab="",main="C  Multivariate Cox (independent)",log="x")
abline(v=1,lty=2,col="grey"); axis(2,at=seq_along(rn),labels=rev(rn),las=1)
for(i in seq_along(rn)){j<-length(rn)-i+1;segments(hr[i,2],j,hr[i,3],j,lwd=2,col="#3C5488");points(hr[i,1],j,pch=19,col="#E64B35")}
dev.off()
cat("Fig4 完成 | 风险评分多因素 Cox HR=", round(sm$conf.int["Risk",1],2), " p=", signif(sm$coefficients["Risk",5],2), "\n")
