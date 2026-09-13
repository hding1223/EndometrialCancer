#install.packages("regplot")

library(survivalROC)
library(pROC)
library(timeROC)
library(nomogramFormula)
library("survival")
library(regplot)
library(rms)

setwd("D:/workspace/TCGA/UCEC/input/clinical")
cli = read.table('clinical_all_items.txt',header=T, sep="\t",check.names=F,row.names = 1)
#cli$age = as.numeric(cli$age)
cli$OS.time = as.numeric(cli$OS.time)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)

train_case = read.table('Train-Case.txt',header=T, sep="\t",check.names=F,row.names = 1)
validation_case = read.table('Validation-Case.txt',header=T, sep="\t",check.names=F,row.names = 1)



all_sameSample = intersect(row.names(cli),row.names(rs))
all_riskscore = rs[all_sameSample,,drop=F]
all_clinical =  cli[all_sameSample,,drop=F]
all_rt = cbind(all_riskscore,all_clinical)

train_sameSample = intersect(row.names(cli),row.names(rs))
train_sameSample = intersect(train_sameSample,row.names(train_case))
train_riskscore = rs[train_sameSample,,drop=F]
train_clinical =  cli[train_sameSample,,drop=F]

validation_sameSample = intersect(row.names(cli),row.names(rs))
validation_sameSample = intersect(validation_sameSample,row.names(validation_case))
validation_riskscore = rs[validation_sameSample,,drop=F]
validation_clinical =  cli[validation_sameSample,,drop=F]

train_rt = cbind(train_riskscore,train_clinical)
validation_rt = cbind(validation_riskscore,validation_clinical)

# Analysis step; see the surrounding code for details.
f5<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=train_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365)
cal5 <- calibrate(f5,cmethod="KM",method="boot",u=365,m=100,B=600)

f8<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=train_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*3)
cal8 <- calibrate(f8,cmethod="KM",method="boot",u=365*3,m=100,B=600)

f10<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=train_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*5)
cal10 <- calibrate(f10,cmethod="KM",method="boot",u=365*5,m=100,B=600)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/verify")
pdf("IRG_train_calibration.pdf",width =8,height=8)
plot(cal5,lwd=2,lty=0,errbar.col=c("#2166AC"),
     xlim = c(0.5,1),ylim=c(0.5,1),
     xlab="Prediced OS (%)",ylab = "Observed OS(%)",
     col=c("#2166AC"),
     cex.lab=1.2,cex.axis=1,cex.main=1.2,cex.sub=0.6)
lines(cal5[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,col=c("#2166AC"),pch=16)
mtext("")

plot(cal8,
     lwd=2,
     lty=0,
     errbar.col=c("#B21828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#B21828"),add = T)

lines(cal8[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#B21828"))
mtext("")

plot(cal10,
     lwd=2,
     lty=0,
     errbar.col=c("#441828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#441828"),add = T)

lines(cal10[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#441828"))
mtext("")

box(lwd=1)
abline(0,1,lty=3,
       lwd=2,
       col=c("#224444")
       )


legend("topleft",legend=c("1-year","3-year","5-year"),
       col=c("#2166AC","#B21828","#441828"),
       lwd=2,
       cex=1.2,
       bty="n")
dev.off()

# Analysis step; see the surrounding code for details.
f5<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365)
cal5 <- calibrate(f5,cmethod="KM",method="boot",u=365,m=50,B=600)

f8<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*3)
cal8 <- calibrate(f8,cmethod="KM",method="boot",u=365*3,m=50,B=600)

f10<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*5)
cal10 <- calibrate(f10,cmethod="KM",method="boot",u=365*5,m=50,B=600)


pdf("IRG_validation_calibration.pdf",width =8,height=8)
plot(cal5,lwd=2,lty=0,errbar.col=c("#2166AC"),
     xlim = c(0.5,1),ylim=c(0.5,1),
     xlab="Prediced OS (%)",ylab = "Observed OS(%)",
     col=c("#2166AC"),
     cex.lab=1.2,cex.axis=1,cex.main=1.2,cex.sub=0.6)
lines(cal5[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,col=c("#2166AC"),pch=16)
mtext("")

plot(cal8,
     lwd=2,
     lty=0,
     errbar.col=c("#B21828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#B21828"),add = T)

lines(cal8[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#B21828"))
mtext("")

plot(cal10,
     lwd=2,
     lty=0,
     errbar.col=c("#441828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#441828"),add = T)

lines(cal10[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#441828"))
mtext("")

box(lwd=1)
abline(0,1,lty=3,
       lwd=2,
       col=c("#224444")
)


legend("topleft",legend=c("1-year","3-year","5-year"),
       col=c("#2166AC","#B21828","#441828"),
       lwd=2,
       cex=1.2,
       bty="n")
dev.off()

# Analysis step; see the surrounding code for details.
f5<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365)
cal5 <- calibrate(f5,cmethod="KM",method="boot",u=365,m=100,B=600)

f8<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*3)
cal8 <- calibrate(f8,cmethod="KM",method="boot",u=365*3,m=100,B=600)

f10<- cph(formula = Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T,na.action=na.delete,time.inc=365*5)
cal10 <- calibrate(f10,cmethod="KM",method="boot",u=365*5,m=100,B=600)


pdf("IRG_all_calibration.pdf",width =8,height=8)
plot(cal5,lwd=2,lty=0,errbar.col=c("#2166AC"),
     xlim = c(0.5,1),ylim=c(0.5,1),
     xlab="Prediced OS (%)",ylab = "Observed OS(%)",
     col=c("#2166AC"),
     cex.lab=1.2,cex.axis=1,cex.main=1.2,cex.sub=0.6)
lines(cal5[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,col=c("#2166AC"),pch=16)
mtext("")

plot(cal8,
     lwd=2,
     lty=0,
     errbar.col=c("#B21828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#B21828"),add = T)

lines(cal8[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#B21828"))
mtext("")

plot(cal10,
     lwd=2,
     lty=0,
     errbar.col=c("#441828"),
     xlim = c(0.5,1),ylim=c(0.5,1),col=c("#441828"),add = T)

lines(cal10[,c('mean.predicted','KM')],
      type = 'b',
      lwd =1,
      pch=16,
      col=c("#441828"))
mtext("")

box(lwd=1)
abline(0,1,lty=3,
       lwd=2,
       col=c("#224444")
)


legend("topleft",legend=c("1-year","3-year","5-year"),
       col=c("#2166AC","#B21828","#441828"),
       lwd=2,
       cex=1.2,
       bty="n")
dev.off()
