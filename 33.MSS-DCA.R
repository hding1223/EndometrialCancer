#remotes::install_github('yikeshu0611/ggDCA') 

library(survivalROC)
library(pROC)
library(timeROC)
library(nomogramFormula)
library("survival")
library(regplot)
library(rms)
library(ggDCA)

setwd("D:/workspace/TCGA/UCEC/input/clinical")
cli = read.table('clinical_all_items.txt',header=T, sep="\t",check.names=F,row.names = 1)
cli$OS.time = as.numeric(cli$OS.time)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)

train_case = read.table('Train-Case.txt',header=T, sep="\t",check.names=F,row.names = 1)
validation_case = read.table('Validation-Case.txt',header=T, sep="\t",check.names=F,row.names = 1)


#cli$age = as.numeric(cli$age)
cli$OS.time = as.numeric(cli$OS.time)

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

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/verify")
pdf("IRG_train_DCA-1.pdf",width =8,height=6)
f_cph <-cph(Surv(OS.time, OS) ~ riskscore,data=train_rt,x=T,y=T,surv=T)
dca_training <- dca(f_cph,times=c(365,3*365,5*365))
ggplot(dca_training)
dev.off()

pdf("IRG_train_DCA-3.pdf",width =8,height=6)
f_cph <-cph(Surv(OS.time, OS) ~ riskscore,data=validation_rt,x=T,y=T,surv=T)
dca_training <- dca(f_cph,times=c(365,3*365,5*365))
ggplot(dca_training)
dev.off()

pdf("IRG_train_DCA-5.pdf",width =8,height=6)
f_cph <-cph(Surv(OS.time, OS) ~ riskscore,data=all_rt,x=T,y=T,surv=T)
dca_training <- dca(f_cph,times=c(365,3*365,5*365))
ggplot(dca_training)
dev.off()