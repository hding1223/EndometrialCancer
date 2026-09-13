#install.packages("regplot")

library("survival")
library(regplot)
library(rms)

setwd("D:/workspace/TCGA/UCEC/input/clinical") 
#cliFile="clinical_above30_ml.txt"
cliFile="clinical_all_items.txt"
cli = read.table(cliFile,header=T, sep="\t",check.names=F,row.names = 1)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")      
rs = read.table('final_mode_riskscore.txt',header=F, sep="\t",check.names=F,row.names = 1)

cli$age = as.numeric(cli$age)
cli$OS.time = as.numeric(cli$OS.time/365)

sameSample = intersect(row.names(cli),row.names(rs))
riskscore = rs[sameSample,,drop=F]
clinical =  cli[sameSample,,drop=F]
colnames(riskscore) = 'riskscore'

rt = cbind(riskscore,clinical)
res.cox=coxph(Surv(OS.time, OS) ~ riskscore,data=rt)

nom1 = regplot(res.cox,
               clickable = F,
               title="",
               points=TRUE,
               droplines=TRUE,
               observation=rt[50,],
               rank="sd",
               failtime=c(1,3,5),
               prfail=F)