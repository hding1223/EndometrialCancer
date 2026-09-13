#install.packages("regplot")

library("survival")
library(survminer)
library(rms)
library(nomogramFormula)

setwd("D:/workspace/TCGA/UCEC/input")

data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)
data = data[rowMeans(data)>1,]

geneFile = 'MSS_ML_SELECTED_GENE.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# Immune-infiltration analysis step.
sameGene=intersect(row.names(data), as.vector(gene[,1]))

data=data[sameGene,]
# Analysis step; see the surrounding code for details.
data=t(data)

rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))

setwd("D:/workspace/TCGA/UCEC/input/clinical")
#cliFile="clinical_above30_ml.txt"
cliFile="clinical_all_items.txt"
cli = read.table(cliFile,header=T, sep="\t",check.names=F,row.names = 1)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")  
rs = read.table('final_mode_riskscore.txt',header=F, sep="\t",check.names=F,row.names = 1)
cli$age = as.numeric(cli$age)

sameSample = intersect(row.names(cli),row.names(rs))
sameSample = intersect(sameSample,row.names(data))
riskscore = rs[sameSample,,drop=F]
clinical =  cli[sameSample,,drop=F]
data = data[sameSample,,drop=F]
colnames(riskscore) = 'riskscore'

rt = cbind(riskscore,clinical)
rt = cbind(rt,data)

index <- sort(sample(nrow(rt),nrow(rt)*0.7))
train <- rt[index,]
test <- rt[-index,]
res=coxph(Surv(OS.time, OS) ~riskscore,data=train)
ddist <- datadist(train)
options(datadist="ddist")

f_cph <- cph(Surv(OS.time, OS) ~GAL+PGR+THRB,x=T,y=T,surv=T,data=train)
#f_cph <- cph(Surv(OS.time, OS) ~riskscore,x=T,y=T,surv=T,data=train)
summary(f_cph)
ddist <- datadist(train)
options(datadist="ddist")
med <- Quantile(f_cph)
surv <- Survival(f_cph)

#nomo <- nomogram(f_cph,fun=list(function(x) surv(365,x),
#                                function(x) surv(365*3,x),
#                                function(x) surv(365*5,x)),
#                 funlabel = c("1-year Survival Probability",
#                              "3-year Survival Probability",
#                              "5-year Survival Probability"))

plot(nomogram(f_cph,fun=list(function(x) surv(365,x),
                             function(x) surv(365*3,x),
                             function(x) surv(365*5,x)),
              funlabel = c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability"))
)