#install.packages("regplot")

library("survival")
library(survminer)
library(rms)
library(nomogramFormula)

setwd("D:/workspace/TCGA/UCEC/input/clinical")
cliFile="clinical_all_items.txt"
cli = read.table(cliFile,header=T, sep="\t",check.names=F,row.names = 1)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
#setwd("D:/workspace/TCGA/UCEC/input")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)

cli$age = as.numeric(cli$age)

sameSample = intersect(row.names(cli),row.names(rs))
riskscore = rs[sameSample,,drop=F]
#colnames(riskscore) = 'riskscore'
clinical =  cli[sameSample,,drop=F]
rt = cbind(riskscore,clinical)

index <- sort(sample(nrow(rt),nrow(rt)*0.7))
train <- rt[index,]
test <- rt[-index,]
res=coxph(Surv(OS.time, OS) ~riskscore,data=train)
#res=coxph(Surv(OS.time, OS) ~age+riskscore+stage+grade,data=train)
ddist <- datadist(train)
options(datadist="ddist")

f_cph <- cph(Surv(OS.time, OS) ~riskscore,x=T,y=T,surv=T,data=train)
#f_cph <- cph(Surv(OS.time, OS) ~age+riskscore+stage+grade,x=T,y=T,surv=T,data=train)
summary(f_cph)
ddist <- datadist(train)
options(datadist="ddist")
med <- Quantile(f_cph)
surv <- Survival(f_cph)

#nomo <- nomogram(f_cph,fun=list(function(x) surv(365,x),
#                                function(x) surv(365*3,x),
#                                function(x) surv(365*5,x)),
#                      funlabel = c("1-year Survival Probability",
#                                   "3-year Survival Probability",
#                                  "5-year Survival Probability"))
png("mss-nomogram.png", width=640, height=480)
nomo <- nomogram(f_cph,fun=list(function(x) surv(365,x),
                              function(x) surv(365*3,x),
                             function(x) surv(365*5,x)),
                 funlabel = c("1-year Survival Probability",
                              "3-year Survival Probability",
                              "5-year Survival Probability"))
plot(nomo)
dev.off()
