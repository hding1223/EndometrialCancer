library("survival")
library("survminer")
library(tidyverse)
#gene = "GAL"
#gene = "PGR"
gene = "THRB"
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
scoreFile <- "final_mode_riskscore.txt"
score <- read.table(scoreFile, header=FALSE, sep="\t", check.names=FALSE, row.names=1)
dimnames = list(rownames(score),colnames(score))
score_data = matrix(as.numeric(as.matrix(score)),nrow=nrow(score),dimnames=dimnames)

# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
cliFile="clinical_above30_ml.txt"
clinical <- read.table(cliFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)

# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(score_data), rownames(clinical))
clinical <- clinical[sameSample, ]
clinical$OS.time = clinical$OS.time/365

# Analysis step; see the surrounding code for details.
setwd("D:/workspace/TCGA/UCEC/input") 
data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)

# Analysis step; see the surrounding code for details.
data=t(data)
#view(data)
rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))

# Analysis step; see the surrounding code for details.
sameSample=intersect(row.names(data),row.names(clinical))
data=data[sameSample,]
clinical = clinical[sameSample,]
rt = cbind(clinical,data)
group = ifelse(rt[,gene]>quantile(rt[,gene],seq(0,1,1/2))[2],"High","Low")
rt[,"group"] = group
length = length(levels(factor(group)))

# Survival/Cox analysis step.
diff=survdiff(Surv(OS.time,OS) ~group,data=rt)
pValue = 1-pchisq(diff$chisq,df=length-1)
pValue = paste0("p=",sprintf("%.04f",pValue))
fit <- survfit(Surv(OS.time,OS) ~group,data=rt)

surPlot = ggsurvplot(fit,
                      data=rt,
                      conf.int=F,
                      pval=pValue,
                      pval.size=6,
                      surv.median.line="hv",
                      legend.labs=c("High level", "Low level"),
                      xlab="Time(years)",
                      ylab="Overall survival",
                      break.time.by = 1,
                      palette=c("red", "blue"),
                      risk.table=F,
                      risk.table.title="",
                      risk.table.height=.25,
                      legend.title=gene
                      )
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/KS-Genes")
pdf(file=paste0("ks-survival-",gene,".pdf"),width=5.5,height=5,onefile=FALSE)
print(surPlot)
dev.off()
