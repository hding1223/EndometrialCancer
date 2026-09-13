
library("survival")
library("survminer")
library(tidyverse)

setwd("D:/workspace/TCGA/UCEC/input")


data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)


# Process the gene data.
moduleGrey=read.table('append/module_grey.txt', header=F, sep="\t", check.names=F, stringsAsFactors=F)[,1]
degs=read.table('append/DEGs_362.txt', header=F, sep="\t", check.names=F, stringsAsFactors=F)[,1]
gene=data.frame(intersect(moduleGrey,degs))
# Process the gene data.
sameGene=intersect(row.names(data), as.vector(gene[,1]))
# Process the gene data.
#view(sameGene)
data=data[sameGene,]
data=t(data)
data=log2(data+1)
view(data)

rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))
#view(data)
clinical = read.table('TCGA_UCEC_clinical_selected_above_32_MMS.txt',header=T, sep="\t",check.names=F,row.names = 1)
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
clinical = clinical[, c("OS","OS.time")]
#clinical$time = clinical$time/365
#view(clinical)

# Process the sample data.
sameSample=intersect(row.names(data),row.names(clinical))
# WGCNA processing step.
trainSamples=read.table('append/trainSamples.txt', header=F, sep="\t", check.names=F, stringsAsFactors=F)[,1]
sameSample=intersect(sameSample,trainSamples)
data=data[sameSample,]
clinical = clinical[sameSample,]
rt = cbind(clinical,data)

#0.05,0.01,0.001
p.value=0.05
outTab = data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  # Cox regression analysis step.
  cox <- coxph(Surv(OS.time,OS) ~ rt[,i],data=rt)
  coxSummary = summary(cox)
  coxP = coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP < p.value){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"]
                 )
    )
  }
}
setwd("D:/workspace/TCGA/UCEC/output/append/0.2.Immune_Related_Genes")

write.table(outTab,'IRG_Gene-OS-uniCox1.txt',sep="\t",quote=F,row.names = F)
rt <- read.table('IRG_Gene-OS-uniCox1.txt',header = T, sep="\t",check.names = F,row.names = 1)
rt <- rt[rt$pvalue<0.05,]
#view(rt)

gene <- rownames(rt)
hr <-sprintf("%.3f",rt$"HR")
hrLow <-sprintf("%.3f",rt$"HR.95L")
hrHigh <-sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.05,"<0.05",sprintf("%.3f",rt$pvalue))

# Plot the corresponding figure.
pdf(file="Gene-IRG-uniCoxforest.pdf",width =7,height=nrow(rt)/13+5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)

par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,"Hazard ratio",cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim=c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col='black',lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) >1,'red','blue')
points(as.numeric(hr),n:1,pch=15,col=boxcolor,cex=1.6)
axis(1)
dev.off()



