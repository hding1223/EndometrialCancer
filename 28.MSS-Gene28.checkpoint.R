#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")
#install.packages("ggplot2")
#install.packages("ggpubr")

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(tidyverse)

pFilter=0.001             
   
expFile='TCGA_UCEC_TPM.txt'
geneFile="IGC-gene.txt"       
setwd("D:/workspace/TCGA/UCEC/input")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

# Immune-infiltration analysis step.
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[c(sameGene),])
data=log2(data+1)


group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)

sameSample=intersect(row.names(rs), row.names(data))
rs=rs[sameSample,,drop=F]
data=data[sameSample,,drop=F]
data=cbind(data, rs)

x=as.numeric(data[,'riskscore'])
outTab=data.frame()
for(i in sameGene){
  y=as.numeric(data[,i])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<pFilter){
		outTab=rbind(outTab, cbind(Query="riskscore", Gene=i, cor, pvalue))
	}
}

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/TME")
write.table(file="IRG-corResult.txt", outTab, sep="\t", quote=F, row.names=F)
dataM=data[c("riskscore", as.vector(outTab[,2]))]
M=cor(dataM)
view(M)
pdf(file="IRG-corpot1.pdf",width=7,height=7)
corrplot(M,
         method = "circle",
         order = "original",
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
         )
dev.off()

pdf(file="IRG-corpot2.pdf",width=8,height=8)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7,
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
