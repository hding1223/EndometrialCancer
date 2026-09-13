
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#???ð?
library(limma)
library(estimate)

expFile="TCGA_UCEC_TPM.txt"       #?????????ļ?
setwd("D:/workspace/TCGA/UCEC/input")       #???ù???Ŀ¼

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

out=rbind(ID=colnames(data),data)
write.table(out,file="IRG_uniq.symbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="IRG_uniq.symbol.txt", 
                  output.f="IRG_commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "IRG_commonGenes.gct",
              output.ds="IRG_estimateScore.gct")

scores=read.table("IRG_estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="IRG_TMEscores.txt", sep="\t", quote=F, col.names=F)


