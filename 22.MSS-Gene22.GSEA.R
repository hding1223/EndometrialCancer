#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#packageVersion(clusterProfiler)
#packageVersion("clusterProfiler")

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)

library(enrichplot)
library(tidyverse)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)


expFile="TCGA_UCEC_TPM.txt"     
#gmtFile="c5.go.v2024.1.Hs.symbols.gmt"
gmtFile="c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"
setwd("D:/workspace/TCGA/UCEC/input")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
rownames(data)=gsub('[.]', '-', rownames(data))

highGroup = subset(rs,rs[,"riskscore"]>median(rs[,"riskscore"]))
lowGroup = subset(rs,rs[,"riskscore"]<=median(rs[,"riskscore"]))

dataL=data[rownames(lowGroup),]
dataH=data[rownames(highGroup),]
#view(dataH)
dataL = t(dataL)
dataH = t(dataH)
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
#view(meanL)
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)
#view(logFC)

gmt=read.gmt(gmtFile)

kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/GSEA") 
write.table(kkTab,file="IRG_GSEA_KEGG.result.txt",sep="\t",quote=F,row.names = F)
	
termNum=5     
if(nrow(kkTab)>=termNum){
	showTerm=row.names(kkTab)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
	pdf(file="IRG_GSEA_KEGG.pdf", width=7.5, height=5.6)
	print(gseaplot)
	dev.off()
}



