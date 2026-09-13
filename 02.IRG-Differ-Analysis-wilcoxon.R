#install.packages("VennDiagram")
#install.packages("pheatmap")
#install.packages("stringr")
#install.packages("ggplot2")

library(VennDiagram)
library(pheatmap)
library(stringr)
library(ggplot2)
library(tidyverse)
library(ggrepel)

setwd("D:/workspace/TCGA/UCEC/input")
data = read.table('TCGA_UCEC_TPM.txt',header=T,sep="\t",check.names = F,row.names = 1)
dimnames = list(rownames(data),colnames(data))
data = matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
#view(data)
data = data[rowMeans(data)>1,]

# Process the gene data.
geneFile = 'IRG_2499.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=data[sameGene,]

#view(data)
# Analysis step; see the surrounding code for details.
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)

group0_count =length(group[group==0])
group1_count =length(group[group==1])
group2_count =length(group[group==2])

group=gsub("2","1",group)

#view(data)
data1 = data[,group==1]
data2 = data[,group==0]
# WGCNA processing step.
trainSamples=read.table('trainSamples.txt', header=F, sep="\t", check.names=F, stringsAsFactors=F)[,1]
tumorPts=gsub('[.]','-',substr(colnames(data2),1,12))
data2=data2[,tumorPts %in% trainSamples,drop=FALSE]
data=cbind(data1,data2)
conNum=ncol(data1)
treatNum=ncol(data2)
Type=c(rep(1,conNum),rep(2,treatNum))
# Process the sample data.
colnames(data)=substr(colnames(data),1,12)
colnames(data)=gsub('[.]', '-', colnames(data))
view(data)

outTab=data.frame()
cut_off_FDR =0.05 # Set or apply the FDR threshold.
cut_off_log2FC =1 # Set or apply the log2 fold-change threshold.
for(i in row.names(data)){
  rt=data.frame(expression=data[i,],Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  fdr=p.adjust(pvalue,method="fdr")
  sig = ifelse(fdr < cut_off_FDR &    # Set the analysis threshold.
                 abs(logFC) >= cut_off_log2FC,  # Analysis step; see the surrounding code for details.
               ifelse(logFC > cut_off_log2FC ,'Up','Down'),'no')
  if((logFC>0 & diffMed >0)|(logFC <0 & diffMed < 0)){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue,sig=sig))
  }
}



#view(outTab)
pValue = outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
log10FDR = -log10(fdr)
outTab=cbind(outTab,log10FDR=log10FDR)

logFCfilter=1
fdrFilter=0.05
outDiff=outTab[(abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
#view(outDiff)
setwd("D:/workspace/TCGA/UCEC/output/append/0.2.Immune_Related_Genes")
table(outDiff$sig)

write.csv(outDiff,'IRG.diff.wilcoxon.csv')
outputdata <- read.csv(file = "IRG.diff.wilcoxon.csv")

# Plot the corresponding figure.
ggplot(data = outputdata, mapping = aes(x =logFC, y=log10FDR)) + # Analysis step; see the surrounding code for details.
  geom_point(aes(color=sig),alpha=0.65, size=2) +  # Analysis step; see the surrounding code for details.
  # Analysis step; see the surrounding code for details.
  geom_vline(xintercept = c(-1,1), linetype = "dashed")+
  # Analysis step; see the surrounding code for details.
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  scale_color_manual(values=c("green","#ff4757", "#d2dae2"))+
  labs(x="log2FC", y="-log10FDR") + # Analysis step; see the surrounding code for details.
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()) +
  # Analysis step; see the surrounding code for details.
  geom_text_repel(data = outputdata, 
                  aes(label = gene,color=sig),
                  size=2.5)


# Plot the corresponding figure.
ggsave(filename="IRG.diff.wilcoxon-Volcano.png")


#geneNum=50
#outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
#diffGeneName=as.vector(outDiff[,1])
#diffLength=length(diffGeneName)
#hmGene=c()
#if(diffLength>(2*geneNum)){
#   hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
#}else{
#   hmGene=diffGeneName
#}
#hmExp=log2(data[hmGene,]+0.01)
#Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
#names(Type)=colnames(data)
#Type=as.data.frame(Type)
#pdf(file="IRG_heatmap.pdf",width=10,height=6.5)
#pheatmap(hmExp,annotaion=Type,color=colorRampPalette(c(rep('blue',5),'white',rep("red",5)))(50),
#         cluster_cols = F,
#         show_rownames = F,
#         scale="row",
#         fontsize = 8,
#         fontsize_row = 5,
#         fontsize_col = 8)
#dev.off()

#pdf(file="IRG_Volcano.pdf",width=5,height=5)
#xMax=6
#yMax=max(-log10(outTab$fdr))+1
#plot(as.numeric(as.vector(outTab$logFC)),-log10(outTab$fdr),xlab = "logFC",ylab = "-log10(fdr)",
#     main="Volcano",ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20,cex=1.2)
#diffSub=subset(outTab,fdr<fdrFilter&as.numeric(as.vector(logFC))>logFCfilter)
#points(as.numeric(as.vector(diffSub$logFC)),-log10(diffSub$fdr),pch=20,col="red",cex=1.5)
#diffSub=subset(outTab,fdr<fdrFilter&as.numeric(as.vector(logFC))<(-logFCfilter))
#points(as.numeric(as.vector(diffSub$logFC)),-log10(diffSub$fdr),pch=20,col="green",cex=1.5)
#abline(v=0,lty=2,lwd=3)
#dev.off()




