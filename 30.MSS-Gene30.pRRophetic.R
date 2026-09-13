#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("limma", "car", "ridge", "preprocessCore", "genefilter", "sva"))
#BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("pRRophetic")
#install.packages("E:/R/package/pRRophetic_0.5.tar.gz", repos = NULL, type = "source")
#???ð?
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
library(tidyverse)

#trace(calcPhenotype, edit = T)
#trace(summarizeGenesByMean, edit = T)
#trace(pmap, edit = T)

set.seed(12345)

pFilter=0.001            #pvalue

expFile='TCGA_UCEC_TPM.txt'   
setwd("D:/workspace/TCGA/UCEC/input") 


data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))


setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)
rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/drugs")
#view(allDrugs)
for(drug in allDrugs){
	possibleError=tryCatch(
    	{senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    	error=function(e) {cat("Error:", e$message, "\n")})
	print("before error")
  if(inherits(possibleError, "error")){next}
	senstivity=senstivity[senstivity!="NaN"]
	senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
	print("after error")
	sameSample=intersect(row.names(rs), names(senstivity))
	rs=rs[sameSample, "Type", drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(rs, senstivity)


	rt$Type=factor(rt$Type, levels=c("Low", "High"))
	type=levels(factor(rt[,"Type"]))
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	

	test=wilcox.test(senstivity~Type, data=rt)
	diffPvalue=test$p.value

	print(diffPvalue)
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt, x="Type", y="senstivity", fill="Type",
					      xlab='RiskScore',
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title='RiskScore',
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}