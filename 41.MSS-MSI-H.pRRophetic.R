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

#pFilter=0.001            #pvalue
pFilter=0.05            #pvalue

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
#rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")

# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
#msi_status <- read.table('TCIA-ClinicalData_MMRd.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- read.table('TCIA-ClinicalData_MSS.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- as.matrix(msi_status)  # Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(rs), rownames(msi_status))
rs <- rs[sameSample, ]
msi_status <- msi_status[sameSample, ]
rs=cbind(msi_status, rs)
#rs = data

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/MSS_Drugs")
#drug = 'Sunitinib'
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
	rs=rs[sameSample, "msi_status", drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(rs, senstivity)
	rt<- data.frame(rt)
	#rt$Type=factor(rt$Type, levels=c("Low", "High"))
	rt$msi_status=factor(rt$msi_status, levels=c("MSS","MSI-H"))
	type=levels(factor(rt[,"msi_status"]))
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	rt$senstivity = as.numeric(rt$senstivity)
	test=wilcox.test(senstivity~msi_status, data=rt)
	diffPvalue=test$p.value

	print(diffPvalue)
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt, x="msi_status", y="senstivity", fill="msi_status",
					      xlab='Sub group',
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title='Sub group',
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}