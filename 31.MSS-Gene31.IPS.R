#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)


setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)
rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")

tciaFile="TCIA.txt"        
setwd("D:/workspace/TCGA/UCEC/input") 

ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(ips), row.names(rs))
ips=ips[sameSample, , drop=F]
rs=rs[sameSample, "Type", drop=F]
data=cbind(ips, rs)

group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c("Low", "High"))
group=levels(factor(data$Type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/TME")
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "Type")]
	gg1=ggviolin(rt, x="Type", y=i, fill = "Type", 
	         xlab="", ylab=i,
	         legend.title="RiskScore",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=4.8, height=4.25)
	print(gg1)
	dev.off()
}
