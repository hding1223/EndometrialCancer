#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)


setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)


# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
msi_status <- read.table('TCIA-ClinicalData_MMRd.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
#msi_status <- read.table('TCIA-ClinicalData_MSS.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- as.matrix(msi_status)  # Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(rs), rownames(msi_status))
rs <- rs[sameSample, ]
msi_status <- msi_status[sameSample, ]
data=cbind(msi_status, rs)
rs = data
#rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")

tciaFile="TCIA.txt"        
setwd("D:/workspace/TCGA/UCEC/input") 

ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(ips), row.names(rs))
ips=ips[sameSample, , drop=F]
#rs=rs[sameSample, "Type", drop=F]
rs=rs[sameSample, , drop=F]
data=cbind(ips, rs)

#group=levels(factor(data$Type))
#group=levels(factor(data$group))
data$msi_status=factor(data$msi_status, levels=c("MSS", "MSI-H_epi", "MSI-H_mut"))
group=levels(factor(data$msi_status))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/MSS_IPS")
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "msi_status")]
	gg1=ggviolin(rt, x="msi_status", y=i, fill = "msi_status", 
	         xlab="", ylab=i,
	         legend.title="Sub group",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, "_MSI-H.pdf"), width=4.8, height=4.25)
	print(gg1)
	dev.off()
}
