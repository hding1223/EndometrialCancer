#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")

library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(tidyverse)


setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)


tmbFile="TMB.txt"        
setwd("D:/workspace/TCGA/UCEC/input")
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
ucec=tmb[tmb$CancerType=="UCEC",1,drop=F]
rownames(ucec)=substr(rownames(ucec),1,12)
rownames(ucec)=gsub('[.]', '-', rownames(ucec))

sameSample=intersect(row.names(rs), row.names(ucec))
rs=rs[sameSample,,drop=F]
ucec=ucec[sameSample,,drop=F]
rt=cbind(rs, ucec)
view(rt)

x=as.numeric(rt[,"riskscore"])
y=log2(as.numeric(rt[,"TMB"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0("riskscore", " expression"))+ylab("Tumor mutation burden")+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/TME")
pdf(file="IRG-TMB-cor.pdf",width=5,height=5)
print(p2)
dev.off()

