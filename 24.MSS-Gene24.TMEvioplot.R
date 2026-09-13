#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")

library(limma)
library(reshape2)
library(ggpubr)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)
rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")

setwd("D:/workspace/TCGA/UCEC/input")
scoreFile="IRG_TMEscores.txt"
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
dimnames=list(rownames(score),colnames(score))
score=matrix(as.numeric(as.matrix(score)), nrow=nrow(score), dimnames = dimnames)


score=score[,1:3]
rownames(score)=substr(rownames(score),1,12)
rownames(score)=gsub('[.]', '-', rownames(score))

sameSample=intersect(row.names(rs), row.names(score))
rs=rs[sameSample,"Type",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(score, rs)
#view(rt)
rt$Type=factor(rt$Type, levels=c("Low", "High"))

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "scoreType", "Score")

p=ggviolin(data, x="scoreType", y="Score", fill = "Type",
           xlab="",
           ylab="TME score",
           legend.title="Group",
           add = "boxplot", add.params = list(color="white"),
           palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")


setwd("D:/workspace/TCGA/UCEC/output/MSS-2/TME")
pdf(file="IRG_Immune_vioplot.pdf", width=6, height=5)
print(p1)
dev.off()
