#install.packages("VennDiagram")
library(VennDiagram)

setwd("D:/workspace/TCGA/UCEC/input/genes")
# WGCNA module-analysis step.
data1 = read.table('Ding1.txt',header=F, sep="\t",check.names=F)
data1 = data1[,1]
# Differential-expression analysis step.
data2 = read.table('Xiao.txt',header=F, sep="\t",check.names=F)
data2 = data2[,1]

data3 = read.table('Sang.txt',header=F, sep="\t",check.names=F)
data3 = data3[,1]

data4 = read.table('Zheng.txt',header=F, sep="\t",check.names=F)
data4 = data4[,1]

data5 = read.table('Ding2.txt',header=F, sep="\t",check.names=F)
data5 = data5[,1]

setwd("D:/workspace/TCGA/UCEC/output/MSS/genes")
venn.diagram(x=list('Ding\'s Model 1'=data1,'Xiao\'s Model'=data2,'Sang\'s Model'=data3,'Zheng\'s Model'=data4),filename = 'process_after1_VN.png',fill=c('dodgerblue','goldenrod1','darkorange1','green'))

venn.diagram(x=list('Ding\'s Model 2'=data5,'Xiao\'s Model'=data2,'Sang\'s Model'=data3,'Zheng\'s Model'=data4),filename = 'process_after2_VN.png',fill=c('dodgerblue','goldenrod1','darkorange1','green'))
