#install.packages("VennDiagram")
library(VennDiagram)

setwd("D:/workspace/TCGA/UCEC/input/Drugs")

#Score_Drugs
data1 = read.table('Score_Drugs.txt',header=F, sep="\t",check.names=F)
data1 = data1[,1]
#MSS_Drugs
data2 = read.table('MSS_Drugs.txt',header=F, sep="\t",check.names=F)
data2 = data2[,1]

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/compare_drugs")

intersect_Drugs  =  intersect(data1,data2)
write.table(intersect_Drugs,file="intersect_Drugs.txt",sep="\t",quote=F,col.name=F)
venn.diagram(x=list('Group2'=data1,'Group1'=data2),filename = 'compare_drugs.png',fill=c('dodgerblue','darkorange1'),width = 2400,  # Analysis step; see the surrounding code for details.
             height = 1600)

