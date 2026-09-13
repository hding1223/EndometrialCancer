#install.packages("VennDiagram")
library(VennDiagram)

setwd("D:/workspace/TCGA/UCEC/input/genes")
# WGCNA module-analysis step.
module_grey = read.table('module_grey.txt',header=F, sep="\t",check.names=F)
module_grey = module_grey[,1]
# Differential-expression analysis step.
UCEC_351 = read.table('IRG_UCEC_351.txt',header=F, sep="\t",check.names=F)
UCEC_351 = UCEC_351[,1]

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/Genes")
venn.diagram(x=list('Module Grey'=module_grey,'DEGs'=UCEC_351),filename = 'wgcna_degs_intersect.png',fill=c('grey','darkorange1'))


intersect_grey_Degs  =  intersect(module_grey,UCEC_351)
write.table(intersect_grey_Degs,file="intersect_grey_Degs_292.txt",sep="\t",quote=F,col.name=F)


# WGCNA module-analysis step.
data1 = read.table('IRG_MSS_Gene-OS-uniCox165.txt',header=F, sep="\t",check.names=F)
data1 = data1[,1]
# Differential-expression analysis step.
data2 = read.table('IRG_DEGs_OS_uniCox_45.txt',header=F, sep="\t",check.names=F)
data2 = data2[,1]

data3 = read.table('MSS_DEGs_diff_gene_5.txt',header=F, sep="\t",check.names=F)
data3 = data3[,1]


setwd("D:/workspace/TCGA/UCEC/output/MSS/genes")
venn.diagram(x=list('WGCNA'=data1,'DEGs Group1'=data2,'DEGs Group2'=data3),filename = 'process_before_VN.png',fill=c('dodgerblue','goldenrod1','darkorange1'))

intersect_data12 =  intersect(data1,data2)
intersect_data123 =  intersect(data12,data3)
intersect_data13 =  intersect(data1,data3)

union_data12 =  union(data1,data2)
union_data123 = union(union_data12,data3)

write.table(data12,file="gene_165_intersect_45.txt",sep="\t",quote=F,col.name=F)
write.table(union_data123,file="gene_165_intersect_172.txt",sep="\t",quote=F,col.name=F)
write.table(intersect_data123,file="intersect_data123.txt",sep="\t",quote=F,col.name=F)
write.table(intersect_data13,file="intersect_data13.txt",sep="\t",quote=F,col.name=F)


