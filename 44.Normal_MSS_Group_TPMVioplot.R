library(VennDiagram)
library(pheatmap)
library(stringr)
library(ggplot2)
library(tidyverse)
library(limma)
library(reshape2)
library(ggpubr)
library(tidyverse)

#gene = "PGR"
#gene = "GAL"
gene = "THRB"

setwd("D:/workspace/TCGA/UCEC")
data = read.table('TCGA_UCEC_TPM.txt',header=T,sep="\t",check.names = F,row.names = 1)
dimnames = list(rownames(data),colnames(data))
data = matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
#view(data)
data = data[rowMeans(data)>1,]
#view(data)
# Analysis step; see the surrounding code for details.
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)

data1 = data[,group==1]
#view(data1)
# Analysis step; see the surrounding code for details.
data1=t(data1)
rownames(data1)=substr(rownames(data1),1,12)
rownames(data1)=gsub('[.]', '-', rownames(data1))
rt_normal = log2(data1[,gene]+1)
view(rt_normal)
# Analysis step; see the surrounding code for details.
new_vector <- rep('Normal', times = 35)
df_normal <- data.frame("TPM" = rt_normal, group = new_vector)
colnames(df_normal) <- c("TPM", "group")  # Analysis step; see the surrounding code for details.


data=t(data)
rownames(data)=substr(rownames(data),1,12)


# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
#msi_status <- read.table('ClinicalData_MSS_Group.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- read.table('TCIA-ClinicalData_MSS_WITHOUTNA.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- as.matrix(msi_status)  # Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(data), rownames(msi_status))
tpm_filtered <- data[sameSample, ]
group_filtered <- msi_status[sameSample, ]

rt = log2(tpm_filtered[,gene]+1)
#view(rt)
df <- data.frame("TPM" = rt, group = group_filtered)
colnames(df) <- c("TPM", "group")  # Analysis step; see the surrounding code for details.
df <- as.data.frame(df)

df <- rbind(df,df_normal)
# Analysis step; see the surrounding code for details.
table(df$group)
cmpgroup=levels(factor(df$group))
df$group = factor(df$group,levels=cmpgroup)
#view(df)
comp = combn(cmpgroup,2)
#view(comp)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]] <- comp[,i]}
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/genes")
p=ggboxplot(df, x = "group", y = "TPM", fill = "group", 
            palette = "jco",
            xlab = "", ylab = paste0(gene," log2(TPM+1)"),
            , add = "jitter",  short.panel.labs = FALSE) 
#  Add p-value
p1 = p + stat_compare_means(comparisons = my_comparisons,label="p.format") #default Wilcoxon
#p2 = p + stat_compare_means(comparisons = my_comparisons,label="p.format",method = "t.test")
pdf(file=paste0(gene,"_Normal_MSS_MSI-H.pdf"),width=5.5,height=5)
print(p1)
dev.off()

