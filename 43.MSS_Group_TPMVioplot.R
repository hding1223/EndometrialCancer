# Load required packages.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")

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
# Analysis step; see the surrounding code for details.
#data=log2(data+1)
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
#view(df)

colnames(df) <- c("TPM", "group")  # Analysis step; see the surrounding code for details.

df <- as.data.frame(df)

# Analysis step; see the surrounding code for details.
table(df$group)
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/genes")
p <- ggboxplot(df, x = "group", y = "TPM", fill = "group", 
               palette = "jco",
               xlab = "", ylab = paste0(gene," log2(TPM+1)"),
               , add = "jitter",  short.panel.labs = FALSE) 
# Analysis step; see the surrounding code for details.
p <-   p + stat_compare_means(method = "wilcox.test",label = "p.format")
pdf(file=paste0(gene,"_MSS_MSI-H.pdf"),width=5.5,height=5)
print(p)
dev.off()

