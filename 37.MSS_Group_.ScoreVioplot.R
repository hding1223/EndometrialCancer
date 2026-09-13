# Load required packages.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")

library(limma)
library(reshape2)
library(ggpubr)

# Read input data.
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma_13_0.01")      
scoreFile <- "final_mode_riskscore.txt"
score <- read.table(scoreFile, header=FALSE, sep="\t", check.names=FALSE, row.names=1)
data <- as.matrix(score)  # Analysis step; see the surrounding code for details.

# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
msi_status <- read.table('ClinicalData_MSS_Group.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
#msi_status <- read.table('TCIA-ClinicalData_MSS.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- as.matrix(msi_status)  # Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(data), rownames(msi_status))
score_filtered <- data[sameSample, ]
group_filtered <- msi_status[sameSample, ]
#view(score_filtered)
#view(group_filtered)
# Analysis step; see the surrounding code for details.
df <- cbind(score_filtered, group_filtered)
colnames(df) <- c("Score", "group")  # Analysis step; see the surrounding code for details.

df <- as.data.frame(df)

# Analysis step; see the surrounding code for details.
table(df$group)
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/MSS_GROUP")
p <- ggboxplot(df, x = "group", y = "Score", fill = "group", 
               palette = "jco",
               xlab = "", ylab = "Risk Score",
               , add = "jitter",  short.panel.labs = FALSE) 
# Analysis step; see the surrounding code for details.
p <-   p + stat_compare_means(method = "wilcox.test",label = "p.format")
pdf(file=paste0("RS_MSS_MSI-H.pdf"),width=5.5,height=5)
print(p)
dev.off()

