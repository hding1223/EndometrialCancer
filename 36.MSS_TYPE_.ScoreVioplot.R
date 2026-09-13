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
#setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma_13_0.01")
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
scoreFile <- "final_mode_riskscore.txt"
score <- read.table(scoreFile, header=FALSE, sep="\t", check.names=FALSE, row.names=1)
dimnames = list(rownames(score),colnames(score))
data = matrix(as.numeric(as.matrix(score)),nrow=nrow(score),dimnames=dimnames)

# Read input data.
setwd("D:/workspace/TCGA/UCEC/input/clinical") 
#msi_status <- read.table('ClinicalData_MSS_Group.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- read.table('TCIA-ClinicalData_MSS.txt', header=TRUE, sep="\t", check.names=FALSE, row.names=1)
msi_status <- as.matrix(msi_status)  # Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sameSample <- intersect(rownames(data), rownames(msi_status))
score_filtered <- data[sameSample, ]
group_filtered <- msi_status[sameSample, ]

typeof(score_filtered)
typeof(group_filtered)
#view(score_filtered)
#view(group_filtered)
# Analysis step; see the surrounding code for details.


df <- data.frame("Score" = score_filtered, group = group_filtered)
typeof(df)

colnames(df) <- c("Score", "group")  # Analysis step; see the surrounding code for details.


# Analysis step; see the surrounding code for details.
table(df$group)

cmpgroup=levels(factor(df$group))
df$group = factor(df$group,levels=cmpgroup)
#view(df)
comp = combn(cmpgroup,2)
#view(comp)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]] <- comp[,i]}

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/MSS_GROUP")
#boxplot=ggboxplot(df,x="group",y="Score",fill="group",xlab="",ylab=paste("Risk Score"),
#                  legend.title="") + stat_compare_means(comparisons = my_comparisons,label="p.signif")

boxplot=ggboxplot(df, x = "group", y = "Score", fill = "group", 
          palette = "jco",
          xlab = "", ylab = "Risk Score",
          , add = "jitter",  short.panel.labs = FALSE)  + stat_compare_means(comparisons = my_comparisons,label="p.format")

pdf(file=paste0("RS_MSS_MSI-H_63.pdf"),width=5.5,height=5)
print(boxplot)
dev.off()


