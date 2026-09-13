#install.packages(c("tidyverse", "magrittr", "WGCNA"))
#install.packages("WGCNA")   # WGCNA is available on CRAN
#BiocManager::install("impute")
#BiocManager::install("ggrepel")
#BiocManager::install("cowplot")
#BiocManager::install("ggthemes")
library(WGCNA)
library("survival")
library("survminer")
library(tidyverse)
library(caret)

setwd("D:/workspace/TCGA/UCEC/input")
data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)

#01-09: tumor; 11-19: normal; 20-29: adjacent normal tissue
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data = data[,group==0]

#geneFile = 'IRG_UCEC_351.txt'
geneFile = 'IRG_2499.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# Identify the intersection between genes in the expression data and immune-related genes
sameGene=intersect(row.names(data), as.vector(gene[,1]))
# Extract the target genes (geneName) and their expression values from the raw expression matrix, then transpose the matrix
#view(sameGene)
data=data[sameGene,]
data=t(data)
#view(data)

rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))
clinical = read.table('TCGA_UCEC_clinical_selected_above_32_MMS.txt',header=T, sep="\t",check.names=F,row.names = 1)

# Identify shared samples
sameSample=intersect(row.names(data),row.names(clinical))
clinical = clinical[sameSample,]
data=data[sameSample,]

data = log2(data+1)

# Randomly split samples stratified by OS; only the expression matrix is used for WGCNA
set.seed(123)
index <- createDataPartition(clinical$OS, p = 0.7, list = FALSE)
trainSamples <- rownames(clinical)[index]
write.table(trainSamples, file = "D:/workspace/TCGA/UCEC/output/append/01.WGCNA/trainSamples.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

train <- data[trainSamples, ]
test  <- data[setdiff(rownames(data), trainSamples), ]

WGCNA1=train

view(WGCNA1)
# Check for missing values and outliers
gsg <- goodSamplesGenes(WGCNA1, verbose = 3);
gsg[["allOK"]]
if(!gsg[["allOK"]]){
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",paste(names(WGCNA1)[!gsg$goodGenes],collapse=",")))
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing genes:",paste(names(WGCNA1)[!gsg$goodSamples],collapse=",")))
  WGCNA1 = WGCNA1[gsg$goodSamples,gsg$goodGenes]
}

setwd("D:/workspace/TCGA/UCEC/output/append/01.WGCNA")
# Perform hierarchical clustering of samples with hclust; dist calculates distances between gene-expression profiles
sampleTree <- hclust(dist(WGCNA1), method = "average")
pdf(file="1_sample_cluster.1.pdf",width=12,height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main = "Sample clustering to detect Outliner", sub="",xlab="",
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
dev.off()

cutHeight = 70
pdf(file="1_sample_cluster.2.pdf",width=12,height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main = "Sample clustering to detect Outliner", sub="",xlab="",
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(h=cutHeight,col="red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
# Check the size of each cluster
table(clust)
# Remove groups or samples that do not meet the criteria based on cutreeStatic results, retaining only eligible samples
keepSamples <- (clust == 1)
WGCNA1 <- WGCNA1[keepSamples, ]
dim(WGCNA1)

# 3. Construct the co-expression network
# Construct the co-expression network using the optimal soft-thresholding power, assign genes to modules, and plot the gene-clustering tree
enableWGCNAThreads()  # Enable multithreading to accelerate WGCNA
powers = c(1:20)
# Use pickSoftThreshold to calculate the scale-free topology model fit and mean connectivity for each power value,
# and identify the optimal power
sft <- pickSoftThreshold(as.matrix(WGCNA1),
                         powerVector = powers, verbose =5)
pdf(file="2.scale_independence.pdf",width=9,height=5)
par(mfrow=c(1,2))

cex1=0.9
plot(sft$fitIndices[,1],  # x-axis: power value
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],  # y-axis: scale-free topology fit
     xlab = "Soft Threshold (power)",  # x-axis label
     ylab = "Scale Free Topology Model Fit,signed R^2",  # y-axis label
     type = "n",  # Do not draw points or lines
     main = paste("Scale independence")) +  # Plot title
  text(sft$fitIndices[,1],  # Add labels to each point
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, 
       cex = cex1,
       col = "steelblue") +
  abline(h = cex1,  # Add a horizontal line
         col = "red")

plot(sft$fitIndices[,1],  # x-axis: power value
     sft$fitIndices[,5],  # y-axis: mean connectivity
     xlab = "Soft Threshold (power)",  # x-axis label
     ylab = "Mean Connectivity",  # y-axis label
     type="n",  # Do not draw points or lines
     main = paste("Mean connectivity")) +  # Plot title
   text(sft$fitIndices[,1],  # Add labels to each point
       sft$fitIndices[,5],
       labels = powers, 
       cex = cex1, 
       col = "steelblue")

dev.off()

save(sft , file = "wgcna_step2.RData")

# Identify the optimal power from the sft results and assign it to powerEstimate
powerEstimate = sft$powerEstimate
adjacency = adjacency(WGCNA1,power=powerEstimate)
powerEstimate

# Calculate the TOM similarity matrix
TOM=TOMsimilarityFromExpr(adjacency)
# Convert the TOM similarity matrix into a dissimilarity matrix
dissTOM=1-TOM

# Gene clustering
geneTree = hclust(as.dist(dissTOM),method = "average");
pdf(file="3_gene_clustering.pdf",width=12,height=9)
plot(geneTree,xlab="",sub="",,main="Gene clustering on TOM-based dissmilarity",
     labels = FALSE,hang = 0.04)
dev.off()

minModuleSize=100
dynamicMods <- cutreeDynamic(dendro = geneTree,distM = dissTOM,
                             deepSplit =2,pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize);  # Adjust minClusterSize as needed

table(dynamicMods)
dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="4_Dynamic_Tree.pdf",width=8,height=6)

# Plot the sample dendrogram and phenotype heatmap; groupLabels adds phenotype labels
plotDendroAndColors(geneTree, 
                    dynamicColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE,hang=0.03,
                    addGuide = TRUE,guideHang = 0.05,
                    main = "Gene dendrogram and Module colors")

dev.off()

MEList = moduleEigengenes(WGCNA1,colors=dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss),method="average")
pdf(file="5_Clustering_module.pdf",,width=8,height=6)
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
MEDissThres=0.3
abline(h=MEDissThres,col="red")
dev.off()

# Calculate the gene co-expression network for the retained samples and obtain the numbers of genes and samples
nGenes <- ncol(WGCNA1)
nSamples <- nrow(WGCNA1)
sameSample2=intersect(row.names(clinical),rownames(WGCNA1))
datTraits = clinical[sameSample2,]

setwd("D:/workspace/TCGA/UCEC/output/append/01.WGCNA")
moduleTraitCor = cor(MEs,datTraits,use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSamples)
pdf(file="7_Module_trait.pdf",width=6.5,height=5.5)
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
#dim(textMatrix)
#names(datTraits)
par(mar=c(5,10,3,3))
labeledHeatmap(Matrix = moduleTraitCor,  # Plot a labeled heatmap
               xLabels = names(datTraits),  # x-axis labels
               yLabels = names(MEs),  # y-axis labels
               ySymbols = names(MEs),  # y-axis symbols
               colorLabels = FALSE,  # Do not display color labels
               colors = blueWhiteRed(50),  # Color range
               textMatrix = textMatrix,  # Display the text matrix
               setStdMargins = FALSE,  # Do not use standard margins
               cex.text = 0.5,  # Text size
               zlim = c(-1,1),  # Color-mapping range
               main = paste("Module-trait relationships"))  # Plot title
dev.off()

moduleColors =dynamicColors
probes = colnames(WGCNA1)
geneInfo0 = data.frame(probes = probes,moduleColor=moduleColors)
geneOrder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder,]
write.table(geneInfo,file = "module_all.txt",sep="\t",row.names = F,quote=F)

for(mod in 1:nrow(table(moduleColors))){
  modules = names(table(moduleColors))[mod]
  probes = colnames(WGCNA1)
  inModule = (moduleColors==modules)
  modGenes = probes[inModule]
  write.table(modGenes,file=paste0("module_",modules,".txt"),sep="\t",row.names = F,quote=F)
}

module = "grey"
Selectedclinical = "MSI_STATUS"
Selectedclinical2 = "MSI_STATUS"
Selectedclinical = as.data.frame(datTraits[,Selectedclinical]);
names(Selectedclinical) = "Selectedclinical";
modNames = substring(names(MEs),3)
WGCNA2 = WGCNA1[rownames(MEs),]
geneModuleMembership = as.data.frame(cor(WGCNA2,MEs,use="p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM",modNames,sep="");
names(MMPvalue) = paste("p.MM",modNames,sep="")
geneTraitSignificance = as.data.frame(cor(WGCNA2,Selectedclinical,use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples));
names(geneTraitSignificance) = paste("GS.",names(Selectedclinical),sep="");
names(GSPvalue) = paste("p.GS",names(Selectedclinical),sep="");

column = match(module,modNames)
moduleGenes = moduleColors==module
outPdf=paste(Selectedclinical2,"_",module,".pdf",sep="")
pdf(file=outPdf,width=7,height=7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab=paste("Module Membership in", module,"module"),
                   ylab=paste("Gene significance for ",Selectedclinical2,sep=""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab=1.2,cex.axis=1.2,col=module)

dev.off()

datMM=cbind(geneModuleMembership[,paste("MM",module,sep="")],geneTraitSignificance)
colnames(datMM)[1]=paste("MM",module,sep="")
geneSigFilter=0.1
moduleSigFilter=0.4
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
datMM=datMM[abs(datMM[,1])>moduleSigFilter,]
write.table(row.names(datMM),file=paste0("hubGenes",module,".txt"),sep="\t",row.names = F,quote=F)
