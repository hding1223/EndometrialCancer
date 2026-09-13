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

setwd("D:/workspace/TCGA/UCEC/input")
data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)

# Analysis step; see the surrounding code for details.
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data = data[,group==0]

geneFile = 'IRG_UCEC_351.txt'
#geneFile = 'IRG_2499.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# Immune-infiltration analysis step.
sameGene=intersect(row.names(data), as.vector(gene[,1]))
# Immune-infiltration analysis step.
#view(sameGene)
data=data[sameGene,]
data=t(data)
#view(data)

rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))
clinical = read.table('TCIA-ClinicalData-update-MSS.txt',header=T, sep="\t",check.names=F,row.names = 1)

# Analysis step; see the surrounding code for details.
sameSample=intersect(row.names(data),row.names(clinical))
clinical = clinical[sameSample,]
data=data[sameSample,]
WGCNA1=log2(data+1)

# Analysis step; see the surrounding code for details.
gsg <- goodSamplesGenes(WGCNA1, verbose = 3);
gsg[["allOK"]]
if(!gsg[["allOK"]]){
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",paste(names(WGCNA1)[!gsg$goodGenes],collapse=",")))
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing genes:",paste(names(WGCNA1)[!gsg$goodSamples],collapse=",")))
  WGCNA1 = WGCNA1[gsg$goodSamples,gsg$goodGenes]
}

setwd("D:/workspace/TCGA/UCEC/output/MSS/WGCNA-DEGs")
# Analysis step; see the surrounding code for details.
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
# Analysis step; see the surrounding code for details.
table(clust)
# Analysis step; see the surrounding code for details.
keepSamples <- (clust == 1)
WGCNA1 <- WGCNA1[keepSamples, ]
dim(WGCNA1)

# Analysis step; see the surrounding code for details.
# WGCNA module-analysis step.
enableWGCNAThreads()  # WGCNA module-analysis step.
powers = c(1:20)
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sft <- pickSoftThreshold(as.matrix(WGCNA1),
                         powerVector = powers, verbose =5)
pdf(file="2.scale_independence.pdf",width=9,height=5)
par(mfrow=c(1,2))

cex1=0.9
plot(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],  # Analysis step; see the surrounding code for details.
     xlab = "Soft Threshold (power)",  # Analysis step; see the surrounding code for details.
     ylab = "Scale Free Topology Model Fit,signed R^2",  # Analysis step; see the surrounding code for details.
     type = "n",  # Plot the corresponding figure.
     main = paste("Scale independence")) +  # Plot the corresponding figure.
  text(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, 
       cex = cex1,
       col = "steelblue") +
  abline(h = cex1,  # Analysis step; see the surrounding code for details.
         col = "red")

plot(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
     sft$fitIndices[,5],  # Analysis step; see the surrounding code for details.
     xlab = "Soft Threshold (power)",  # Analysis step; see the surrounding code for details.
     ylab = "Mean Connectivity",  # Analysis step; see the surrounding code for details.
     type="n",  # Plot the corresponding figure.
     main = paste("Mean connectivity")) +  # Plot the corresponding figure.
  text(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
       sft$fitIndices[,5],
       labels = powers, 
       cex = cex1, 
       col = "steelblue")

dev.off()
                         
save(sft , file = "wgcna_step2.RData")

# Tumor-microenvironment analysis step.
powerEstimate = sft$powerEstimate
adjacency = adjacency(WGCNA1,power=powerEstimate)
powerEstimate

# Analysis step; see the surrounding code for details.
TOM=TOMsimilarityFromExpr(adjacency)
# Analysis step; see the surrounding code for details.
dissTOM=1-TOM

# Analysis step; see the surrounding code for details.
geneTree = hclust(as.dist(dissTOM),method = "average");
pdf(file="3_gene_clustering.pdf",width=12,height=9)
plot(geneTree,xlab="",sub="",,main="Gene clustering on TOM-based dissmilarity",
     labels = FALSE,hang = 0.04)
dev.off()

minModuleSize=2
dynamicMods <- cutreeDynamic(dendro = geneTree,distM = dissTOM,
                       deepSplit =2,pamRespectsDendro = FALSE,
                       minClusterSize = minModuleSize);  # Analysis step; see the surrounding code for details.

table(dynamicMods)
dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="4_Dynamic_Tree.pdf",width=8,height=6)

# Plot the corresponding figure.
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

# Analysis step; see the surrounding code for details.
nGenes <- ncol(WGCNA1)
nSamples <- nrow(WGCNA1)
sameSample2=intersect(row.names(clinical),rownames(WGCNA1))
datTraits = clinical[sameSample2,]

setwd("D:/workspace/TCGA/UCEC/output/MSS/WGCNA-DEGs")
moduleTraitCor = cor(MEs,datTraits,use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nSamples)
pdf(file="7_Module_trait.pdf",width=6.5,height=5.5)
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
#dim(textMatrix)
#names(datTraits)
par(mar=c(5,10,3,3))
labeledHeatmap(Matrix = moduleTraitCor,  # Plot the corresponding figure.
               xLabels = names(datTraits),  # Analysis step; see the surrounding code for details.
               yLabels = names(MEs),  # Analysis step; see the surrounding code for details.
               ySymbols = names(MEs),  # Analysis step; see the surrounding code for details.
               colorLabels = FALSE,  # Analysis step; see the surrounding code for details.
               colors = blueWhiteRed(50),  # Analysis step; see the surrounding code for details.
               textMatrix = textMatrix,  # Analysis step; see the surrounding code for details.
               setStdMargins = FALSE,  # Analysis step; see the surrounding code for details.
               cex.text = 0.5,  # Analysis step; see the surrounding code for details.
               zlim = c(-1,1),  # Analysis step; see the surrounding code for details.
               main = paste("Module-trait relationships"))  # Plot the corresponding figure.
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
Selectedclinical = "msi_status"
Selectedclinical2 = "msi_status"
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
