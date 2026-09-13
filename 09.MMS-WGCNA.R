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

#geneFile = 'IRG_UCEC_351.txt'
geneFile = 'IRG_2499.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# Process the gene data.
sameGene=intersect(row.names(data), as.vector(gene[,1]))
# Process the gene data.
#view(sameGene)
data=data[sameGene,]
data=t(data)
#view(data)

rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))
clinical = read.table('TCIA-ClinicalData-update-MSS.txt',header=T, sep="\t",check.names=F,row.names = 1)

# Process the sample data.
sameSample=intersect(row.names(data),row.names(clinical))
data=data[sameSample,]
WGCNA1=log2(data+1)

gsg <- goodSamplesGenes(WGCNA1, verbose = 3);
gsg[["allOK"]]

clinical = clinical[sameSample,]
#rt = cbind(clinical,data)
# Process the sample data.
sampleTree <- hclust(dist(WGCNA1), method = "average")

# Process the sample data.
traitColors <- numbers2colors(as.numeric(factor(clinical$msi_status)), 
                              colors = rainbow(length(table(clinical$msi_status))), 
                              signed = FALSE)

# Plot the corresponding figure.
plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(clinical),
                    main = "Sample dendrogram and trait heatmap")

# Plot the corresponding figure.
#clust <- cutreeStatic(sampleTree, cutHeight = , minSize = 100)
# Analysis step; see the surrounding code for details.
clust <- cutreeDynamic(dendro = sampleTree,
                              minClusterSize = 20)  # Analysis step; see the surrounding code for details.

# Analysis step; see the surrounding code for details.
table(clust)


# Process the sample data.
keepSamples <- (clust == 1)
WGCNA1 <- WGCNA1[keepSamples, ]

# Process the sample data.
nGenes <- ncol(WGCNA1)
nSamples <- nrow(WGCNA1)
dim(WGCNA1)

setwd("D:/workspace/TCGA/UCEC/output/MSS/WGCNA")
save(WGCNA1, clinical, nGenes, nSamples, file = "wgcna_step1.RData")

# Process the gene data.
powers <- c(c(1:10), 
            seq(from = 12, 
                to = 30,
                by = 2))

# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
sft <- pickSoftThreshold(as.matrix(WGCNA1),
                         powerVector = powers, 
                         networkType = "signed")

# Save or summarize the results.
powerEstimate = sft$powerEstimate

par(mfrow = c(1,2))  # Plot the corresponding figure.
cex1 = 0.9  # Analysis step; see the surrounding code for details.

# Open a PDF device
pdf("power1.pdf", width = 5, height = 6)  # Adjust width and height as needed

plot(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],  # Analysis step; see the surrounding code for details.
     xlab = "Soft Threshold (power)",  # Analysis step; see the surrounding code for details.
     ylab = "Scale Free Topology Model Fit,signed R^2",  # Analysis step; see the surrounding code for details.
     type = "n",  # Analysis step; see the surrounding code for details.
     main = paste("Scale independence")) +  # Plot the corresponding figure.
  text(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, 
       cex = cex1,
       col = "steelblue") +
  abline(h = 0.80,  # Analysis step; see the surrounding code for details.
         col = "red")

# Save the first plot to the PDF
dev.off()  # Close the current PDF device


# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
pdf("power2.pdf", width = 5, height = 6)

plot(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
     sft$fitIndices[,5],  # Analysis step; see the surrounding code for details.
     xlab = "Soft Threshold (power)",  # Analysis step; see the surrounding code for details.
     ylab = "Mean Connectivity",  # Analysis step; see the surrounding code for details.
     type="n",  # Analysis step; see the surrounding code for details.
     main = paste("Mean connectivity")) +  # Plot the corresponding figure.
  text(sft$fitIndices[,1],  # Analysis step; see the surrounding code for details.
       sft$fitIndices[,5],
       labels = powers, 
       cex = cex1, 
       col = "steelblue")

dev.off()


save(sft , file = "wgcna_step2.RData")


# Analysis step; see the surrounding code for details.
# Set the analysis threshold.
enableWGCNAThreads()  # WGCNA processing step.


if(T){ # Process the gene data.
  # Save or summarize the results.
  net = blockwiseModules(
    as.matrix(WGCNA1),  # Process the gene data.
    # Analysis step; see the surrounding code for details.
    power = 8,
    maxBlockSize = nGenes,  # Analysis step; see the surrounding code for details.
    TOMType = "unsigned",  # Analysis step; see the surrounding code for details.
    minModuleSize = 30,  # Process the gene modules.
    reassignThreshold = 0,  # Set the analysis threshold.
    mergeCutHeight = 0.25,  # Set the analysis threshold.
    numericLabels = TRUE,  # Analysis step; see the surrounding code for details.
    pamRespectsDendro = FALSE,  # Plot the corresponding figure.
    saveTOMs = F,  # Analysis step; see the surrounding code for details.
    verbose = 3  # Analysis step; see the surrounding code for details.
  )
  table(net$colors)  # Process the gene data.
}

# Plot the corresponding figure.

pdf("plotDendroAndColors.pdf", width = 10, height = 6)
if(T){  # Analysis step; see the surrounding code for details.
  # Plot the corresponding figure.
  moduleColors=labels2colors(net$colors)
  table(moduleColors)  # Process the gene data.
  
  # Plot the corresponding figure.
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}
dev.off()

save(moduleColors, net, file="wgcna_step3.RData")

library(forcats)  # Analysis step; see the surrounding code for details.
clinical$msi_status <- factor(clinical$msi_status)  # Analysis step; see the surrounding code for details.
levels(clinical$OS.time)  # Analysis step; see the surrounding code for details.

pdf("Module-trait_relationships.pdf", width = 10, height = 8)

if(T){ 
  nGenes = ncol(WGCNA1)  # Process the gene data.
  nSamples = nrow(WGCNA1)  # Process the sample data.
  design <- model.matrix(~0+clinical$msi_status)  # Analysis step; see the surrounding code for details.
  colnames(design)= levels(clinical$msi_status)  # Analysis step; see the surrounding code for details.
  MES0 <- moduleEigengenes(WGCNA1,moduleColors)$eigengenes  # Process the gene modules.
  MEs = orderMEs(MES0)  # Process the gene modules.
  moduleTraitCor <- cor(MEs,design,use = "p")  # Process the gene modules.
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  # Analysis step; see the surrounding code for details.
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  # Plot the corresponding figure.
  dim(textMatrix)=dim(moduleTraitCor)  # Analysis step; see the surrounding code for details.
  par(mar=c(6, 8.5, 3, 3))  # Plot the corresponding figure.
  labeledHeatmap(Matrix = moduleTraitCor,  # Plot the corresponding figure.
                 xLabels = colnames(design),  # Analysis step; see the surrounding code for details.
                 yLabels = names(MEs),  # Analysis step; see the surrounding code for details.
                 ySymbols = names(MEs),  # Analysis step; see the surrounding code for details.
                 colorLabels = FALSE,  # Analysis step; see the surrounding code for details.
                 colors = blueWhiteRed(50),  # Analysis step; see the surrounding code for details.
                 textMatrix = textMatrix,  # Analysis step; see the surrounding code for details.
                 setStdMargins = FALSE,  # Analysis step; see the surrounding code for details.
                 cex.text = 0.5,  # Analysis step; see the surrounding code for details.
                 zlim = c(-1,1),  # Analysis step; see the surrounding code for details.
                 main = paste("Module-trait relationships"))  # Plot the corresponding figure.
}
dev.off()
save(MEs, file = "wgcna_step4.RData")


pdf("verboseScatterplot.pdf", width = 8, height = 8)

if(T){
  modNames = substring(names(MEs), 3)  # Process the gene modules.
  geneModuleMembership = as.data.frame(cor(WGCNA1, MEs,
                                           use = "p",method = "spearman"))  # Process the gene data.
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))  # Analysis step; see the surrounding code for details.
  names(geneModuleMembership) = paste("MM", modNames, sep="")  # Analysis step; see the surrounding code for details.
  names(MMPvalue) = paste("p.MM", modNames, sep="")  # Analysis step; see the surrounding code for details.
  
  geneTraitSignificance <- as.data.frame(cor(WGCNA1,as.matrix(clinical$msi_status),use = "p"))  # Process the gene data.
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))  # Analysis step; see the surrounding code for details.
  names(geneTraitSignificance)<- paste("GS.",names(clinical$msi_status),sep = "")  # Analysis step; see the surrounding code for details.
  names(GSPvalue)<-paste("GS.",names(clinical$msi_status),sep = "")  # Analysis step; see the surrounding code for details.
  
  selectModule<-c("purple")  # Process the gene modules.
  par(mfrow=c(ceiling(length(selectModule)/2),2))  # Plot the corresponding figure.
  for(module in selectModule){
    column <- match(module,selectModule)  # Process the gene modules.
    print(module)  # Process the gene modules.
    moduleGenes <- moduleColors==module  # Process the gene data.
    
    # Plot the corresponding figure.
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", module, "module"),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
}
dev.off()
# Plot the corresponding figure.


pdf("TOMplot.pdf", width = 15, height = 15)

if(T){
  # WGCNA processing step.
  geneTree = net$dendrograms[[1]]
  # Analysis step; see the surrounding code for details.
  TOM=TOMsimilarityFromExpr(WGCNA1,power=6)
  # Analysis step; see the surrounding code for details.
  dissTOM=1-TOM
  # Analysis step; see the surrounding code for details.
  plotTOM = dissTOM^7
  # Analysis step; see the surrounding code for details.
  diag(plotTOM)=NA
  # Plot the corresponding figure.
  pdf("6_allgene_Network-heatmap.pdf",width = 100,height = 100)
  TOMplot(plotTOM,geneTree,moduleColors,main="Network heapmap plot of all genes")
  
  # Process the gene data.
  nSelect =200
  set.seed(10)
  select=sample(nGenes,size = nSelect)
  # Process the gene data.
  selectTOM = dissTOM[select,select]
  # Process the gene data.
  selectTree = hclust(as.dist(selectTOM),method = "average")
  # Process the gene data.
  selectColors = moduleColors[select]
  # Process the gene data.
  plotDiss=selectTOM^7
  # Analysis step; see the surrounding code for details.
  diag(plotDiss)=NA
  # Plot the corresponding figure.
  TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap of selected gene")
}

dev.off()

# Plot the corresponding figure.
# Plot the corresponding figure.

pdf("plotEigengeneNetworks.pdf", width = 10, height = 10)

if(T){ 
  MEs=moduleEigengenes(WGCNA1,moduleColors)$eigengenes  # Process the gene modules.
  MET = orderMEs(cbind(MEs,clinical$msi_status))  # Process the gene modules.
  par(cex = 0.9)  # Analysis step; see the surrounding code for details.
  # Plot the corresponding figure.
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90,excludeGrey = FALSE)  # Plot the corresponding figure.
}

dev.off()

allModules <- unique(moduleColors)

# Process the gene modules.
for (module in allModules) {
  # Analysis step; see the surrounding code for details.
  if (TRUE) {
    
    # Analysis step; see the surrounding code for details.
    TOM=TOMsimilarityFromExpr(WGCNA1,power=14)
    
    probes = colnames(WGCNA1)  # Process the gene data.
    inModule = (moduleColors == module)  # Process the gene data.
    modProbes = probes[inModule]  # Process the gene data.
    head(modProbes)  # Process the gene data.
    modTOM = TOM[inModule, inModule]  # Process the gene data.
    dimnames(modTOM) = list(modProbes, modProbes)  # Analysis step; see the surrounding code for details.
    # Process the gene data.
    nTop = 100  # Process the gene data.
    IMConn = softConnectivity(WGCNA1[, modProbes])  # Process the gene data.
    top = (rank(-IMConn) <= nTop)  # Process the gene data.
    filterTOM = modTOM[top, top]  # Process the gene data.
    
    # for cytoscape
    cyt = exportNetworkToCytoscape(filterTOM,
                                   edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                   nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0.02,
                                   nodeNames = modProbes[top], 
                                   nodeAttr = moduleColors[inModule][top])  # Analysis step; see the surrounding code for details.
  }
}
