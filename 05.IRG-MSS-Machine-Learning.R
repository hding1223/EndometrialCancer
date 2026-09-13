#if(!requireNamespace("BiocManager", quietly = TRUE)) BiocManager::install("BiocManager")
#depens <- c('GSEBase','GSVA','cancerclass','mixOmics','sparrow','sva','ComplexHeatmap')
#for(i in 1:length(depens)){
#   depen <- depens[i]
#   if(!requireNamespace(depen,quietly = TRUE)) BiocManager::install(depen,update = FALSE)
#}

#if(!requireNamespace("CoxBoost",quietly = TRUE))
#   devtools::install_github("binderh/CoxBoost")

#if(!requireNamespace("fastAdaboost",quietly = TRUE))
#  devtools::install_github("souravc83/fastAdaboost")

#if(!requireNamespace("Mime",quietly = TRUE))
#  devtools::install_github("l-magnificence/Mime")

#install.packages("MatrixModels")
#install.packages("Hmisc")
#install.packages("polspline")
#install.packages("multcomp")
#install.packages("survival")
#install.packages("Design")
#install.packages("MASS")
#install.packages("survcomp")
#BiocManager::install("limma")

#install.packages("caret")
#devtools::install_version("Matrix", version = "1.6.0")
#update.packages("Matrix", version = "1.6.0")
#BiocManager::install("rms", force = TRUE)
#install.packages("xgboost")

library("xgboost")
library("survival")
library(caret)
# Analysis step; see the surrounding code for details.
library(Mime1)
library(limma)
library(tidyverse)

library(rms)

names_list <- c("Train","Validation")

setwd("D:/workspace/TCGA/UCEC/input")
# Analysis step; see the surrounding code for details.
clinical = read.table('clinical_above30_ml.txt',header=T, sep="\t",check.names=F,row.names = 1)
#clinical = read.table('clinical_above30_ml_MMS.txt',header=T, sep="\t",check.names=F,row.names = 1)
#clinical$OS.time = clinical$OS.time/365

data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)
data = data[rowMeans(data)>1,]

setwd("D:/workspace/TCGA/UCEC/input/genes")
#geneFile = 'IRG_OS_45.txt'
#geneFile = 'gene_intersect_42.txt'
#geneFile = 'gene_union_172.txt'
#geneFile = 'gene_unicox_63.txt'
#geneFile = 'IRG.diff.limma.MSS.txt'
#geneFile = 'IRG.diff.limma.wilcoxon.MSS.txt'
geneFile = 'IRG.diff.limma_ml3.MSS.txt'
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# Process the gene data.
sameGene=intersect(row.names(data), as.vector(gene[,1]))

data=data[sameGene,]
# Analysis step; see the surrounding code for details.
data=t(data)
rownames(data)=substr(rownames(data),1,12)
rownames(data)=gsub('[.]', '-', rownames(data))

# Process the sample data.
sameSample=intersect(row.names(data),row.names(clinical))
data=data[sameSample,]
data = log2(data+1)
clinical = clinical[sameSample,]
rt = cbind(rownames(data),clinical,data)
colnames(rt)[1] = "ID"
#view(rt)
# Analysis step; see the surrounding code for details.
set.seed(123)
index <- createDataPartition(rt$OS, p = 0.7, list = FALSE)
train <- rt[index, ]
test <- rt[-index, ]
#view(train)


list_train_vali_Data = list(train,test)
names(list_train_vali_Data) = c("Train","Validation")

#res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data[[1]],
#                                              candidate_genes = sameGene,
#                                              mode = "all",nodesize =5,seed = 5201314 )


#view(list_train_vali_Data[[1]])
#res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data[[1]],
#                        list_train_vali_Data = list_train_vali_Data,
#                        unicox.filter.for.candi = T,
#                        unicox_p_cutoff = 0.01,
#                        candidate_genes = sameGene,
#                        mode = 'all',
#                        nodesize=5,seed = 123456)

#setwd("D:/workspace/TCGA/UCEC/output/MSS/ML-45_0.001")
#setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma.wilcoxon_20_0.05")
#setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma_13_0.01")
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")

res = readRDS("res.rds")
final_model = "StepCox[forward] + GBM"
mlres =res[["ml.res"]]
stepcoxGBM = mlres[[final_model]]


# Analysis step; see the surrounding code for details.
setwd("D:/workspace/TCGA/UCEC/output/append/05.cindex")
write.table(res[["Cindex.res"]],file="Cindex.res.txt",sep="\t",quote=F,col.name=T,row.names = F)
write.table(res[["Sig.genes"]],file="unicox.txt",sep="\t",quote=F,col.name=T,row.names = F)


#C-index
pdf(file="1.cindex_dis_all.pdf",width=8,height=15)
cindex_dis_all(res,validate_set=names(list_train_vali_Data)[-1],order=names(list_train_vali_Data),width=0.35)
dev.off()

#final_model = "survival - SVM"
#final_model = "StepCox[forward] + Ridge"
#final_model = "StepCox[both] + survival−SVM"
#final_model = "StepCox[forward] + Enet[α=0.1]"
#final_model = "StepCox[forward] + GBM"
#final_model = "StepCox[both] + Ridge"
# Analysis step; see the surrounding code for details.
write.table(res[["riskscore"]][[final_model]][["Train"]],file="riskscore.Train.txt",sep="\t",quote=F,col.name=T,row.names = F)
write.table(res[["riskscore"]][[final_model]][["Validation"]],file="riskscore.Validation.txt",sep="\t",quote=F,col.name=T,row.names = F)


RS = data.frame()
for(i in seq_along(res[["riskscore"]])){
  for(l in c(1:length(list_train_vali_Data))){
    rtd = cbind(res[["riskscore"]][[i]][[l]][,c(1,4)],rep(names(res[["riskscore"]][i]),length(rownames(res[["riskscore"]][[i]][[l]][,c(1,4)]))))
    colnames(rtd)[3]="model"
    RS= rbind(RS,rtd)
  }
}
write.table(RS,file="riskscore.txt",sep="\t",quote=F,col.name=T,row.names = F)

# Calculate or evaluate the concordance index.
# Bootstrap resampling and confidence-interval calculation.
# Calculate or evaluate the concordance index.
# Analysis step; see the surrounding code for details.
# Bootstrap resampling and confidence-interval calculation.

valData <- res[["riskscore"]][[final_model]][["Validation"]]
# Analysis step; see the surrounding code for details.
valDF <- data.frame(
  time   = as.numeric(valData[[grep("time|Time", colnames(valData), value=TRUE)[1]]]),
  status = as.numeric(valData[[grep("^OS$|status|event", colnames(valData), value=TRUE)[1]]]),
  rs     = as.numeric(valData[[ncol(valData)]]),   # Analysis step; see the surrounding code for details.
  ID     = as.character(valData$ID),
  stringsAsFactors = FALSE)

# Plot the corresponding figure.
# Plot the corresponding figure.
# Analysis step; see the surrounding code for details.
# Calculate or evaluate the concordance index.
valDF$rs <- -valDF$rs

nBootCI  <- 1000
seedBoot <- 20240   # Analysis step; see the surrounding code for details.
set.seed(seedBoot)
bootCidx <- numeric(nBootCI)
for(b in seq_len(nBootCI)){
  idxB   <- sample(seq_len(nrow(valDF)), size=nrow(valDF), replace=TRUE)
  dBoot  <- valDF[idxB, ]
  # Bootstrap resampling and confidence-interval calculation.
  if(length(unique(dBoot$status)) < 2) { bootCidx[b] <- NA; next }
  sc     <- survival::concordance(Surv(time, status) ~ rs, data=dBoot)
  bootCidx[b] <- as.numeric(sc$concordance)
}

ciLo <- quantile(bootCidx, 0.025, na.rm=TRUE)
ciHi <- quantile(bootCidx, 0.975, na.rm=TRUE)
valCidxFull <- as.numeric(survival::concordance(Surv(time, status) ~ rs, data=valDF)$concordance)

bootSummary <- data.frame(
  model        = final_model,
  dataset      = "Validation",
  nVal         = nrow(valDF),
  cindex_point = round(valCidxFull, 4),
  ci95_lower   = round(ciLo, 4),
  ci95_upper   = round(ciHi, 4),
  nBoot        = nBootCI,
  seed         = seedBoot)
write.table(bootSummary, file="Validation_Cindex_bootstrap95CI.txt",
            sep="\t", quote=F, row.names=F)

pdf(file="Validation_Cindex_bootstrap.pdf", width=6, height=5)
hist(bootCidx, breaks=8, col="#2c7fb8", border="white",
     main=paste0("Bootstrap C-index (", final_model, ", Validation)"),
     xlab="C-index", ylab="Frequency")
abline(v=valCidxFull, col="black", lty=2, lwd=2)
abline(v=c(ciLo, ciHi), col="red", lty=3, lwd=2)
legend("topleft", legend=c(sprintf("Point est. = %.3f", valCidxFull),
                           sprintf("95%% CI = [%.3f, %.3f]", ciLo, ciHi)),
       col=c("black","red"), lty=c(2,3), lwd=2, bty="n", cex=0.9)
dev.off()

cat(sprintf("\n[Validation C-index] %.4f (bootstrap 95%% CI: %.4f - %.4f, nBoot=%d)\n",
            valCidxFull, ciLo, ciHi, nBootCI))

# =============================================================================
# Analysis step; see the surrounding code for details.
# =============================================================================
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Plot the corresponding figure.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Set the analysis threshold.
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Calculate or evaluate the concordance index.
# Cox regression analysis step.
# -----------------------------------------------------------------------------
# Analysis step; see the surrounding code for details.
#   Input  : normalized expression matrix (log2 TPM+1) of GAL, PGR, THRB
#   Load   : trained StepCox[forward] + GBM object (fit, best)
#   Predict: RS = predict(fit, newdata, n.trees = best, type = "link")
#   Output : risk score per patient
#   Cutoff : predefined Training-set median (for external validation)
# =============================================================================

# Analysis step; see the surrounding code for details.
hp <- data.frame(
  parameter = c("distribution", "n.trees.grown", "cv.folds", "shrinkage",
                "interaction.depth", "n.minobsinnode", "best.trees",
                "selected.genes", "unicox.p.cutoff", "stepcox.direction",
                "n.train", "n.events.train", "n.validation"),
  value = c("coxph", "10000", "10",
            format(stepcoxGBM$fit$shrinkage), format(stepcoxGBM$fit$interaction.depth),
            format(stepcoxGBM$fit$n.minobsinnode), format(stepcoxGBM$best),
            paste(stepcoxGBM$fit$var.names, collapse = "/"), "0.01", "forward",
            nrow(train), sum(train$OS == 1), nrow(test)))
write.table(hp, "model_hyperparameters.txt", sep="\t", quote=F, row.names=F)

saveRDS(list(model_name = final_model,
             fit        = stepcoxGBM$fit,     # Analysis step; see the surrounding code for details.
             best       = stepcoxGBM$best,    # Analysis step; see the surrounding code for details.
             var.names  = stepcoxGBM$fit$var.names,
             hyperparams= hp),
        "StepCoxForward_GBM_model_object.rds")

# Analysis step; see the surrounding code for details.
genes3   <- stepcoxGBM$fit$var.names          # "GAL" "PGR" "THRB"
rsStore  <- list(Train = res[["riskscore"]][[final_model]][["Train"]],
                 Validation = res[["riskscore"]][[final_model]][["Validation"]])
chkTab <- data.frame()
for(ds in c("Train", "Validation")){
  dd  <- if(ds == "Train") train else test                       # Analysis step; see the surrounding code for details.
  rsR <- as.numeric(predict(stepcoxGBM$fit, newdata = dd[, genes3],
                            n.trees = stepcoxGBM$best, type = "link"))
  chk <- merge(rsStore[[ds]][, c("ID","RS")],
               data.frame(ID = dd$ID, RS_rebuilt = rsR), by = "ID")
  chkTab <- rbind(chkTab, data.frame(
    dataset = ds, n = nrow(chk),
    max_abs_diff = max(abs(chk$RS - chk$RS_rebuilt)),
    pearson_cor  = cor(chk$RS, chk$RS_rebuilt),
    reproducible = all(abs(chk$RS - chk$RS_rebuilt) < 1e-6)))
}
write.table(chkTab, "RiskScore_reproducibility_check.txt",
            sep="\t", quote=F, row.names=F)
print(chkTab)   # Analysis step; see the surrounding code for details.

# Set the analysis threshold.
# Plot the corresponding figure.
# Analysis step; see the surrounding code for details.
rsStore$Train$RS      <- -rsStore$Train$RS
rsStore$Validation$RS <- -rsStore$Validation$RS
cutTrainMed <- median(rsStore$Train$RS)
grpTab <- data.frame()
for(ds in c("Train","Validation")){
  rs <- rsStore[[ds]]
  grpTab <- rbind(grpTab, data.frame(
    dataset = ds, cutoff = "training-set median", cutoff_value = round(cutTrainMed, 4),
    nHigh = sum(rs$RS >  cutTrainMed), nLow = sum(rs$RS <= cutTrainMed)))
}
write.table(grpTab, "riskgroup_training_median.txt", sep="\t", quote=F, row.names=F)
# Plot the corresponding figure.
# Analysis step; see the surrounding code for details.
