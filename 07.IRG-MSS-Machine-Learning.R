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
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data[[1]],
                        list_train_vali_Data = list_train_vali_Data,
                        unicox.filter.for.candi = T,
                        unicox_p_cutoff = 0.01,
                        candidate_genes = sameGene,
                        mode = 'all',
                        nodesize=5,seed = 123456)

#setwd("D:/workspace/TCGA/UCEC/output/MSS/ML-45_0.001")
#setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma.wilcoxon_20_0.05")
#setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/limma_13_0.01")
setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")

res = readRDS("res.rds")
final_model = "StepCox[forward] + GBM"
mlres =res[["ml.res"]]
stepcoxGBM = mlres[[final_model]]

# Cox regression analysis step.
# Analysis step; see the surrounding code for details.
# Process the gene data.
str(stepcoxGBM, max.level = 1)

# Cox regression analysis step.
genes_in_model <- stepcoxGBM$fit$var.names
cat("入选基因数:", length(genes_in_model), "\n"); print(genes_in_model)

# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
ri <- gbm::relative.influence(stepcoxGBM$fit, n.trees = stepcoxGBM$best)
relInf <- data.frame(gene = names(ri), relative_influence = as.numeric(ri),
                     stringsAsFactors = FALSE)
relInf <- relInf[order(-relInf$relative_influence), ]
row.names(relInf) <- NULL
write.table(relInf, "stepcoxGBM_gene_relInfluence.txt", sep="\t", quote=F, row.names=F)
pdf("stepcoxGBM_gene_relInfluence.pdf", width=7, height=max(4, nrow(relInf)*0.3+2))
barplot(rev(relInf$relative_influence), names.arg=rev(relInf$gene),
        horiz=TRUE, col="#2c7fb8", las=1,
        xlab="Relative influence", main="StepCox[forward]+GBM")
dev.off()

# Analysis step; see the surrounding code for details.
stepcoxGBM$best                                      # Analysis step; see the surrounding code for details.
stepcoxGBM$fit$n.trees; stepcoxGBM$fit$shrinkage     # Analysis step; see the surrounding code for details.
stepcoxGBM$fit$interaction.depth; stepcoxGBM$fit$distribution$name

# Plot the corresponding figure.
# Process the gene data.
pdf("stepcoxGBM_gene_PDP.pdf", width=6, height=5)
for(g in relInf$gene[1:min(10, nrow(relInf))])
  print(plot(stepcoxGBM$fit, i.var=g, n.trees=stepcoxGBM$best,
             main=paste("PDP:", g, "vs. cumulative hazard"),
             xlab=paste(g, "expression (log2(TPM+1))"),
             ylab="Partial dependence (predicted cumulative hazard)"))
dev.off()

# Process the gene data.
tt <- train; colnames(tt) <- gsub("-",".",colnames(tt))
rs  <- res[["riskscore"]][[final_model]][["Train"]]        # Analysis step; see the surrounding code for details.
m   <- merge(rs, tt, by = "ID")
dirCor <- sapply(genes_in_model, function(g) cor(m[[g]], m$RS))
relInf$cor_with_RS <- round(dirCor[match(relInf$gene, names(dirCor))], 3)
write.table(relInf, "stepcoxGBM_gene_params.txt", sep="\t", quote=F, row.names=F)

# Cox regression analysis step.
sc <- mlres[["StepCox[forward]"]]                    # Analysis step; see the surrounding code for details.
stepCoxTab <- data.frame(
  gene   = names(coef(sc)),
  coef   = round(coef(sc), 4),
  HR     = round(exp(coef(sc)), 3),
  pvalue = signif(summary(sc)$coefficients[, "Pr(>|z|)"], 3))
write.table(stepCoxTab, "stepCox_forward_gene_coef.txt", sep="\t", quote=F, row.names=F)

# Analysis step; see the surrounding code for details.
#saveRDS(res,"res.rds")

# Analysis step; see the surrounding code for details.
write.table(res[["Cindex.res"]],file="Cindex.res.txt",sep="\t",quote=F,col.name=T,row.names = F)
write.table(res[["Sig.genes"]],file="unicox.txt",sep="\t",quote=F,col.name=T,row.names = F)


# Analysis step; see the surrounding code for details.
#res = readRDS("res.rds")
#mlres =res[["ml.res"]]
#stepcox = mlres[["StepCox[forward] + Ridge"]]
#riskscore = res[["riskscore"]]
#riskscore_stepcox = riskscore[["StepCox[forward] + Ridge"]]

#C-index
pdf(file="1.cindex_dis_all.pdf",width=8,height=15)
cindex_dis_all(res,validate_set=names(list_train_vali_Data)[-1],order=names(list_train_vali_Data),width=0.35)
dev.off()

#final_model = "survival - SVM"
#final_model = "StepCox[forward] + Ridge"
#final_model = "StepCox[both] + survival−SVM"
#final_model = "StepCox[forward] + Enet[α=0.1]"
#final_model = "StepCox[forward] + GBM"
final_model = "StepCox[both] + Ridge"
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

# Process the gene data.
#res[["ml.res"]]

pdf(file="cindex_dis_final_model.pdf")
cindex_dis_select(res,
                  model=final_model,
                  order= c("Train","Validation"))

#pdf(file="Survival_StepCox[forward] + Ridge.pdf")
#pdf(file="Survival_survival_SVM.pdf")
#pdf(file="StepCox[forward]_Enet.pdf")
#pdf(file="Survival_StepCox[forward] + GBM.pdf")
pdf(file="StepCox[both] + Ridge.pdf")
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = final_model,dataset = names_list[i],
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
dev.off()

# Analysis step; see the surrounding code for details.
for(i in c(1:length(list_train_vali_Data))){
  pdf(file=paste0("3.survival.",names(list_train_vali_Data)[i],".pdf"),width=8,height=8)
  rs_sur(res,model_name = final_model,dataset = names(list_train_vali_Data)[i],
         median.line ="hv",
         cutoff = 0.7,
         conf.int = T,
         xlab = "Day",pval.coord = c(1000,0.9)
         )
  dev.off()
}

# Calculate or plot AUC values.
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode='all',AUC_time = 1,
                             auc_cal_method = "KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode='all',AUC_time = 3,
                             auc_cal_method = "KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode='all',AUC_time = 5,
                             auc_cal_method = "KM")

# Calculate or plot AUC values.
for(i in c(1,3,5)){
  pdf(file=paste0("4.auc_dis_all.",i,"y.pdf"),width=8,height=15)
  auc_dis_all(eval(parse(text = paste0("all.auc.",i,"y"))),
              dataset = names(list_train_vali_Data),
              validate_set = names(list_train_vali_Data)[-1],
              order = names(list_train_vali_Data),
              width = 0.35,
              year = i
              )
  dev.off()
}

#pdf(file="auc_135_survival_SVM.pdf",width=8,height=4)
#pdf(file="auc_135_StepCox[forward]_Enet.pdf",width=8,height=4)
#pdf(file="auc_135_StepCox[forward]_Ridge.pdf",width=8,height=4)
#pdf(file="auc_135_StepCox[forward]_GBM.pdf",width=8,height=4)
pdf(file="auc_135_StepCox[both]_Ridge.pdf",width=8,height=4)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name=final_model,
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
dev.off()

# Plot ROC curves.
for(i in c(1,3,5)){
  pdf(file=paste0("5.roc_vis.",i,"y.pdf"),width=8,height=8)
  roc_vis(eval(parse(text = paste0("all.auc.",i,"y"))),
              model_name = final_model,
              dataset = names(list_train_vali_Data),
              order = names(list_train_vali_Data),
              anno_position = c(0.65,0.55),
              year = i
  )
  dev.off()
}