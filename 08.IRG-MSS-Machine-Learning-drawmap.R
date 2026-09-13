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
#install.packages("forestploter")
library("survival")
library(caret)
# Analysis step; see the surrounding code for details.
library(Mime1)
library(limma)
library(tidyverse)

library(rms)

setwd("D:/workspace/TCGA/UCEC/output/Immune_Related_Genes")

names_list <- c("Train","Validation")

# Analysis step; see the surrounding code for details.
res = readRDS("res.rds")
final_model = "StepCox[forward] + plsRcox"

pdf(file="cindex_dis_final_model.pdf")
cindex_dis_select(res,
                  model=final_model,
                  order= c("Train","Validation"))
dev.off()

pdf(file="Survival_StepCox_plsRcox.pdf")
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + plsRcox",dataset = names_list[i],
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
dev.off()

#pdf(file="metamodel_StepCox_plsRcox.pdf")
#unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "StepCox[forward] + plsRcox",type ='categorical')
#metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
#meta_unicox_vis(metamodel,dataset=names_list)
#dev.off()

pdf(file="auc_135_StepCox_plsRcox.pdf",width=8,height=4)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + plsRcox",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
dev.off()

