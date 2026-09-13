# =============================================================================
# 04.stability_pipeline.R
# WGCNA processing step.
# WGCNA processing step.
# Cox regression analysis step.
#
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# WGCNA parameter configuration.
# Process the gene data.
# Save or summarize the results.
# Bootstrap resampling and confidence-interval calculation.
#
# -----------------------------------------------------------------------------
# Cox regression analysis step.
# Analysis step; see the surrounding code for details.
# -----------------------------------------------------------------------------
#
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Calculate or evaluate the concordance index.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# WGCNA processing step.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Bootstrap resampling and confidence-interval calculation.
# Bootstrap resampling and confidence-interval calculation.
#
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
#    - RSF：rfsrc(ntree=1000, nodesize=5, splitrule="logrank")；
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Cox regression analysis step.
# Analysis step; see the surrounding code for details.
# WGCNA parameter configuration.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Set or apply the FDR threshold.
# Cox regression analysis step.
#
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Save or summarize the results.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# WGCNA processing step.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
#
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Cox regression analysis step.
# Calculate or evaluate the concordance index.
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Process the gene data.
# Process the gene data.
# Process the gene data.
# Analysis step; see the surrounding code for details.
#
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Set or apply the FDR threshold.
# Analysis step; see the surrounding code for details.
# Set or apply the FDR threshold.
# Analysis step; see the surrounding code for details.
# Set or apply the FDR threshold.
# Analysis step; see the surrounding code for details.
# Set or apply the FDR threshold.
# Analysis step; see the surrounding code for details.
# Set the analysis threshold.
# Analysis step; see the surrounding code for details.
# Process the sample data.
#
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Analysis step; see the surrounding code for details.
# Analysis step; see the surrounding code for details.
# Process the gene data.
# Process the gene data.
# WGCNA processing step.
# WGCNA processing step.
# Analysis step; see the surrounding code for details.
# Process the sample data.
#
# Analysis step; see the surrounding code for details.
# =============================================================================

library(WGCNA)
library(survival)
library(caret)
library(tidyverse)
library(pheatmap)

# Analysis step; see the surrounding code for details.
nSplit        <- 50        # Analysis step; see the surrounding code for details.
pTrain        <- 0.7       # Analysis step; see the surrounding code for details.
nBoot         <- 0         # Bootstrap resampling and confidence-interval calculation.
seedBaseSplit <- 1L        # Analysis step; see the surrounding code for details.
seedBaseBoot  <- 10000L    # Bootstrap resampling and confidence-interval calculation.

# WGCNA parameter configuration.
netType       <- "unsigned"  # Analysis step; see the surrounding code for details.
minModSize    <- 100         # Analysis step; see the surrounding code for details.
deepSplitVal  <- 2
mergeCutH     <- 0           # Analysis step; see the surrounding code for details.
powers        <- 1:20

# Set the analysis threshold.
degFCthr      <- 1         # Set or apply the log2 fold-change threshold.
degFDRthr     <- 0.05      # Set or apply the FDR threshold.

# Cox regression analysis step.
coxPthr       <- 0.05      # Cox regression analysis step.
freqThr       <- 0.8       # Set the analysis threshold.

# Analysis step; see the surrounding code for details.
useParallel   <- TRUE
nCores        <- min(4L, parallel::detectCores() - 1L)

# Analysis step; see the surrounding code for details.
inputDir  <- "D:/workspace/TCGA/UCEC/input"
outRoot   <- "D:/workspace/TCGA/UCEC/output/append/stability"

# Analysis step; see the surrounding code for details.
setwd(inputDir)

# Process the sample data.
tpm <- read.table('TCGA_UCEC_TPM.txt', header=T, sep="\t", check.names=F, row.names=1)
dimnames <- list(rownames(tpm), colnames(tpm))
tpm <- matrix(as.numeric(as.matrix(tpm)), nrow=nrow(tpm), dimnames=dimnames)

# Process the sample data.
grp <- sapply(strsplit(colnames(tpm), "\\-"), "[", 4)
grp <- sapply(strsplit(grp, ""), "[", 1)
isTumor  <- grp == "0"
isNormal <- grp %in% c("1", "2")   # Analysis step; see the surrounding code for details.

# Process the sample data.
irg <- read.table('IRG_2499.txt', header=F, sep="\t", check.names=F, stringsAsFactors=F)[,1]
bgGenes <- intersect(rownames(tpm), irg)
rawAll <- tpm[bgGenes, , drop=FALSE]
rawAll <- rawAll[rowMeans(rawAll) > 1, , drop=FALSE]
bgGenes <- rownames(rawAll)

# Process the sample data.
patientOf <- function(cn) gsub('[.]', '-', substr(cn, 1, 12))
collapseToPatient <- function(mat){
  pts <- patientOf(colnames(mat))
  up  <- unique(pts)
  out <- t(sapply(up, function(p) rowMeans(mat[, pts == p, drop=FALSE])))
  rownames(out) <- up
  t(out)   # Process the gene data.
}
rawTumor  <- collapseToPatient(rawAll[, isTumor,  drop=FALSE])  # Process the gene data.
rawNormal <- collapseToPatient(rawAll[, isNormal, drop=FALSE])  # Process the gene data.
logTumor  <- log2(rawTumor + 1)                                 # WGCNA processing step.
rawNormal <- rawNormal[bgGenes, , drop=FALSE]                   # Process the gene data.

# Analysis step; see the surrounding code for details.
cliAll <- read.table('TCGA_UCEC_clinical_selected_above_32_MMS.txt',
                     header=T, sep="\t", check.names=F, row.names=1)
cliAll <- cliAll[, c("OS", "OS.time")]

# Analysis step; see the surrounding code for details.
cohortPts <- intersect(colnames(logTumor), rownames(cliAll))
cohortCli <- cliAll[cohortPts, ]
logTumor  <- logTumor[, cohortPts, drop=FALSE]
rawTumor  <- rawTumor[, cohortPts, drop=FALSE]

cat(sprintf("[预处理] 基因背景 %d | 正常对照 %d | 队列患者 %d | OS事件 %d\n",
            length(bgGenes), ncol(rawNormal), length(cohortPts), sum(cohortCli$OS == 1)))

# Analysis step; see the surrounding code for details.
run_one <- function(task){
  # task: list(runId, mode, seed)
  library(WGCNA); library(survival); library(caret)
  disableWGCNAThreads()   # Analysis step; see the surrounding code for details.
  set.seed(task$seed, kind="Mersenne-Twister", normal.kind="Inversion")

  # Analysis step; see the surrounding code for details.
  if(task$mode == "split"){
    idx <- createDataPartition(cohortCli$OS, p=pTrain, list=FALSE)
    trainPts <- rownames(cohortCli)[idx]
  } else {
    trainPts <- sample(rownames(cohortCli), size=nrow(cohortCli), replace=TRUE)
  }

  # WGCNA processing step.
  exprT <- t(logTumor[, trainPts, drop=FALSE])   # Process the sample data.
  gsg <- goodSamplesGenes(exprT, verbose=0)
  exprT <- exprT[gsg$goodSamples, gsg$goodGenes, drop=FALSE]

  # WGCNA processing step.
  # Analysis step; see the surrounding code for details.
  sft <- pickSoftThreshold(as.matrix(exprT), powerVector=powers, verbose=0)
  power <- sft$powerEstimate
  if(is.na(power)){   # Analysis step; see the surrounding code for details.
    fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    hit <- which(fit >= 0.8)[1]
    power <- ifelse(is.na(hit), 6, sft$fitIndices[hit, 1])
  }
  adj  <- adjacency(exprT, power=power)
  tom  <- TOMsimilarityFromExpr(adj)
  dissTOM <- 1 - tom
  geneTree <- hclust(as.dist(dissTOM), method="average")
  dynMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=deepSplitVal,
                           pamRespectsDendro=FALSE, minClusterSize=minModSize, verbose=0)
  dynColors <- labels2colors(dynMods)
  if(mergeCutH > 0){
    merged <- mergeCloseModules(exprT, dynColors, cutHeight=mergeCutH, verbose=0)
    moduleColors <- merged$colors
  } else {
    moduleColors <- dynColors   # Analysis step; see the surrounding code for details.
  }
  greyGenes <- colnames(exprT)[moduleColors == "grey"]

  # Set or apply the FDR threshold.
  tumRaw <- rawTumor[, trainPts, drop=FALSE]
  conMean <- rowMeans(rawNormal); treatMean <- rowMeans(tumRaw)
  logFC  <- log2(treatMean) - log2(conMean)
  diffMed <- apply(tumRaw, 1, median) - apply(rawNormal, 1, median)
  degP <- setNames(rep(NA_real_, length(bgGenes)), bgGenes)
  for(g in seq_along(bgGenes)){
    degP[g] <- tryCatch(
      wilcox.test(as.numeric(rawNormal[g, ]), as.numeric(tumRaw[g, ]))$p.value,
      error = function(e) NA_real_)
  }
  degFDR <- p.adjust(degP, method="BH")
  degSig <- !is.na(degFDR) & degFDR < degFDRthr & abs(logFC) >= degFCthr &
            ((logFC > 0 & diffMed > 0) | (logFC < 0 & diffMed < 0))
  degGenes <- bgGenes[degSig]

  # Process the gene modules.
  interGenes <- intersect(greyGenes, degGenes)

  # Cox regression analysis step.
  exprCox <- t(logTumor[, trainPts, drop=FALSE])   # Process the sample data.
  cliCox  <- cohortCli[trainPts, ]                 # Bootstrap resampling and confidence-interval calculation.
  nEvents <- sum(cliCox$OS == 1)

  emptyRes <- data.frame(gene=character(0), HR=numeric(0), HR.L=numeric(0),
                         HR.H=numeric(0), p=numeric(0), zph.p=numeric(0),
                         fdr=numeric(0), selected=logical(0), stringsAsFactors=FALSE)
  geneRes <- emptyRes
  if(length(interGenes) > 0 && nEvents >= 5){
    for(g in interGenes){
      df <- data.frame(OS=cliCox$OS, OS.time=cliCox$OS.time, x=as.numeric(exprCox[, g]))
      fit <- tryCatch(coxph(Surv(OS.time, OS) ~ x, data=df), error=function(e) NULL)
      if(is.null(fit)) next
      cs <- summary(fit)
      zphP <- tryCatch(cox.zph(fit)$table["x", "p"], error=function(e) NA_real_)
      geneRes <- rbind(geneRes, data.frame(
        gene=g,
        HR=cs$conf.int[1, "exp(coef)"],
        HR.L=cs$conf.int[1, "lower .95"],
        HR.H=cs$conf.int[1, "upper .95"],
        p=cs$coefficients[1, "Pr(>|z|)"],
        zph.p=zphP, stringsAsFactors=FALSE))
    }
    if(nrow(geneRes) > 0){
      geneRes$fdr <- p.adjust(geneRes$p, method="BH")
      geneRes$selected <- geneRes$p < coxPthr
    }
  }
  # Save or summarize the results.
  if(nrow(geneRes) > 0){
    geneRes <- cbind(runId=task$runId, mode=task$mode, seed=task$seed, geneRes,
                     stringsAsFactors=FALSE)
  } else {
    geneRes <- cbind(data.frame(runId=character(0), mode=character(0), seed=integer(0),
                                stringsAsFactors=FALSE), emptyRes)
  }

  runInfo <- data.frame(runId=task$runId, mode=task$mode, seed=task$seed,
                        nSamples=length(trainPts), nEvents=nEvents,
                        power=power, nModules=length(unique(moduleColors)),
                        nGrey=length(greyGenes), nDEG=length(degGenes),
                        nInter=length(interGenes),
                        nCox=sum(geneRes$selected), stringsAsFactors=FALSE)

  list(geneRes=geneRes, runInfo=runInfo,
       sets=list(grey=greyGenes, deg=degGenes, inter=interGenes,
                 cox=geneRes$gene[geneRes$selected]))
}

# Analysis step; see the surrounding code for details.
tasks <- c(lapply(seq_len(nSplit), function(k)
  list(runId=sprintf("split_%03d", k), mode="split", seed=seedBaseSplit + k - 1L)),
  if(nBoot > 0) lapply(seq_len(nBoot), function(k)
    list(runId=sprintf("boot_%03d", k), mode="boot", seed=seedBaseBoot + k - 1L)))

dir.create(outRoot, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outRoot, "runs"), recursive=TRUE, showWarnings=FALSE)

cat(sprintf("[执行] 共 %d 轮任务（%d split + %d boot），并行核数 %d\n",
            length(tasks), nSplit, nBoot, if(useParallel) nCores else 1))

res <- if(useParallel){
  cl <- parallel::makeCluster(nCores)
  on.exit(parallel::stopCluster(cl), add=TRUE)
  # Analysis step; see the surrounding code for details.
  parallel::clusterExport(cl, c("run_one","cohortCli","logTumor","rawTumor","rawNormal","bgGenes",
                                "pTrain","netType","powers","minModSize","deepSplitVal",
                                "mergeCutH","degFCthr","degFDRthr","coxPthr"),
                          envir=environment())
  parallel::parLapply(cl, tasks, function(t) tryCatch(run_one(t),
                          error=function(e) list(error=conditionMessage(e), task=t)))
} else {
  lapply(tasks, function(t) tryCatch(run_one(t),
                          error=function(e) list(error=conditionMessage(e), task=t)))
}

# Analysis step; see the surrounding code for details.
errRuns <- Filter(function(r) !is.null(r$error), res)
if(length(errRuns) > 0){
  cat(sprintf("[警告] %d 轮失败，首个错误信息：%s\n",
              length(errRuns), errRuns[[1]]$error))
  cat("失败轮次：", paste(sapply(errRuns, function(r) r$task$runId), collapse=","), "\n")
}
res <- Filter(function(r) is.null(r$error), res)
if(length(res) == 0){
  stop("所有轮次均失败：请根据上方首个错误信息排查（并行模式可先设 useParallel=FALSE 调试单轮）。")
}

# Save or summarize the results.
long    <- do.call(rbind, lapply(res, `[[`, "geneRes"))
runInfo <- do.call(rbind, lapply(res, `[[`, "runInfo"))
setsAll <- setNames(lapply(res, `[[`, "sets"), sapply(res, function(r) r$runInfo$runId))
if(is.null(long) || nrow(long) == 0){
  print(runInfo)
  stop("所有轮次均未产生交集基因的Cox结果：请查看上方 run_summary（可能事件数过少或 grey∩DEG 为空）。")
}

write.table(long,    file.path(outRoot, "summary_long.tsv"),  sep="\t", quote=F, row.names=F)
write.table(runInfo, file.path(outRoot, "run_summary.tsv"),   sep="\t", quote=F, row.names=F)
saveRDS(setsAll, file.path(outRoot, "stage_sets.rds"))
writeLines(capture.output(sessionInfo()), file.path(outRoot, "sessionInfo.txt"))

# Analysis step; see the surrounding code for details.
spLong <- long[long$mode == "split", ]
spRuns <- runInfo[runInfo$mode == "split", ]

freqDf <- spLong %>%
  group_by(gene) %>%
  summarise(freq   = mean(selected),
            nRun   = n(),
            medHR  = median(HR, na.rm=TRUE),
            hrLo   = quantile(HR, 0.025, na.rm=TRUE),
            hrHi   = quantile(HR, 0.975, na.rm=TRUE),
            dirCons= mean(sign(HR) == sign(median(HR, na.rm=TRUE)), na.rm=TRUE),
            phViol = mean(zph.p < 0.05, na.rm=TRUE),
            medP   = median(p, na.rm=TRUE),
            .groups="drop") %>%
  arrange(desc(freq), medP)
write.table(freqDf, file.path(outRoot, "gene_frequency.tsv"), sep="\t", quote=F, row.names=F)

robust <- freqDf %>% filter(freq >= freqThr)
write.table(robust, file.path(outRoot, "robust_genes.tsv"), sep="\t", quote=F, row.names=F)

# Bootstrap resampling and confidence-interval calculation.
if(nBoot > 0){
  bootCi <- long[long$mode == "boot", ] %>%
    group_by(gene) %>%
    summarise(bootHRmed=median(HR, na.rm=TRUE),
              bootLo=quantile(HR, 0.025, na.rm=TRUE),
              bootHi=quantile(HR, 0.975, na.rm=TRUE), .groups="drop")
  write.table(bootCi, file.path(outRoot, "bootstrap_HR_CI.tsv"), sep="\t", quote=F, row.names=F)
}

# Analysis step; see the surrounding code for details.
jaccard <- function(a, b){ u <- union(a, b); if(length(u) == 0) 1 else length(intersect(a, b)) / length(u) }
pairJ <- function(sets){
  if(length(sets) < 2) return(numeric(0))   # Analysis step; see the surrounding code for details.
  idx <- combn(seq_along(sets), 2)
  sapply(seq_len(ncol(idx)), function(k) jaccard(sets[[idx[1, k]]], sets[[idx[2, k]]]))
}
splitSets <- setsAll[spRuns$runId]
jGrey <- pairJ(lapply(splitSets, `[[`, "grey"))
jDeg  <- pairJ(lapply(splitSets, `[[`, "deg"))
jCox  <- pairJ(lapply(splitSets, `[[`, "cox"))
jDf <- rbind(data.frame(stage="WGCNA grey module", j=jGrey),
             data.frame(stage="Differential genes (DEG)", j=jDeg),
             data.frame(stage="Cox-selected genes", j=jCox))
jMean <- jDf %>% group_by(stage) %>% summarise(m=mean(j), .groups="drop")

cat("\n========== 稳定性结果 ==========\n")
cat(sprintf("稳健基因（选择频率≥%.1f）：%d 个\n", freqThr, nrow(robust)))
print(jMean)
cat(sprintf("轮均事件数 %.1f | 轮均grey %.0f | 轮均DEG %.0f | 轮均交集 %.0f | 轮均Cox入选 %.1f\n",
            mean(spRuns$nEvents), mean(spRuns$nGrey), mean(spRuns$nDEG),
            mean(spRuns$nInter), mean(spRuns$nCox)))

# Plot the corresponding figure.
theme_pub <- theme_minimal(base_size=12) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        panel.grid.minor=element_blank())

# Plot the corresponding figure.
topN <- min(30, nrow(freqDf))
topGenes <- freqDf$gene[seq_len(topN)]
p1 <- freqDf %>%
  filter(gene %in% topGenes) %>%
  mutate(gene=factor(gene, levels=rev(topGenes)),
         dir=ifelse(medHR > 1, "Risk (HR>1)", "Protective (HR<1)")) %>%
  ggplot(aes(x=gene, y=freq, fill=dir)) +
  geom_col(width=0.7, color="grey30", linewidth=0.2) +
  geom_hline(yintercept=freqThr, linetype="dashed", color="red", linewidth=0.6) +
  coord_flip() + scale_fill_manual(values=c("Risk (HR>1)"="#d7301f", "Protective (HR<1)"="#2c7fb8")) +
  labs(x=NULL, y="Selection frequency (fraction of 50 runs with p<0.05)", fill=NULL,
       title=sprintf("Cox gene selection frequency Top%d (red line=%.1f)", topN, freqThr)) +
  theme_pub + theme(legend.position="bottom")
ggsave(file.path(outRoot, "1_selection_frequency.pdf"), p1, width=7, height=8)
ggsave(file.path(outRoot, "1_selection_frequency.png"), p1, width=7, height=8, dpi=300)

# Plot the corresponding figure.
if(nrow(jDf) > 0){
  p2 <- ggplot(jDf, aes(x=j, fill=stage)) +
    geom_histogram(bins=30, alpha=0.8, color="white", linewidth=0.2) +
    geom_vline(data=jMean, aes(xintercept=m), linetype="dashed", color="black", linewidth=0.6) +
    geom_text(data=jMean, aes(x=m, y=Inf, label=sprintf("mean %.2f", m)),
              vjust=1.5, hjust=-0.05, size=3.5) +
    facet_wrap(~stage, ncol=1, scales="free_y") +
    scale_fill_manual(values=c("#4daf4a", "#ff7f00", "#984ea3")) +
    labs(x="Pairwise Jaccard similarity across runs", y="Number of run pairs", title="Cross-run stability of gene sets at each stage") +
    theme_pub + theme(legend.position="none")
  ggsave(file.path(outRoot, "2_jaccard_stability.pdf"), p2, width=7, height=8)
  ggsave(file.path(outRoot, "2_jaccard_stability.png"), p2, width=7, height=8, dpi=300)
}

# Plot the corresponding figure.
if(nrow(robust) > 0){
  p3 <- robust %>%
    mutate(gene=factor(gene, levels=rev(gene)),
           dir=ifelse(medHR > 1, "Risk (HR>1)", "Protective (HR<1)")) %>%
    ggplot(aes(x=medHR, y=gene, color=dir)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40") +
    geom_errorbarh(aes(xmin=hrLo, xmax=hrHi), height=0.2, color="grey50", linewidth=0.5) +
    geom_point(size=3) +
    scale_color_manual(values=c("Risk (HR>1)"="#d7301f", "Protective (HR<1)"="#2c7fb8")) +
    labs(x="Hazard ratio (cross-run median, whiskers=2.5%-97.5% empirical CI)", y=NULL,
         color=NULL, title=sprintf("Robust genes (frequency>=%.1f): %d genes", freqThr, nrow(robust))) +
    theme_pub + theme(legend.position="bottom")
  ggsave(file.path(outRoot, "3_robust_forest.pdf"), p3, width=8,
         height=max(4, nrow(robust)*0.35 + 2))
  ggsave(file.path(outRoot, "3_robust_forest.png"), p3, width=8,
         height=max(4, nrow(robust)*0.35 + 2), dpi=300)
}

# Plot the corresponding figure.
selMat <- matrix(0L, nrow=length(topGenes), ncol=nrow(spRuns),
                 dimnames=list(topGenes, spRuns$runId))
for(i in seq_len(nrow(spLong))){
  if(spLong$selected[i] && spLong$gene[i] %in% topGenes)
    selMat[spLong$gene[i], spLong$runId[i]] <- 1L
}
pdf(file.path(outRoot, "4_selection_heatmap.pdf"), width=9, height=7)
pheatmap(selMat, color=c("grey92", "#d7301f"), cluster_cols=FALSE, cluster_rows=FALSE,
         show_colnames=FALSE, fontsize_row=8, border_color=NA,
         legend_breaks=c(0,1), legend_labels=c("Not selected","p<0.05"),
         main="Cox selection status of top genes across 50 runs (red = selected)")
dev.off()
png(file.path(outRoot, "4_selection_heatmap.png"), width=9, height=7, units="in", res=300)
pheatmap(selMat, color=c("grey92", "#d7301f"), cluster_cols=FALSE, cluster_rows=FALSE,
         show_colnames=FALSE, fontsize_row=8, border_color=NA,
         legend_breaks=c(0,1), legend_labels=c("Not selected","p<0.05"),
         main="Cox selection status of top genes across 50 runs (red = selected)")
dev.off()

# Cox regression analysis step.
p5 <- runInfo %>%
  filter(mode == "split") %>%
  mutate(runIdx=as.integer(str_extract(runId, "\\d+"))) %>%
  pivot_longer(c(nEvents, nGrey, nDEG, nInter, nCox), names_to="metric", values_to="value") %>%
  mutate(metric=factor(metric, levels=c("nEvents","nGrey","nDEG","nInter","nCox"),
                       labels=c("OS events","grey genes","DEG count","Intersection count","Cox selected"))) %>%
  ggplot(aes(x=runIdx, y=value)) +
  geom_col(fill="#5690b5", width=0.7) +
  facet_wrap(~metric, ncol=1, scales="free_y") +
  labs(x="Run index", y=NULL, title="Per-run summary over 50 runs") + theme_pub
ggsave(file.path(outRoot, "5_run_summary.pdf"), p5, width=8, height=9)
ggsave(file.path(outRoot, "5_run_summary.png"), p5, width=8, height=9, dpi=300)

cat(sprintf("\n[完成] 结果输出至：%s\n", outRoot))
