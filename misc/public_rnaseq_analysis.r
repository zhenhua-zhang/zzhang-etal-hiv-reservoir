#!/usr/bin/env Rscript


if (!require(DESeq2)) {
  install.packages("DESeq2")
}


#' Do a differential expression analysis using DESeq2 package
#' @param rcm read counts matrix
#' @param dsm design matrix
#' @param desform design formula
#' @param tagsmp samples will be analysized
#' @param contvec contrast vector
#' @param lfccoef lfc coefficient
#' @param mincounts minimum counts
#' @return a list including raw results, shrinked results, and deseq matrix.
#' @example 
#' diffExpr(readCountMatrix, designMatrix)
diffExpr <- function(rcm, dsm, desform=NULL, tagsmp=NULL, contvec=NULL, lfccoef=NULL, lfctype="apeglm", mincounts=10) {
  rcm <- rcm[, tagsmp] # load read counts
  dsm <- dsm[tagsmp, ]
  dds <- DESeqDataSetFromMatrix(countData = rcm, colData = dsm, design = desform) # Make deseq data set
  dds <- dds[which(rowSums(counts(dds)) >= mincounts), ] # Filtering
  dds <- DESeq(dds) # Differential expression analysis
  res <- results(dds, contrast=contvec)
  lfcRes <- lfcShrink(dds, coef=lfccoef, type=lfctype)
  
  return(list(lfcRes=lfcRes, rawRes=res, dds=dds))
}

mklfcCeof <- function(vec) {
  return(paste(vec[1], vec[2], "vs", vec[3], sep="_"))
}



# The data are from a RNA-seq study by title:
# Deep sequencing of HIV-infected cells: insights into nascent transcription and
# host-directed therapy

# The GEO accessiont number is GSE53993

setwd("~/Documents/projects/wp_hiv_reservoir/inputs/publicRNASeqResults/GSE53993")
dsm <- read.csv("/GSE53993_Total_RNA_RPKM_Matrix_sample_annotation.txt", sep="\t", row.names = 1)
dsm$condition <- factor(dsm$condition, levels = c("control_H12", "treated_H12", "control_H24", "treated_H24"))

rcm <- read.csv("/GSE53993_Total_RNA_RPKM_Matrix.txt", sep="\t", row.names = "ID")
target_samples <- base::sapply(rownames(dsm), function(e) {paste0("X", e)})
contrast_vec <- c("condition", "control_H12", "treated_H12")
lfccoef <- "condition_treated_H12_vs_control_H12"
desform <- formula("~ condition")

dea <- diffExpr(rcm, dsm, desform = desform, tagsmp = target_samples, contvec = contrast_vec, lfccoef = lfccoef)

write.csv(as.data.frame(dea$rawRes[order(dea$rawRes$pvalue), ]), file="condition_treated_H12_vs_control_H12.csv")
plotMA(dea$lfcRes, ylim=c(-2, 2))



# meaningless without repeat.
# Second data set are from RNA-seq study by title:
# Dynamics of HIV Latency and Reactivation in a Primary CD4+ T Cell Model 
# The GEO accession number is GSE95297
setwd("~/Documents/projects/wp_hiv_reservoir/inputs/publicRNASeqResults/GSE95297")
dsm <- read.csv("GSE95297_sample_annotation.txt", sep="\t", row.names = 1)
rcm <- read.csv("GSE95297_rawcounts.txt", sep="\t", row.names = "geneID")

desform <- formula("~ condition")

# target_samples <- c(paste0("HIVW", c("0", "2")), paste0("mockW", c("0", "2")))
# contrast_vec <- c("condition", "mockW02", "HIVW02")
target_samples <- c(paste0("HIVW", c("4", "6")), paste0("mockW", c("4", "6")))
contrast_vec <- c("condition", "mockW46", "HIVW46")
# target_samples <- c(paste0("HIVW", c("8", "10")), paste0("mockW", c("8", "10")))
# contrast_vec <- c("condition", "mockW810", "HIVW810")

lfccoef <- mklfcCeof(contrast_vec)
dea <- diffExpr(rcm, dsm, desform=desform, tagsmp=target_samples, contvec=contrast_vec, lfccoef=lfccoef, lfctype="ashr")

plotMA(dea$rawRes, ylim=c(-10, 10))
write.csv(as.data.frame(dea$rawRes[order(dea$rawRes$pvalue), ]), file=paste0(lfccoef, ".csv"))
