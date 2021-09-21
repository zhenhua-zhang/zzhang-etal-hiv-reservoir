#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 23, 2021
# Updated: Aug 23, 2021
options(error=traceback)

# Estimate the DE genes in development of HIV latency.

library(DESeq2)
library(magrittr)
library(ggsignif)
library(tidyverse)
library(data.table)

# Workspace
wkdir <- "../../outputs/bmc_revision/ver_alpha/pub_rnaseq"
rnaseq_file <- str_glue("../../inputs/publicRNASeqResults/GSE127468/GSE127468_ExpressionMatrix_counts.csv")

# Load read counts.
rnaseq_tab <- data.table::fread(rnaseq_file) %>%
  dplyr::mutate(GeneID=str_remove(GeneID, pattern="\\.[0-9]+$")) # Remove isoform information.

ensg_sym_map <- c("ENSG00000174915"="PTDSS2", "ENSG00000185507"="IRF7",
                  "ENSG00000177030"="DEAF1",  "ENSG00000023191"="RNH1")

tp_set <- paste0(rep("day", 4), c(3, 7, 9, 14))
deg_est <- NULL
for (cur_tp in tp_set) {
  # Construct a data.frame with RNA-seq samples information, e.g. condition, time-point.
  samples <- colnames(rnaseq_tab) %>%
    tibble(.name_repair=function(e) {if_else(e==".", "Samples", "Uknown")}) %>%
    dplyr::filter(!(Samples %in% c("GeneID"))) %>%
    tidyr::separate(Samples, into=c("Sample", "Condition", "TimePoint"), sep="_", remove=F) %>%
    as.data.frame() %>%
    (function(x) {rownames(x) <- x$Samples; x}) %>%
    dplyr::select(-Samples) %>%
    dplyr::mutate(Condition=factor(Condition, levels=c("Infected", "Mock")),
                  TimePoint=factor(TimePoint, levels=tp_set)) %>%
    dplyr::filter(TimePoint==cur_tp)

  # Convert the read counts table into a matrix with colnames (samples) and rownames (gene).
  rnaseq_mat <- rnaseq_tab %>%
    dplyr::select(c(ends_with(cur_tp), GeneID)) %>%
    as.data.frame() %>%
    (function(x) {rownames(x) <- x$GeneID; x}) %>%
    dplyr::select(-one_of("GeneID")) %>%
    as.matrix()

  # Construct a DESeqDataSet object.
  dds <- DESeqDataSetFromMatrix(rnaseq_mat, colData=samples, design=~Condition)

  # Remove some low quality genes.
  keep <- base::rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq(dds)
  res <- results(dds, contrast=c("Condition", "Infected", "Mock"))
  deg_est <- res[names(ensg_sym_map), ] %>%
    as.data.frame() %>%
    dplyr::mutate(ensgid=rownames(.),
                  gene_sym=ensg_sym_map[ensgid],
                  Day=cur_tp,
                  exp={
                    vst(dds, blind=F) %>%
                      assay() %>%
                      (function(x) {
                         x <- x[ensgid, ]
                         colnames(x) <- colnames(x) %>%
                           str_remove(paste0("_", cur_tp)) %>% unlist()
                         return(as.data.table(x)) })
                  }) %>%
    as.data.table() %>%
    dplyr::bind_rows(deg_est)
}

deg_est %>%
  fwrite("../../outputs/bmc_revision/ver_alpha/pub_rnaseq/gene-exp_mock-vs-infected_de-est.csv",
         sep=",", row.names=F)

g_data <- deg_est %>%
  dplyr::mutate(fdr=p.adjust(pvalue, method="fdr")) %>%
  tidyr::pivot_longer(cols=starts_with("exp."),
                      names_to="Condition",
                      values_to="Expression") %>%
  dplyr::mutate(Condition=str_extract(Condition, "Mock|Infected"),
                Day=factor(Day, levels=c(tp_set))) %>%
  dplyr::relocate(ensgid, gene_sym, .before="baseMean") %>%
  dplyr::select(-padj)

g <- g_data %>%
  ggplot(aes(x=Day, y=Expression, color=Condition)) +
  geom_boxplot() +
  geom_point(size=1, position=position_dodge2(0.2)) +
  facet_grid(.~gene_sym) +
  theme_classic()

ggsave("../../outputs/bmc_revision/ver_alpha/pub_rnaseq/gene-exp_mock-vs-infected.pdf",
       plot=g, width=8, height=4)
