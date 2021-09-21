#!/usr/bin/env Rscript
library(qqman)
library(tidyverse)
library(data.table)

gwas_ss_file <- "../../outputs/bmc_revision/ver_alpha/200hiv_vl_cc.assoc.logistic"
gwas_ss_tab <- fread(gwas_ss_file, verbose=F, showProgress=F)

pdf("../../outputs/bmc_revision/ver_alpha/200hiv_vl_cc.manhattan_plot.pdf",
    width=16, height=9, onefile=F)
gwas_ss_tab %>%
  dplyr::filter(P<0.1) %>%
  qqman::manhattan(chr="CHR", bp="BP", p="P")
dev.off()

pdf("../../outputs/bmc_revision/ver_alpha/200hiv_vl_cc.qq_plot.pdf",
    width=7, height=7, onefile=F)
gwas_ss_tab %$% P %>%
  qqman::qq()
dev.off()
