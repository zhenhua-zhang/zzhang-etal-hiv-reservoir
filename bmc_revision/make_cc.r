#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 18, 2021
# Updated: Aug 23, 2021

library(magrittr)
library(tidyverse)
library(data.table)

# High viral load individuals (hvli)
hvli <- c("X1001", "X1003", "X1004", "X1022", "X1041", "X1209")
phtp_file <- "../../outputs/hiv_res/preprocess/hiv-reservoir.proc_phtp.tsv"
phtp_tab <- data.table::fread(phtp_file, verbose=F, showProgress=F)

# Minimum RNA:DNA ratio
rdr <- phtp_tab %>%
  dplyr::filter(V1=="RNAvsDNA_CD4LOG") %>%
  dplyr::select(all_of(hvli), -V1) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min=min(.), avg=mean(c_across()))
#    X1001  X1003  X1004 X1022 X1041 X1209   min    avg
#   -0.459 -0.616 -0.745 -1.02 -1.17 -1.01 -1.17 -0.885


# 1 is control, 2 is case.
cc_tab <- phtp_tab %>%
  dplyr::filter(V1=="RNAvsDNA_CD4LOG") %>%
  dplyr::select(-V1) %>%
  unlist() %>%
  (function(x) {
     data.table(FID=0, IID=names(x), rd_ratio=x) %>%
       dplyr::mutate(traits_avg=as.numeric(rd_ratio>=rdr$avg)+1,
                     traits_avg=case_when(is.na(traits_avg) ~ -9, T ~ traits_avg),
                     traits_min=as.numeric(rd_ratio>=rdr$min)+1,
                     traits_min=case_when(is.na(traits_min) ~ -9, T ~ traits_min))}) %T>%
  data.table::fwrite("../../outputs/bmc_revision/inputs/cc_pheno.txt", sep=" ",
                     row.names=F, quote=F)

cc_tab %>%
  summarise(traits_avg_n_ctrl=sum(traits_avg==1),
            traits_avg_n_case=sum(traits_avg==2),
            traits_min_n_ctrl=sum(traits_min==1),
            traits_min_n_case=sum(traits_min==2))
#   traits_avg_n_ctrl traits_avg_n_case traits_min_n_ctrl traits_min_n_case
#                 114                86               150                50
