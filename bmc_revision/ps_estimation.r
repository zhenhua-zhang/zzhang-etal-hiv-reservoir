#!/usr/bin/env Rscript
options(error=traceback)
library(preprocessCore)
library(data.table)
library(tidyverse)
library(magrittr)
library(ggpubr)

ps_name <- c("PS(14:0/15:0)",
             "PS(15:0/20:3(5Z,8Z,11Z))",
             "PS(16:0/20:3(8Z,11Z,14Z))",
             "PS(18:3(9Z,12Z,15Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))",
             "PS(15:0/24:0)",
             "PS(DiMe(11,3)/DiMe(9,3))",
             "PS(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/24:1(15Z))")
ps_idx <- c(1936, 2007, 2020, 2047, 2050, 2055, 2104)
ps_hmdbid <- c("HMDB0112273", "HMDB0112329", "HMDB0012360", "HMDB0012417",
               "HMDB0112340", "HMDB0061555", "HMDB0112880")
names(ps_name) <- paste0("ionIdx_", ps_idx, "_log10")
names(ps_hmdbid) <- names(ps_name)


itfile_path <- "../../inputs/phenotype/metabolome"
ittab <- NULL
for (itfile in list.files(itfile_path, pattern="*.csv")) {
  ittab <- fread(str_glue("{itfile_path}/{itfile}")) %>%
    dplyr::bind_cols(ittab)
}


# Rename the tables and filter PS measurements.
# There are 174 individuals included in the metabolome results.
ps_ittab <- ittab %>%
  dplyr::select(-contains("V3"), ionIdx...1, ionMz...2, contains("200HIV")) %>%
  dplyr::rename_with(.fn=function(x) str_replace(x, " / 200HIV_[0-9] ", "_"),
                     .cols=contains("HIV")) %>%
  dplyr::rename("ionIdx"="ionIdx...1", "ionMz"="ionMz...2") %>%
  dplyr::filter(ionIdx %in% ps_idx) %>%
  tidyr::pivot_longer(cols=contains("_X"), values_to="intensity", names_to="individuals") %>%
  dplyr::mutate(individuals=str_extract(individuals, "X[0-9]{4,4}"),
                ionIdx=as.character(paste0("ionIdx_", ionIdx))) %>%
  dplyr::group_by(ionIdx, individuals) %>%
  dplyr::summarise(intensity=mean(intensity)) %>%
  tidyr::pivot_wider(id_cols=ionIdx, names_from=individuals, values_from=intensity) %>%
  dplyr::rename("mol_trait"="ionIdx") %>%
  dplyr::ungroup()


# Estimate the correlation between RNA:DNA ratio and PS
caxna_file <- "../../outputs/hiv_res/preprocess/hiv-reservoir.proc_phtp.tsv"
caxna_tab <- fread(caxna_file) %>%
  dplyr::rename("mol_trait"="V1")

#
chosen_cols <- intersect(colnames(caxna_tab), colnames(ps_ittab))
data_tab <- caxna_tab %>%
  dplyr::select(one_of(chosen_cols)) %>%
  dplyr::bind_rows(ps_ittab %>% dplyr::select(chosen_cols)) %>%
  tidyr::pivot_longer(cols=starts_with("X"), names_to="ind.", values_to="msm") %>%
  tidyr::pivot_wider(id_cols="ind.", names_from="mol_trait", values_from="msm") %>%
  (function(e) { # Quantiles normalization
     ion_qnorm <- e %>%
       dplyr::select(starts_with("ionIdx_")) %>%
       as.matrix() %>%
       normalize.quantiles() %>%
       as.data.frame()
     colnames(ion_qnorm) <- colnames(e)[5:11]
     cbind(e[, 1:4], ion_qnorm) }) %>%
  dplyr::mutate(ionIdx_1936_log10=log10(ionIdx_1936),
                ionIdx_2007_log10=log10(ionIdx_2007),
                ionIdx_2020_log10=log10(ionIdx_2020),
                ionIdx_2047_log10=log10(ionIdx_2047),
                ionIdx_2050_log10=log10(ionIdx_2050),
                ionIdx_2055_log10=log10(ionIdx_2055),
                ionIdx_2104_log10=log10(ionIdx_2104)) %>%
  dplyr::select(-c(ionIdx_1936, ionIdx_2007, ionIdx_2020, ionIdx_2047,
                   ionIdx_2050, ionIdx_2055, ionIdx_2104))

# Add genotype information
snp_file <- "../../outputs/bmc_revision/inputs/ps_snp.vcf"
snp_tab <- fread(snp_file, sep="\t") %>% 
  dplyr::select(-c(QUAL, FILTER, INFO, FORMAT)) %>%
  dplyr::rename("CHROM"="#CHROM") %>%
  dplyr::rename_with(.cols=starts_with("v_"), .fn=function(x) str_extract(x, "X[0-9]{4,4}")) %>%
  apply(2, FUN=function(x) {str_split(x, ":", simplify=T)[, 1]}) %>%
  as.data.table() %>%
  tidyr::pivot_longer(cols=starts_with("X"), names_to="ind.", values_to="genotype") %>%
  dplyr::mutate(genotype=case_when(genotype %in% c("0|0")        ~ paste0(REF, REF),
                                   genotype %in% c("1|0", "0|1") ~ paste0(REF, ALT),
                                   genotype %in% c("1|1")        ~ paste0(ALT, ALT)),
                mol_trait=ID,
                genotype=) %>%
  dplyr::select(-c(CHROM, ID, POS, REF, ALT)) %>%
  tidyr::pivot_wider(id_cols=ind., names_from=mol_trait, values_from=genotype)

data_tab <- dplyr::inner_join(snp_tab, data_tab, by="ind.")
%T>%
  fwrite("../../outputs/bmc_revision/ver_alpha/200hiv_metabolites/ps_vs_res-traits_data.csv",
         sep=",", row.names=F, quote=F)

cor_tab <- data_tab %>%
  dplyr::summarise(cor={
                     cur_colnames <- colnames(.)
                     res_trait <- cur_colnames[grepl("_CD4LOG", cur_colnames)]
                     ps_ionidx <- cur_colnames[grepl("_log10", cur_colnames)]
                     
                     cordf <- NULL
                     for (ionidx in ps_ionidx) {
                       for (trait in res_trait) {
                         cur_msm <- cur_data() %>%
                           dplyr::select(one_of(trait, ionidx)) %>%
                           dplyr::rename(c("trait"=all_of(trait),
                                           "ionidx"=all_of(ionidx))) %>%
                           as.data.table()

                         ctest <- cor.test(x=cur_msm$trait, y=cur_msm$ionidx,
                                           method="spearman", exact=F)
                         cordf %<>%
                           dplyr::bind_rows(data.frame(p_value=ctest$p.value,
                                                       rho=ctest$estimate,
                                                       trait=trait,
                                                       ionidx=ionidx,
                                                       ps_name=ps_name[ionidx],
                                                       hmdbid=ps_hmdbid[ionidx]))
                       }
                     }
                     cordf }) %>%
  as.data.table() %>%
  dplyr::rename_with(.fn=function(x) {str_remove(x, "cor.")})

cor_tab %>%
  fwrite("../../outputs/bmc_revision/ver_alpha/200hiv_metabolites/ps_vs_res-traits_spearman-cor.csv",
         sep=",", row.names=F, quote=F)

signif_ionidx <- cor_tab %>%
  dplyr::filter(trait=="RNAvsDNA_CD4LOG" & p_value<0.05) %$%
  ionidx

cur_data <- data_tab %>%
  dplyr::select(ind., starts_with("rs"), RNAvsDNA_CD4LOG, ends_with("_log10")) %>%
  tidyr::pivot_longer(cols=ends_with("_log10"), names_to="ionidx",
                      values_to="inten_log10") %>%
  dplyr::filter(ionidx %in% signif_ionidx) %>%
  dplyr::mutate(ps_hmdbid=ps_hmdbid[ionidx])

g_cor <- cur_data %>%
  ggplot(aes(x=RNAvsDNA_CD4LOG, y=inten_log10)) +
  geom_point() +
  geom_smooth(method="glm") +
  facet_grid(.~ps_hmdbid) +
  theme_classic() +
  labs(x=quote(~Log[2]~"(RNA:DNA ratio)"), y=quote(~Log[10]~"(PS intensity)"))

g_box_rs7 <- cur_data %>%
  dplyr::mutate(rs7113204=factor(rs7113204, levels=c("GG", "GA", "AA"))) %>%
  ggplot(aes(x=rs7113204, y=inten_log10)) +
  geom_boxplot() +
  facet_grid(.~ps_hmdbid) +
  theme_classic() +
  labs(x="Genotype of rs7113204", y=quote(~Log[10]~"(PS intensity)"))

g_box_rs2 <- cur_data %>%
  dplyr::mutate(rs2613996=factor(rs2613996, levels=c("AA", "AG", "GG"))) %>%
  ggplot(aes(x=rs2613996, y=inten_log10)) +
  geom_boxplot() +
  facet_grid(.~ps_hmdbid) +
  theme_classic() +
  labs(x="Genotype of rs2613996", y=quote(~Log[10]~"(PS intensity)"))

g <- ggarrange(g_cor, g_box_rs7, g_box_rs2, ncol=1, nrow=3, labels=c("A", "B", "C"))

ggsave(str_glue("../../outputs/bmc_revision/ver_alpha/200hiv_metabolites/ps_vs_res-traits_dotplot.pdf"),
         plot=g, width=7, height=10)
