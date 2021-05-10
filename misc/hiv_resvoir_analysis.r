#!/usr/bin/env Rscript

library(dplyr)
library(Hmisc)
library(ggpubr)
library(ggpmisc)
library(ggplot2)
library(data.table)

pjdir <- "~/Documents/projects/wp_hiv_reservoir"
setwd(paste0(pjdir, "/outputs/hiv_res/preprocess"))

cvrt_flnm <- "hiv-reservoir.proc_cvrt.tsv"
cvrt_dttb <- fread(cvrt_flnm)

phtp_flnm <- "hiv-reservoir.proc_phtp.tsv"
phtp_dttb <- fread(phtp_flnm)

name_vec <- c("gender", "age", "CD4_NADIR", "HIV_DURATION", "RNAHIV_CD4LOG", "DNAHIV_CD4LOG", "RNAvsDNA_CD4LOG")
names(name_vec) <- paste0("V", 1:7)
trait_dtfm <- rbindlist(list(cvrt_dttb, phtp_dttb)) %>%
  select(starts_with("X")) %>%
  transpose() %>%
  rename_with(function(e) {return(name_vec[e])}) %>%
  mutate(classes=1, gender=as.factor(gender)) %>%
  filter(!is.na(gender)) %>%
  as.data.frame()

#
## Plot boxplot per gender
#
# wk_trait <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "gender")
# d_bxpt <- reshape2::melt(trait_dtfm[, wk_trait])
# 
# dp_name <- c("CA HIV DNA", "CA HIV RNA", "RNA:DNA ratio", "Gender")
# names(dp_name) <- wk_trait
# 
# TODO: This will give the p-values
# g_bxpl <- ggboxplot(d_bxpt, x="gender", y="value", color="gender", add="jitter", facet.by="variable") +
#   stat_compare_means(label="p.format")

# g_bxpl <- ggplot(data=d_bxpt) + theme_bw() +
#   geom_boxplot(aes(x=gender, y=value, fill=gender)) +
#   geom_jitter(aes(x=gender, y=value), alpha=0.5, size=1) +
#   compare_means() +
#   stat_compare_means(method="t.test") +
#   facet_wrap(~variable, nrow=1, scales="free", labeller=labeller(variable=dp_name)) +
#   scale_fill_discrete(name="Gender", labels=c("Female", "Male")) +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.background=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.x=element_blank()) +
#   labs(x="", y="Trait level (log10)")
# 
# ggsave("Sup-fig1.png", plot=g_bxpl)
# ggsave("Sup-fig1.pdf", plot=g_bxpl)

#
## plot violin plot
#
wk_trait <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "age", "CD4_NADIR", "HIV_DURATION")
d_bxpt <- reshape2::melt(trait_dtfm[, wk_trait])

dp_name <- c("HIV-1 DNA", "HIV-1 RNA", "Ratio", "Age", "CD4 nadir", "HIV dur.")
names(dp_name) <- wk_trait

g_bxpl <- ggplot(data=d_bxpt) + theme_bw() +
  geom_violin(aes(x=1, y=value)) +
  geom_jitter(aes(x=1, y=value), alpha=0.5, size=1) +
  facet_wrap(~variable, nrow=1, scales="free", labeller=labeller(variable=dp_name)) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="", y="")


#
## Correlation analysis by Spearman's rho
#
wk_trait <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "age", "CD4_NADIR", "HIV_DURATION")
d_bxpt <- reshape2::melt(trait_dtfm[, wk_trait])

dp_name <- c("CA HIV-1 DNA", "CA HIV-1 RNA", "RNA:DNA ratio", "Age", "CD4 nadir", "HIV-1 duration")
names(dp_name) <- wk_trait

cor_mat <- rcorr(as.matrix(trait_dtfm[, wk_trait]), type="spearman")
cor_p_mat <- cor_mat$P
cor_p_mat[upper.tri(cor_p_mat)] <- NA

cor_r_mat <- cor_mat$r
cor_r_mat[lower.tri(cor_r_mat)] <- NA

cor_p_mat <- as.data.frame(reshape2::melt(cor_p_mat, na.rm=TRUE))
cor_r_mat <- as.data.frame(reshape2::melt(cor_r_mat, na.rm=TRUE))
cor_p_mat$value <- p.adjust(cor_p_mat$value, method="fdr", n=15)

myfmt <- function(v) {
  ifelse(v<1e-16, "<1e-16", ifelse(v<0.001, format(v, digits=3, scientific=T), round(v, 3)))
}

g_crpl <- ggplot() + theme_minimal() +
  geom_tile(aes(x=Var1, y=Var2, fill=value), data=cor_r_mat) +
  geom_text(aes(x=Var1, y=Var2, label=round(value, 3)), data=cor_r_mat, color="black", size=4) +
  geom_text(aes(x=Var1, y=Var2, label=myfmt(value)), data=cor_p_mat, color="black", size=4) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limit=c(-1, 1), name="Spearman's Rho") +
  scale_x_discrete(labels=dp_name) +
  scale_y_discrete(labels=dp_name) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        legend.position="top",
        legend.direction="horizontal",
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=90, hjust=0.99, vjust=0.1, size=12),
        axis.text.y=element_text(size=12))


#
## CA HIV RNA v.s. total HIV DNA
#
rna_vs_dna <- c("RNAHIV_CD4LOG", "DNAHIV_CD4LOG")
g_pnpl <- ggplot(aes(x=DNAHIV_CD4LOG, y=RNAHIV_CD4LOG), data=trait_dtfm[, rna_vs_dna]) + 
  theme_classic() +
  geom_point() +
  geom_smooth(method='lm', formula = y~x) +
  stat_cor() +
  labs(x="CA HIV-1 DNA (log10)", y="CA HIV-1 RNA (log10)") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

ggarrange(NULL, g_bxpl, g_crpl, g_pnpl, nrow=2, ncol=2) %>%
  ggexport(filename="Fig1.png", width=4100, height=4000, res=300)

ggarrange(NULL, g_bxpl, g_crpl, g_pnpl, nrow=2, ncol=2) %>%
  ggexport(filename="Fig1.pdf", width=13, height=12, res=300)
