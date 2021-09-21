#!/bin/bash
#SBATCH --mem 10G
#SBATCH --time 0:59:0
#SBATCH --output /groups/umcg-wijmenga/tmp01/users/umcg-zzhang/projects/wp_hiv_reservoir/outputs/bmc_revision/ver_alpha/cc_assoc-%j.log
#SBATCH --job-name plink_hiv_res_bmc
#SBATCH --cpus-per-task 24

# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 18, 2021
# Updated: Aug 18, 2021

# Estimate the genetic effects by splitting the cohort into two groups by the
# RNA:DNA ratio. 

# R1Q1:
# According to table 1, from 207 individuals, 201 presented a plasma viral
# load <40 cps/ml - thus, 6 individuals presented viral load superior to
# 40 cps/ml and their RNA:DNA ratio is informative or suggestive of active
# replication

# Basic settings
set -Eeu -o pipefail
wkdir=~/Documents/projects/wp_hiv_reservoir

# We select 207 individuals used in the original analysis and group them into
# two groups (6 vs 201), this is done by editing .fam file.

# 1. Association analysis by plink 1.9
module load PLINK/1.9-beta6-20190617

## 1.1 Version
plink --version
# PLINK v1.90b6.10 64-bit (17 Jun 2019)

## 1.2 Case/control (splitted by viral load 40 copies/ml) association analysis
plink \
  --bfile $wkdir/inputs/genotype/plink/200hiv_filtered \
  --keep $wkdir/outputs/bmc_revision/inputs/samples_to_proc.csv \
  --pheno $wkdir/outputs/bmc_revision/inputs/cc_pheno.txt \
  --mpheno 2 \
  --covar $wkdir/outputs/bmc_revision/inputs/hiv-reservoir.proc_cvrt.csv \
  --covar-name AGE,CD4_NADIR,HIV_DURATION \
  --logistic mperm=5000 \
  --adjust qq-plot \
  --ci 0.95 \
  --parameters 1 \
  --mperm-save \
  --threads 10 \
  --seed 131072 \
  --out $wkdir/outputs/bmc_revision/ver_alpha/200hiv_vl_cc

