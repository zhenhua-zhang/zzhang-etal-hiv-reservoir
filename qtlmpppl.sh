#!/bin/bash
set -Ee -o pipefail

# hiv reservoir
export OUTDATED_IGNORE=1

# pjdir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/wp_hiv_reservoir
pjdir=${HOME}/Documents/projects/wp_hiv_reservoir

hstnm=$(hostname)
if [[ ${hstnm} =~ calculon ]] || [[ ${hstnm} =~ umcg-node ]]; then
    source /apps/modules/modules.bashrc
    module purge
    module load Python/3.6.3-foss-2015b
    source "${pjdir}/scripts/.env/bin/activate"
fi

mkdir -p "${pjdir}"/{preprocess,qtlmapping/{boxplot,per_trait}}

# Preprocess
preproc=0
[[ ${preproc} -eq 1 ]] \
    && python ../qtlmpppl/src/preprocess.py -c ../configs/hiv_reservoir.json \
    || echo "Skip pre-process."

# QTL mapping
qtlmp=0
if [[ ${qtlmp} -eq 1 ]]; then
    [[ ${hstnm} =~ calculon || ${hstnm} =~ umcg-node ]] \
        && module purge \
        && module load R/3.5.1-foss-2015b-bare

    Rscript ../qtlmpppl/src/matrix_eqtl.r \
    --gntp-file "${pjdir}"/outputs/genotype/dosage/MatrixEQTL/200HIV_dosage.gz \
    --gntp-info-file "${pjdir}"/outputs/genotype/dosage/MatrixEQTL/200HIV_variantInfo.gz \
    --phtp-file "${pjdir}"/outputs/hiv_res/preprocess/hiv-reservoir.proc_phtp.tsv \
    --cvrt-file "${pjdir}"/outputs/hiv_res/preprocess/hiv-reservoir.proc_cvrt.tsv \
    --save-pref "${pjdir}"/outputs/hiv_res/qtlmapping/
else
    echo "Skip QTL mapping."
fi


[[ ${hstnm} =~ calculon || ${hstnm} =~ umcg-node ]] \
    && module purge \
    && module load Python/3.6.3-foss-2015b \
    && source "${pjdir}"/scripts/.env/bin/activate

# Boxplot
test_run=0
if [[ ${test_run} -eq 1 ]]; then
    dosage_file=${pjdir}/inputs/genotype/dosage/MatrixEQTL/per_chrom/chr11_dosage.gz
    info_file=${pjdir}/inputs/genotype/dosage/MatrixEQTL/per_chrom/chr11_variantInfo.gz
else
    dosage_file=${pjdir}/inputs/genotype/dosage/MatrixEQTL/200HIV_dosage.gz
    info_file=${pjdir}/inputs/genotype/dosage/MatrixEQTL/200HIV_variantInfo.gz
fi

draw_bp=1
[[ ${draw_bp} -eq 1 ]] \
    && python ../qtlmpppl/src/boxplot.py \
    -p "${pjdir}"/outputs/hiv_res/preprocess/hiv-reservoir.proc_phtp.tsv \
    -c "${pjdir}"/outputs/hiv_res/preprocess/hiv-reservoir.proc_cvrt.tsv \
    -g "${dosage_file}" \
    -i "${info_file}" \
    -P RNAHIV_CD4LOG DNAHIV_CD4LOG RNAvsDNA_CD4LOG \
    -G rs7817589 rs7113255 rs7113204 rs2613996 rs1015164 rs12366210 \
    -o "${pjdir}"/outputs/hiv_res/qtlmapping/boxplot \
    || echo "Skip boxplot."
