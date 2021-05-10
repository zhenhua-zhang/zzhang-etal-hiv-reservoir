#!/bin/bash
# set -Eeu -o pipefail


erro() {
    echo -e "[E]: $* Exiting..." >&2
}

warn() {
    echo -e "[W]: $*" >&2
}

info() {
    echo -e "[I]: $*" >&2
}

echoUsage() {
    cat <<EOF
Usage:
    ./multi_traits_by_gemma.sh [OPTIONS]
EOF
}

echoHelp() {
    cat <<EOF
multi_traits_by_gemma.sh, a wrapper for GEMMA

Usage:
    ./multi_traits_by_gemma.sh [OPTIONS]

Help:
    -h, --help  Print this help and exit.

Example:
    # 1. Show help message
    ./multi_traits_by_gemma.sh -h
    ./multi_traits_by_gemma.sh --help
EOF
}

opts=$(getopt -l "test-model:,usage,help" -o "m:,u,h" -a -- "$@")
eval set -- "$opts"
while true; do
    case "$1" in
        -h|--help)
            echoHelp
            exit 0 ;;
        -u|--usage)
            echoUsage
            exit 0 ;;
        -m|--test-model)
            case $2 in
                lm|lmm|bslmm)
                    testModel=$2 ;;
                *)
                    warn 'Unknown or unsupported test-model, using default instead.'
                    testModel=ulm ;;
            esac
            shift ;;
        --)
            shift && break ;;
        *)
            erro "Unknown sub-command '$1'"
            echoUsage
            exit 0 ;;
    esac
    shift
done

#
## If not given, about to use the default.
#
testModel=${testModel:=ulm}

currentHost=$(hostname)
if [[ $currentHost =~ umcg || $currentHost =~ calculon || $currentHost =~ boxy ]]; then
    module load BCFtools 1>&2 2>/dev/null
fi

vcfpath=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/wd_hiv_genotype/annotated
fmtstring='%ID, %ALT, %REF, [%DS, ]\n'
headerFmtRegex='s/(\[.{1,3}\])|:DS|v_.{1,4}_//g'   # Sed regex string

#
## BCFtools::Genotypes
#
bimbamFile=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/wd_hiv_genotype/dosage/GEMMA/200hiv-genotype.bimbam
if [[ ! -e $bimbamFile ]]; then
    echo "Generating genotype data in BIMBAM format"
    count=0
    for vcfFile in "$vcfpath"/*.vcf.gz; do
        if [[ $count -eq 0 ]]; then
            bcftools query -H -f "$fmtstring" "$vcfFile" \
                | sed -E -e 's/, $//g' -e "$headerFmtRegex" \
                >> "$bimbamFile"
        else
            bcftools query -f "$fmtstring" "$vcfFile" \
                | sed -e 's/, $//g' \
                >> "$bimbamFile"
        fi

        (( count++ ))
    done
else
    echo "Found genotype data, skipping..."
fi

annotationFile=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/wd_hiv_genotype/dosage/GEMMA/200hiv-annotation.tsv
if [[ ! -e $annotationFile ]]; then
    echo "Generating SNP annotation for genotype data in BIMBAM format"
    for vcfFile in "$vcfpath"/*.vcf.gz; do
        bcftools query -f '%ID, %POS, %CHROM\n' "$vcfFile" >> "$annotationFile"
    done
fi

rawPhenotypeFile=~/Documents/projects/wp_hiv_reservoir/inputs/phenotype/hiv_reservoir/20190524_HIVreservoir_GENT_withRNADNARatio.csv
rawCovariateFile=~/Documents/projects/wp_hiv_reservoir/inputs/phenotype/trait/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv
phenotypeFile=200hiv-phenotype.tsv
covariateFile=200hiv-covariate.tsv

if [[ -e $phenotypeFile ]]; then rm $phenotypeFile; fi
if [[ -e $covariateFile ]]; then rm $covariateFile; fi
for sampleId in $(grep --max-count 1 ID $bimbamFile | sed 's/, /\n/g' | grep ^X); do
    # Covariates
    hasCovariate=$(grep -c --max-count 1 "$sampleId" $rawCovariateFile)
    if [[ $hasCovariate -eq 1 ]]; then
        grep -w "$sampleId" $rawCovariateFile \
            | awk -F',' '{print $4"\t"$5"\t"$13"\t"$14}' \
            | sed -e 's/\t\t/\tNA\t/g' -e 's/^\t/NA\t/g' -e 's/\t$/\tNA/g' \
            >> $covariateFile
    else
        echo -e "NA\tNA\tNA\tNA" >> $covariateFile
    fi

    # Phenotype
    hasPhenotype=$(grep -c --max-count 1 "$sampleId" $rawPhenotypeFile)
    if [[ $hasPhenotype -eq 1 ]]; then
        grep -w "$sampleId" $rawPhenotypeFile \
            | awk -F',' '{print $3"\t"$5"\t"$7}' \
            | sed -e 's/\t\t/\tNA\t/g' -e 's/^\t/NA\t/g' -e 's/\t$/\tNA/g' \
            >> $phenotypeFile
    else
        echo -e "NA\tNA\tNA" >> $phenotypeFile
    fi
done


#
## GEMMA::RelatenessMatrix
#
if [[ $currentHost =~ umcg || $currentHost =~ calculon || $currentHost =~ boxy ]]; then
    module load GEMMA 1>&2 2>/dev/null
fi

# After generating the relateness file, a ".cXX.txt" will be added.
relatenessFile=200hiv-relateness
if [[ ! -e $relatenessFile".cXX.txt" ]]; then
    echo "Generating relateness data..."
    gemma -g $bimbamFile -p $phenotypeFile -gk -o $relatenessFile
    mv output/$relatenessFile".cXX.txt" ./
else
    echo "Found relateness data, skipping..."
fi

#
## GEMMA::MultivariateLinearMixedModel
#
outputPrefix=200hiv-qtlm_gemma-$testModel
case $testModel in
    lm)
        echo "Quatitive traits locus mapping analysis by univariate linear model ..."
        gemma -g $bimbamFile -p $phenotypeFile -a $annotationFile -c $covariateFile \
            -lm 4 \
            -o $outputPrefix
        ;;
    lmm)
        echo "Quatitive traits locus mapping analysis by multivariate linear mixed model ..."
        gemma -g $bimbamFile -p $phenotypeFile -a $annotationFile -c $covariateFile \
            -k $relatenessFile".cXX.txt" \
            -lmm 4 \
            -n 1 2 \
            -o $outputPrefix
        ;;
    bslmm)
        echo "Quatitive traits locus mapping analysis by Bayesian sparse linear mixed model ..."
        gemma -g $bimbamFile -p $phenotypeFile -a $annotationFile -c $covariateFile \
            -bslmm 1 \
            -o $outputPrefix
        ;;
    *)
        erro "Unknown association test model."
esac

# vim: set ft=bash:
