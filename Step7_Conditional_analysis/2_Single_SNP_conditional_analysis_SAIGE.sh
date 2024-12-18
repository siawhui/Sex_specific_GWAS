#!/bin/bash

# File containing the output of GCTA (independent significant SNPs)
files=(conditional_analysis/results/*.jma.cojo)

# Loop through each file
for file in "${files[@]}"; do
    file_name=$(basename "$file")

    # Extract cohort
    if [[ "$file_name" == *"ALL"* ]]; then
        cohort="ALL"
    elif [[ "$file_name" == *"GOFA_KFO"* ]]; then
        cohort="GOFA_KFO"
    elif [[ "$file_name" == *"BCN_GCAT"* ]]; then
        cohort="BCN_GCAT"
    fi

    # Extract study name
    study="${file_name#${cohort}_}"
    study="${study%%.GWAS*}"
    
    # Print results (or use them in your script)
    echo "File: $file"
    echo "Cohort: $cohort"
    echo "Study: $study"
    echo "-----------"

    tail -n +2 "$file" | while IFS=$'\t' read -r col1 col2; do
        # Extract snp, chromosome and position from col2
        snp=$(echo "$col2" | awk '{print $1}')
        chr=$(echo "$col2" | cut -d':' -f1)
        pos=$(echo "$col2" | cut -d':' -f2)
        ref=$(echo "$col2" | cut -d':' -f3)
        alt=$(echo "$col2" | cut -d':' -f4 | cut -d$'\t' -f1)

        # Print the extracted values
        echo "SNP: $snp, Chromosome: $chr, Position: $pos, Ref: $ref, Alt: $alt"

        # Set up
        dir="GWAS/GWAS_${cohort}/${study}_conditional_analysis"
        vcfDir="${cohort}.INFO0.5_0.8.MAF0.001.HWE1e-12"
        rdaFile="GWAS/GWAS_${cohort}/${study}/${cohort}_${study}"

        # Create directory
        mkdir -p "$dir"

        outputPrefix="$dir/${chr}_${pos}_${ref}_${alt}"

        Rscript /fast/AG_Lee/Siaw/Sex_specific_GWAS/R_scripts/step2_SPAtests.R        \
            --vcfFile=$vcfDir/chr$chr.MAF0.001.INFO0.5_0.8.HWE1e-12.recode.vcf.gz       \
            --vcfFileIndex=$vcfDir/chr$chr.MAF0.001.INFO0.5_0.8.HWE1e-12.recode.vcf.gz.csi       \
            --vcfField=DS       \
            --SAIGEOutputFile=$outputPrefix.txt \
            --chrom=$chr       \
            --minMAF=0 \
            --minMAC=20 \
            --GMMATmodelFile=$rdaFile.rda \
            --varianceRatioFile=$rdaFile.varianceRatio.txt   \
            --LOCO=TRUE \
            --is_Firth_beta=TRUE    \
            --pCutoffforFirth=0.01 \
            --is_output_moreDetails=TRUE  \
            --condition=$snp

    done

done