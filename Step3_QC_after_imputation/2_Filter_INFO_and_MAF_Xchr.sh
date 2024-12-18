#!/bin/bash

# Select thresholds
INFO_HIGH=0.8 # for rare variants (MAF <= MAF_HIGH)
INFO_LOW=0.5 # for common variants (MAF > MAF_HIGH)
HWE=1e-12 # only in controls
MAF_LOW=0.001
MAF_HIGH=0.01

datasets=("GERMAN_X" "SPAIN_X" "ALL_X")

for dataset in "${datasets[@]}"; do

    input_dir="$dataset.vcfs"
    output_dir="$dataset.INFO${INFO_LOW}_${INFO_HIGH}.MAF$MAF_LOW.HWE$HWE"
    temp_dir="$dataset.INFO${INFO_LOW}_${INFO_HIGH}.MAF$MAF_LOW.HWE$HWE/temp"
    pheno="phenotype.txt"

    mkdir -p $output_dir
    mkdir -p $temp_dir

    # Prepare phenotype and IDs files
    awk 'BEGIN { OFS = "\t"; print "ID", "gender", "FA" } 
    NR > 1 {
        if ($5 == "original_n4230") {
            id = "0_" $7 "_" $7  # ID format for cohort1
        } else {
            id = $1 "_" $2  # ID format for other samples
        }
        print id, $3, $4  # Print ID, sex (column 3), and FA (column 4)
    }' $pheno > pheno.txt

    awk '{print $1}' pheno.txt > IDs.txt
    awk '$3== 0 {print $1"\tcontrols"}' pheno.txt > IDs_controls.txt

    # Filter INFO and MAF for all variants, remove pseudoautosomal regions (PARs)
    bcftools view -i "((INFO>=$INFO_HIGH && AC >= $MAF_LOW * AN && AC < $MAF_HIGH * AN) || (INFO>=$INFO_LOW && AC >= $MAF_HIGH * AN)) && AC != AN" -t ^X:60001-2699520,X:154931044-155260560 -S IDs.txt --force-samples $input_dir/X.vcf.gz | \
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' -Ov -o $temp_dir/X.TEMP1.vcf

    # Filter HWE in controls on the combined VCF
    bcftools +fill-tags $temp_dir/X.TEMP1.vcf -- -S IDs_controls.txt -t HWE | \
    bcftools view -i "HWE >= $HWE" -Ov -o $temp_dir/X.TEMP2.vcf

    # Keep variants that pass the HWE filter from the combined filtered set
    vcftools --vcf $temp_dir/X.TEMP1.vcf --positions $temp_dir/X.TEMP2.vcf --recode --recode-INFO-all --out $output_dir/chrX.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE

    plink2 --vcf $output_dir/chrX.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf --make-bed --out $output_dir/chrX.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE

    rm -r $temp_dir
done