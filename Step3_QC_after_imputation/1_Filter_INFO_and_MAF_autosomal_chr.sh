#!/bin/bash

# Select thresholds
INFO_HIGH=0.8 # for rare variants (MAF <= MAF_HIGH)
INFO_LOW=0.5 # for common variants (MAF > MAF_HIGH)
HWE=1e-12 # only in controls
MAF_LOW=0.001
MAF_HIGH=0.01

IDs_controls="PATH_TO_FILE_CONTAINING_ALL_CONTROL_IDS"

datasets=("GERMAN" "SPAIN" "ALL")

for dataset in "${datasets[@]}"; do

	input_dir="$dataset.vcfs"  # folder downloaded from Sanger Imputation Server
	output_dir="$dataset.INFO${INFO_LOW}_${INFO_HIGH}.MAF$MAF_LOW.HWE$HWE"
	temp_dir="$dataset.INFO${INFO_LOW}_${INFO_HIGH}.MAF$MAF_LOW.HWE$HWE/temp"

	mkdir -p $output_dir
	mkdir -p $temp_dir

	for chr in {1..22}; do

	    # Filter INFO and MAF for all variants
	    bcftools view -i "(INFO>=$INFO_HIGH && AC >= $MAF_LOW * AN && AC < $MAF_HIGH * AN) || (INFO>=$INFO_LOW && AC >= $MAF_HIGH * AN)" $input_dir/$chr.vcf.gz | \
	    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' -Ov -o $temp_dir/$chr.TEMP1.vcf

	    # Filter HWE in controls on the combined VCF
	    bcftools +fill-tags $temp_dir/$chr.TEMP1.vcf -- -S $IDs_controls -t HWE | \
	    bcftools view -i "HWE >= $HWE" -Ov -o $temp_dir/$chr.TEMP2.vcf
	    
	    # Keep variants that pass the HWE filter from the combined filtered set
	    vcftools --vcf $temp_dir/$chr.TEMP1.vcf --positions $temp_dir/$chr.TEMP2.vcf --recode --recode-INFO-all --stdout | \
	    bgzip -c > $output_dir/chr$chr.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf.gz

	    # Create index file
	    tabix --csi -p vcf $output_dir/chr$chr.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf.gz

	done

	# Concatenate all chromosomes
	bcftools concat $output_dir/chr{1..22}.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf.gz -Oz -o $output_dir/$dataset.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf.gz

	# Generate PLINK files
	plink2 --vcf $output_dir/$dataset.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE.recode.vcf.gz --make-bed --out $output_dir/$dataset.MAF$MAF_LOW.INFO${INFO_LOW}_${INFO_HIGH}.HWE$HWE

	rm -r $temp_dir

done
