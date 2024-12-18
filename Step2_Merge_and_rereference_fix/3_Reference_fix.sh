#!/bin/bash

########## GERMAN #########

# Remove identical or closely related individuals
plink2 --bfile merged_GOFA_KFO --remove merged_GOFA_KFO.fail-IBD.IDs --make-bed --out clean_GOFA_KFO

# Generate and reformat the frequency file
plink2 --bfile clean_GOFA_KFO --freq --out clean_GOFA_KFO_freq
awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' clean_GOFA_KFO_freq.afreq > clean_GOFA_KFO_freq.frq

# Fix merged file to reference
perl $perlFile -b clean_GOFA_KFO.bim -f clean_GOFA_KFO_freq.frq -r $HRC -h -o fix_HRC

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' fix_HRC/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' fix_HRC/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' fix_HRC/Run-plink.sh
source fix_HRC/Run-plink.sh

# Prepare zipped VCF file to impute
plink2 --bfile fix_HRC/clean_GOFA_KFO-updated --recode vcf --make-bed --out fix_HRC/clean_GOFA_KFO_final
pbgzip -c fix_HRC/clean_GOFA_KFO_final.vcf > fix_HRC/clean_GOFA_KFO_final.vcf.gz


########## SPAIN: BCN + GCAT #########

# Remove identical or closely related individuals
plink2 --bfile merged_BCN_GCAT --remove merged_BCN_GCAT.fail-IBD.IDs --make-bed --out clean_BCN_GCAT

# Generate and reformat the frequency file
plink2 --bfile clean_BCN_GCAT --freq --out clean_BCN_GCAT_freq
awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' clean_BCN_GCAT_freq.afreq > clean_BCN_GCAT_freq.frq

# Fix merged file to reference
perl $perlFile -b clean_BCN_GCAT.bim -f clean_BCN_GCAT_freq.frq -r $HRC -h -o fix_HRC_Bar

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' fix_HRC_Bar/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' fix_HRC_Bar/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' fix_HRC_Bar/Run-plink.sh
source fix_HRC_Bar/Run-plink.sh

# Prepare zipped VCF file to impute
plink2 --bfile fix_HRC_Bar/clean_BCN_GCAT-updated --recode vcf --make-bed --out fix_HRC_Bar/clean_BCN_GCAT_final
pbgzip -c fix_HRC_Bar/clean_BCN_GCAT_final.vcf > fix_HRC_Bar/clean_BCN_GCAT_final.vcf.gz


########## ALL: GOFA_KFO + BCN_GCAT #########

# Remove identical or closely related individuals
plink2 --bfile merged_ALL --remove merged_ALL.fail-IBD.IDs --make-bed --out clean_ALL

# Generate and reformat the frequency file
plink2 --bfile clean_ALL --freq --out clean_ALL_freq
awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' clean_ALL_freq.afreq > clean_ALL_freq.frq

# Fix merged file to reference
perl $perlFile -b clean_ALL.bim -f clean_ALL_freq.frq -r $HRC -h -o fix_HRC_ALL

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' fix_HRC_ALL/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' fix_HRC_ALL/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' fix_HRC_ALL/Run-plink.sh
source fix_HRC_ALL/Run-plink.sh

# Prepare zipped VCF file to impute
plink2 --bfile fix_HRC_ALL/clean_ALL-updated --recode vcf --make-bed --out fix_HRC_ALL/clean_ALL_final
pbgzip -c fix_HRC_ALL/clean_ALL_final.vcf > fix_HRC_ALL/clean_ALL_final.vcf.gz
