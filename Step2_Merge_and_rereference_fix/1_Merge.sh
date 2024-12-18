#!/bin/bash

perlFile='PATH_TO_PERL_SCRIPT'
HRC='PATH_TO_HRC_REFERENCE'

########## GOFA (three batches) #########
# Merge n4230 with n285
plink --bfile GOFA_n4230 --bmerge GOFA_n285.bed GOFA_n285.bim GOFA_n285.fam --make-bed --out GOFA_n4230_n285

# Merge n95
plink --bfile GOFA_n4230_n285 --bmerge GOFA_n95.bed GOFA_n95.bim GOFA_n95.fam --make-bed --out GOFA_n4610

# Get rid of splited chromosome
plink2 --bfile GOFA_n4610 --make-pgen --sort-vars --out GOFA_n4610_sort

# Only SNPs that are genotyped on all arrays (carefully select --geno by reviewing the missingness file)
plink2 --pfile GOFA_n4610_sort --geno 0.01 --make-bed --out GOFA

########## GERMAN: GOFA + KFO #########

plink --bfile GOFA --bmerge KFO.bed KFO.bim KFO.fam --make-bed --out GOFA_KFO
plink2 --bfile GOFA_KFO --make-pgen --sort-vars --out GOFA_KFO_sort
plink2 --pfile GOFA_KFO_sort --make-bed --out GOFA_KFO_sort2
plink2 --bfile GOFA_KFO_sort2 --missing --out GOFA_KFO_sort2.missing
plink2 --bfile GOFA_KFO_sort2 --geno 0.001 --make-bed --out merged_GOFA_KFO  # IMPORTANT: carefully select --geno by reviewing the missingness file

########## SPAIN: BCN + GCAT #########

plink --bfile GCAT.clean-updated --bmerge BCN.bed BCN.bim BCN.fam --make-bed --out BCN_GCAT
plink2 --bfile BCN_GCAT --make-pgen --sort-vars --out BCN_GCAT_sort
plink2 --pfile BCN_GCAT_sort --make-bed --out BCN_GCAT_sort2
plink2 --bfile BCN_GCAT_sort2 --missing --out BCN_GCAT_sort2.missing
plink2 --bfile BCN_GCAT_sort2 --geno 0.02 --make-bed --out merged_BCN_GCAT  # IMPORTANT: carefully select --geno by reviewing the missingness file

########## ALL: GOFA_KFO + BCN_GCAT #########

plink --bfile fix_HRC/merged_GOFA_KFO_final \
      --bmerge fix_HRC_Bar/merged_BCN_GCAT_final.bed fix_HRC_Bar/merged_BCN_GCAT_final.bim fix_HRC_Bar/merged_BCN_GCAT_final.fam \
      --make-bed --out ALL
plink2 --bfile ALL --make-pgen --sort-vars --out ALL_sort
plink2 --pfile ALL_sort --make-bed --out ALL_sort2
plink2 --bfile ALL_sort2 --missing --out ALL_sort2.missing
plink2 --bfile ALL_sort2 --geno 0.1 --make-bed --out merged_ALL  # IMPORTANT: carefully select --geno by reviewing the missingness file
