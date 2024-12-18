#!/bin/bash

# Prepare phenotype file
awk 'BEGIN { OFS = "\t"; print "FID", "IID", "SEX", "PHENOTYPE" } 
NR > 1 {
    if ($5 == "original_n4230") {
        fid = "0"
        iid = "0_" $7 "_" $7  # ID format for cohort1
    } else {
        fid = "0"
        iid = $1 "_" $2  # ID format for other samples
    }
    sex = ($3 == "Male") ? 1 : 2
    pheno = ($4 == 0) ? 1 : 2
    print fid, iid, sex, pheno
}' /fast/AG_Lee/Siaw/Sex_specific_GWAS/datasets/phenotypes/pheno_full.txt > pheno_update.txt

# Prepare PCs
pheno="phenotype.txt"

# ALL
pca_file="GWAS/GWAS_ALL/all/pheno_PCs.txt"
awk 'BEGIN { OFS = "\t"} 
NR == FNR {  # Process the first file (pheno file)
    id = $1 "_" $2
    if ($5 == "original_n4230") {
        iid[id] = "0_" $7 "_" $7  # IID format for GOFA_n4230
        fid[id] = "0"
    } else {
        iid[id] = $1 "_" $2  # IID format for other cohorts
        fid[id] = "0"
    }
    next
} 
{  # Process the second file (pca_file)
    id2 = $1
    if (id2 in iid) {
        print fid[id2], iid[id2], $2, $3, $4, $5, $6, $7
    }
}' $pheno $pca_file > ALL_pheno_PCs.txt

# GERMAN
pca_file="GWAS/GWAS_GOFA_KFO/all/pheno_PCs.txt"
awk 'BEGIN { OFS = "\t"} 
NR == FNR {  # Process the first file (pheno file)
    id = $1 "_" $2
    if ($5 == "original_n4230") {
        iid[id] = "0_" $7 "_" $7  # IID format for original GOFA
        fid[id] = "0"
    } else {
        iid[id] = $1 "_" $2 
        fid[id] = "0"
    }
    next
} 
{  # Process the second file (pca_file)
    id2 = $1
    if (id2 in iid) {
        print fid[id2], iid[id2], $2, $3, $4, $5, $6, $7, $8, $9
    }
}' $pheno $pca_file > GERMAN_pheno_PCs.txt

# SPAIN
pca_file="GWAS/GWAS_BCN_GCAT/all/pheno_PCs.txt"
awk 'BEGIN { OFS = "\t"} 
NR == FNR {  # Process the first file (pheno file)
    id = $1 "_" $2
    iid[id] = $1 "_" $2 
    fid[id] = "0"
    next
} 
{  # Process the second file (pca_file)
    id2 = $1
    if (id2 in iid) {
        print fid[id2], iid[id2], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
    }
}' $pheno $pca_file > SPAIN_pheno_PCs.txt

#========== ALL ==========#
bfile="PATH_TO_MERGED_PLINK_FILE"
plink2 --bfile $bfile --update-sex pheno_update.txt col-num=3 --make-bed --out ALL_updated_sex
plink2 --bfile ALL_updated_sex --pheno pheno_update.txt --pheno-col-nums 4 --make-bed --export ped --out ALL_chrX
plink --file ALL_chrX --make-bed --out ALL_chrX_plink1
../bin/xwas --bfile ALL_chrX_plink1 --xwas --strat-sex --fishers --stouffers --xbeta --multi-xchr-model --covar ALL_pheno_PCs.txt --out xwas_ALL

#========== GERMAN ==========#
# remove PAR region from german cohort
bfile="PATH_TO_GERMAN_PLINK_FILE"
plink2 --bfile $bfile --update-sex pheno_update.txt col-num=3 --split-par b37 --make-bed --out GERMAN_updated_sex
plink2 --bfile GERMAN_updated_sex --pheno pheno_update.txt --pheno-col-nums 4 --chr X --make-bed --export ped --out GERMAN_chrX
plink --file GERMAN_chrX --make-bed --out GERMAN_chrX_plink1
../bin/xwas --bfile GERMAN_chrX_plink1 --xwas --strat-sex --fishers --stouffers --xbeta --multi-xchr-model --covar GERMAN_pheno_PCs.txt --out xwas_GERMAN

#========== SPAIN ==========#
bfile="PATH_TO_SPAIN_PLINK_FILE"
plink2 --bfile $$bfile --update-sex pheno_update.txt col-num=3 --make-bed --out SPAIN_updated_sex
plink2 --bfile SPAIN_updated_sex --pheno pheno_update.txt --pheno-col-nums 4 --make-bed --export ped --out SPAIN_chrX
plink --file SPAIN_chrX --make-bed --out SPAIN_chrX_plink1
../bin/xwas --bfile SPAIN_chrX_plink1 --xwas --strat-sex --fishers --stouffers --xbeta --multi-xchr-model --covar SPAIN_pheno_PCs.txt --out xwas_SPAIN

#========= Prepare FUMA tables ==========#
# Extract the header (first line) and process the rest of the data
#=== ALL ===#
awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $7, $8}' xwas_ALL.02.xstrat.logistic > xwas_ALL.02.xstrat_male.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $8 != "NA" {print $1, $3, $4, $5, $7, $8}' xwas_ALL.02.xstrat.logistic | sort -k6 -g >> xwas_ALL.02.xstrat_male.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $9, $10}' xwas_ALL.02.xstrat.logistic > xwas_ALL.02.xstrat_female.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $10 != "NA" {print $1, $3, $4, $5, $9, $10}' xwas_ALL.02.xstrat.logistic | sort -k6 -g >> xwas_ALL.02.xstrat_female.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $13, $14}' xwas_ALL.02.xstrat.logistic > xwas_ALL.02.xstrat_all_Stouffer.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $14 != "NA" {print $1, $3, $4, $5, $13, $14}' xwas_ALL.02.xstrat.logistic | sort -k6 -g >> xwas_ALL.02.xstrat_all_Stouffer.logistic

#=== GERMAN ===#
awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $7, $8}' xwas_GERMAN.02.xstrat.logistic > xwas_GERMAN.02.xstrat_male.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $8 != "NA" {print $1, $3, $4, $5, $7, $8}' xwas_GERMAN.02.xstrat.logistic | sort -k6 -g >> xwas_GERMAN.02.xstrat_male.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $9, $10}' xwas_GERMAN.02.xstrat.logistic > xwas_GERMAN.02.xstrat_female.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $10 != "NA" {print $1, $3, $4, $5, $9, $10}' xwas_GERMAN.02.xstrat.logistic | sort -k6 -g >> xwas_GERMAN.02.xstrat_female.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $13, $14}' xwas_GERMAN.02.xstrat.logistic > xwas_GERMAN.02.xstrat_all_Stouffer.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $14 != "NA" {print $1, $3, $4, $5, $13, $14}' xwas_GERMAN.02.xstrat.logistic | sort -k6 -g >> xwas_GERMAN.02.xstrat_all_Stouffer.logistic

#=== SPAIN ===#
awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $7, $8}' xwas_SPAIN.02.xstrat.logistic > xwas_SPAIN.02.xstrat_male.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $8 != "NA" {print $1, $3, $4, $5, $7, $8}' xwas_SPAIN.02.xstrat.logistic | sort -k6 -g >> xwas_SPAIN.02.xstrat_male.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $9, $10}' xwas_SPAIN.02.xstrat.logistic > xwas_SPAIN.02.xstrat_female.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $10 != "NA" {print $1, $3, $4, $5, $9, $10}' xwas_SPAIN.02.xstrat.logistic | sort -k6 -g >> xwas_SPAIN.02.xstrat_female.logistic

awk 'BEGIN {OFS="\t"} NR==1 {print $1, $3, $4, $5, $13, $14}' xwas_SPAIN.02.xstrat.logistic > xwas_SPAIN.02.xstrat_all_Stouffer.logistic
awk 'BEGIN {OFS="\t"} NR == 1 {next} $14 != "NA" {print $1, $3, $4, $5, $13, $14}' xwas_SPAIN.02.xstrat.logistic | sort -k6 -g >> xwas_SPAIN.02.xstrat_all_Stouffer.logistic
