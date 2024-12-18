#!/bin/bash

mkdir -p heritability
mkdir -p heritability/temp

gcta="PATH_TO_EXECUTABLE_GCTA"
bfile="PATH_TO_PLINK_FILE"

# Prepare phenotype files
pheno="/fast/AG_Lee/Siaw/Sex_specific_GWAS/datasets/phenotypes/pheno_full.txt"
awk 'BEGIN { OFS = "\t"} NR > 1 {print "0", $1 "_" $2, $4}' $pheno > heritability/pheno.txt
awk 'BEGIN { OFS = "\t"} NR > 1 {print "0", $1 "_" $2, $3}' $pheno > heritability/sex.txt
awk 'BEGIN { OFS = "\t"} NR > 1 && $3 == "Female" {print "0", $1 "_" $2, $4}' $pheno > heritability/pheno_female.txt
awk 'BEGIN { OFS = "\t"} NR > 1 && $3 == "Male" {print "0", $1 "_" $2, $4}' $pheno > heritability/pheno_male.txt

# Prepare covariate files
awk 'BEGIN { OFS = "\t"} NR > 1 {print "0", $1, $2, $3, $4, $5, $6, $7, $8, $9}' /GWAS/GWAS_GOFA_KFO/all/pheno_PCs.txt > heritability/PCs.txt

# GREML in imputed data (https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata)
# Step 1: Segment based LD score
$gcta --bfile $bfile --maf 0.05 --ld-score-region 200 --out heritability/herit

# Step 2: Stratify the SNPs by the LD scores of individuals SNPs in R
Rscript GREML_STEP2_stratify_SNPs_by_LD.R heritability/herit.score.ld heritability/herit

# Step 3: Making GRM using SNPs stratified into different groups
for i in $(seq 1 4); do
    $gcta --bfile $bfile --autosome --extract heritability/herit_snp_group${i}.txt --maf 0.05 --max-maf 0.1 --make-grm-bin --out heritability/herit_group${i}_maf0.1_grm --thread-num 4;
done;
for i in $(seq 1 4); do
    $gcta --bfile $bfile --autosome --extract heritability/herit_snp_group${i}.txt --maf 0.1 --max-maf 0.2 --make-grm-bin --out heritability/herit_group${i}_maf0.2_grm --thread-num 4;
done;
for i in $(seq 1 4); do
    $gcta --bfile $bfile --autosome --extract heritability/herit_snp_group${i}.txt --maf 0.2 --max-maf 0.3 --make-grm-bin --out heritability/herit_group${i}_maf0.3_grm --thread-num 4;
done;
for i in $(seq 1 4); do
    $gcta --bfile $bfile --autosome --extract heritability/herit_snp_group${i}.txt --maf 0.3 --max-maf 0.4 --make-grm-bin --out heritability/herit_group${i}_maf0.4_grm --thread-num 4;
done;
for i in $(seq 1 4); do
    $gcta --bfile $bfile --autosome --extract heritability/herit_snp_group${i}.txt --maf 0.4 --max-maf 0.5 --make-grm-bin --out heritability/herit_group${i}_maf0.5_grm --thread-num 4;
done;

# Step 4: REML analysis with multiple GRMs
maf_values=("0.1" "0.2" "0.3" "0.4" "0.5")
groups=("group1" "group2" "group3" "group4")

# Loop over each group and MAF value to create the necessary entries
for group in "${groups[@]}"; do
  for maf in "${maf_values[@]}"; do
    echo "heritability/herit_${group}_maf${maf}_grm" >> heritability/herit_multi_GRMs.txt
  done
done

$gcta --reml --mgrm heritability/herit_multi_GRMs.txt --pheno heritability/pheno.txt --prevalence 0.027 --reml-est-fix --covar heritability/sex.txt --qcovar heritability/PCs.txt --out heritability/GREML_herit_pre0.027

mv heritability/herit* heritability/temp/













datasets=("GERMAN" "SPAIN" "ALL")

for dataset in "${datasets[@]}"; do

    # Dynamically construct the bfile path
    varname="bfile_${dataset}"
    bfile="${!varname}"

    # Step 1: Segment based LD score
    $gcta --bfile $bfile --maf 0.05 --ld-score-region 200 --out heritability/${dataset}

    # Step 2: Stratify the SNPs by the LD scores of individuals SNPs in R
    Rscript GREML_STEP2_stratify_SNPs_by_LD.R heritability/${dataset}.score.ld heritability/${dataset}

    # check wc heritability/${dataset} (less SNPs than 7652558 [ALL])

    # Step 3: Making GRM using SNPs stratified into different groups
    for i in $(seq 1 4); do
        $gcta --bfile $bfile --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.05 --max-maf 0.1 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.1_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile $bfile --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.1 --max-maf 0.2 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.2_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile $bfile --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.2 --max-maf 0.3 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.3_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile $bfile --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.3 --max-maf 0.4 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.4_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile $bfile --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.4 --max-maf 0.5 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.5_grm --thread-num 4;
    done;

    # Step 4: REML analysis with multiple GRMs
    maf_values=("0.1" "0.2" "0.3" "0.4" "0.5")
    groups=("group1" "group2" "group3" "group4")

    # Loop over each group and MAF value to create the necessary entries
    for group in "${groups[@]}"; do
      for maf in "${maf_values[@]}"; do
        echo "heritability/${dataset}_${group}_maf${maf}_grm" >> heritability/${dataset}_multi_GRMs.txt
      done
    done

    $gcta --reml --mgrm heritability/${dataset}_multi_GRMs.txt --pheno heritability/pheno.txt --prevalence 0.027 --reml-est-fix --covar heritability/sex.txt --qcovar heritability/${dataset}_PCs.txt --out heritability/GREML_${dataset}_pre0.027

    mv heritability/${dataset}* heritability/temp/

done

# GERMAN sex-specific heritability
plink2 --bfile ../imputed_data_Sanger/GOFA_KFO.INFO0.5_0.8.MAF0.001.HWE1e-12/GOFA_KFO.MAF0.001.INFO0.5_0.8.HWE1e-12 --keep heritability/pheno_female.txt --make-bed --out heritability/GERMAN_female
plink2 --bfile ../imputed_data_Sanger/GOFA_KFO.INFO0.5_0.8.MAF0.001.HWE1e-12/GOFA_KFO.MAF0.001.INFO0.5_0.8.HWE1e-12 --keep heritability/pheno_male.txt --make-bed --out heritability/GERMAN_male

# Prepare covariate files
awk 'BEGIN { OFS = "\t"} NR > 1 {print "0", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ../GWAS/GWAS_GOFA_KFO/female/pheno_PCs.txt > heritability/GERMAN_female_PCs.txt
awk 'BEGIN { OFS = "\t"} NR > 1 {print "0", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ../GWAS/GWAS_GOFA_KFO/male/pheno_PCs.txt > heritability/GERMAN_male_PCs.txt

datasets=("GERMAN_female" "GERMAN_male")

for dataset in "${datasets[@]}"; do

    # Step 1: Segment based LD score
    $gcta --bfile heritability/${dataset} --maf 0.05 --ld-score-region 200 --out heritability/${dataset}

    # Step 2: Stratify the SNPs by the LD scores of individuals SNPs in R
    Rscript GREML_STEP2_stratify_SNPs_by_LD.R heritability/${dataset}.score.ld heritability/${dataset}

    # Step 3: Making GRM using SNPs stratified into different groups
    for i in $(seq 1 4); do
        $gcta --bfile heritability/${dataset} --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.05 --max-maf 0.1 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.1_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile heritability/${dataset} --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.1 --max-maf 0.2 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.2_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile heritability/${dataset} --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.2 --max-maf 0.3 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.3_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile heritability/${dataset} --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.3 --max-maf 0.4 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.4_grm --thread-num 4;
    done;
    for i in $(seq 1 4); do
        $gcta --bfile heritability/${dataset} --autosome --extract heritability/${dataset}_snp_group${i}.txt --maf 0.4 --max-maf 0.5 --make-grm-bin --out heritability/${dataset}_group${i}_maf0.5_grm --thread-num 4;
    done;

    # Step 4: REML analysis with multiple GRMs
    maf_values=("0.1" "0.2" "0.3" "0.4" "0.5")
    groups=("group1" "group2" "group3" "group4")

    # Loop over each group and MAF value to create the necessary entries
    for group in "${groups[@]}"; do
      for maf in "${maf_values[@]}"; do
        echo "heritability/${dataset}_${group}_maf${maf}_grm" >> heritability/${dataset}_multi_GRMs.txt
      done
    done

    if [ "$dataset" == "GERMAN_female" ]; then
        $gcta --reml --mgrm heritability/${dataset}_multi_GRMs.txt \
              --pheno heritability/pheno_female.txt \
              --prevalence 0.027 \
              --reml-est-fix \
              --qcovar heritability/${dataset}_PCs.txt \
              --out heritability/GREML_${dataset}_pre0.027
    else
        $gcta --reml --mgrm heritability/${dataset}_multi_GRMs.txt \
              --pheno heritability/pheno_male.txt \
              --prevalence 0.027 \
              --reml-est-fix \
              --qcovar heritability/${dataset}_PCs.txt \
              --out heritability/GREML_${dataset}_pre0.027
    fi

    mv heritability/${dataset}* heritability/temp/

done