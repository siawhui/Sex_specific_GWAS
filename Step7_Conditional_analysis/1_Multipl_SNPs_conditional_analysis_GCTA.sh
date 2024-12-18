#!/bin/bash

# GCTA path
gcta="gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"

# Parameters
# p-value < 5e-8
# R2 < 0.2 (FUMA: 0.1)
# Window = 1Mb (1000kb) (FUMA: 250kb)

# reference
# GCTA doesn't suggest 1000G due to small sample size (<4000)
# Use the full dataset as LD reference
ref_path="GWAS/GWAS_ALL/all/GWAS_input"

#====== 6 analysis in ALL ======#
analysis=("all" "female" "male" "male_cases_n566" "female_cases_vs_all" "male_cases_vs_all")
sample_size=(5857 4888 5291 4888 5857 5857)

mkdir -p conditional_results 

# Loop through vector1 using the index
for ((i=0; i<6; i++)); do
    ana=${analysis[i]}
    N=${sample_size[i]}

    GWAS_path="GWAS/GWAS_ALL/$ana"
    
    awk -v N=$N 'BEGIN {OFS="\t"} NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"; next} {print $3, $5, $4, $7, $9, $10, $13, N}' $GWAS_path/ALL_$ana.GWAS.txt > $GWAS_path/ALL_$ana.GWAS_GCTA.txt

    $gcta --bfile $ref_path --cojo-file $GWAS_path/ALL_$ana.GWAS_GCTA.txt --cojo-slct --cojo-p 5e-8 --cojo-collinear 0.2 --cojo-wind 1000 --out conditional_results/ALL_$ana.GWAS_cond.txt

done

#====== 4 analysis in GOFA_KFO ======#
analysis=("all" "female" "male" "male_cases_n480")
sample_size=(4753 3917 4273 3917)

# Loop through vector1 using the index
for ((i=0; i<4; i++)); do
    ana=${analysis[i]}
    N=${sample_size[i]}

    GWAS_path="/fast/AG_Lee/Siaw/Sex_specific_GWAS/GWAS/GWAS_GOFA_KFO/$ana"
    
    awk -v N=$N 'BEGIN {OFS="\t"} NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"; next} {print $3, $5, $4, $7, $9, $10, $13, N}' $GWAS_path/GOFA_KFO_$ana.GWAS.txt > $GWAS_path/GOFA_KFO_$ana.GWAS_GCTA.txt

    $gcta --bfile $ref_path --cojo-file $GWAS_path/GOFA_KFO_$ana.GWAS_GCTA.txt --cojo-slct --cojo-p 5e-8 --cojo-collinear 0.2 --cojo-wind 1000 --out conditional_results/GOFA_KFO_$ana.GWAS_cond.txt

done

#====== 3 analysis in BCN_GCAT ======#
analysis=("all" "female" "male")
sample_size=(1104 971 1018)

# Loop through vector1 using the index
for ((i=0; i<3; i++)); do
    ana=${analysis[i]}
    N=${sample_size[i]}

    GWAS_path="/fast/AG_Lee/Siaw/Sex_specific_GWAS/GWAS/GWAS_BCN_GCAT/$ana"
    
    awk -v N=$N 'BEGIN {OFS="\t"} NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"; next} {print $3, $5, $4, $7, $9, $10, $13, N}' $GWAS_path/BCN_GCAT_$ana.GWAS.txt > $GWAS_path/BCN_GCAT_$ana.GWAS_GCTA.txt

    $gcta --bfile $ref_path --cojo-file $GWAS_path/BCN_GCAT_$ana.GWAS_GCTA.txt --cojo-slct --cojo-p 5e-8 --cojo-collinear 0.2 --cojo-wind 1000 --out conditional_results/BCN_GCAT_$ana.GWAS_cond.txt

done