#!/bin/bash

cohort="GERMAN"  # SPAIN, ALL

datadir="$cohort.INFO0.5_0.8.MAF0.001.HWE1e-12"

######### Prepare data for XWAS ##########

# Get genotype column only
bcftools query -f '%ID\t[%GT\t]\n' $datadir/chrX.MAF0.001.INFO0.5_0.8.HWE1e-12.recode.vcf > chrX_dosage.txt

# Transform the genotype to (0/1/2) for female and (0/1) for male
awk 'NR>=1 {for(i=2; i<=NF; i++) {if ($i == "0|0") $i = "0"; else if ($i == "0|1" || $i == "1|0") $i = "1"; else if ($i == "1|1") $i = "2"}; print $0}' OFS="\t" chrX_dosage.txt > chrX_dosage_transformed.txt

# Remove Trailing Delimiters
sed 's/[ \t]*$//' chrX_dosage_transformed.txt > chrX_dosage_transformed_fixed.txt

# Transpose the table
datamash transpose < chrX_dosage_transformed_fixed.txt > chrX_dosage_transposed.txt

# Add sample IDs
bcftools query -l $datadir/chrX.MAF0.001.INFO0.5_0.8.HWE1e-12.recode.vcf > sample_ids.txt

echo "ID" > sample_ids_final.txt  # Create a new output file with the header
cat sample_ids.txt >> sample_ids_final.txt  # Append sample IDs to the output file

paste sample_ids_final.txt chrX_dosage_transposed.txt > chrX_dosage_final.txt

# Modify fam file
awk 'BEGIN { OFS = "\t"; print "FID", "IID", "SEX", "PHENOTYPE" } 
NR > 1 {
    if ($5 == "original_n4230") {
        fid = "0"
        iid = "0_" $7 "_" $7  # ID format for GOFA_n4230
    } else {
        fid = "0"
        iid = $1 "_" $2  # ID format for other cohorts
    }
    sex = ($3 == "Male") ? 1 : 2  # male=1, female=2
    pheno = ($4 == 0) ? 1 : 2  # control=1, cases=2
    print fid, iid, sex, pheno
}' phenotype.txt > pheno_update.txt

# Update sex info
plink2 --bfile chrX.MAF0.001.INFO0.5_0.8.HWE1e-12 --update-sex pheno_update.txt col-num=3 --make-bed --out updated_sex

# Update phenotype info
plink2 --bfile updated_sex --pheno pheno_update.txt --pheno-col-nums 4 --make-bed --out chrX

################## Include PCs as covariates ########################

# Use PCs from Autosomal GWAS
pca_file="pheno_PCs.txt"

# IMPORTANT: change the number of PCs in header and the number of columns printed
awk 'BEGIN { OFS = "\t"; print "ID", "gender", "FA", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6" } 
NR == FNR {  # Process the first file (pheno file)
    iid = $1 "_" $2
    if ($5 == "original_n4230") {
        id[iid] = "0_" $7 "_" $7  # ID format for cohort1
    } else {
        id[iid] = $1 "_" $2       # ID format for other samples
    }
    sex[iid] = $3
    fa[iid] = $4
    next
} 
{  # Process the second file (pca_file)
    iid = $1
    if (iid in id) {
        print id[iid], sex[iid], fa[iid], $2, $3, $4, $5, $6, $7
    }
}' $pheno $pca_file > pheno_PCs.txt

# Run analysis
Rscript run_XCMAX4.R chrX_dosage_final.txt pheno_PCs.txt results.txt

# FUMA format
cut -f 1-6 results.txt > results_FUMA.txt