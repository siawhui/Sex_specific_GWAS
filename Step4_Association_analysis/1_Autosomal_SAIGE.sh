#!/bin/bash

# Install SAIGE from guix

cohort="GERMAN"  # SPAIN, ALL

############## GWAS ALL ########################

# Set up
dir="PATH_TO_OUTPUT_DIRECTORY"
vcfFile="PATH_TO_VCF_FILE"
plinkFile="PATH_TO_PLINK_FILE"
phenoFile="phenotype.txt"
phenoCol="FA"
qCovarColList="sex"
outputPrefix="$dir/${cohort}_all"
sampleIDColinphenoFile="IID"
nThreads="32"

# Create directory
mkdir -p "$dir"

# Get IDs with pheno
awk '{print 0"\t"$1}' $phenoFile > $dir/IDs.txt

# Run PCA
plink2 --bfile $plinkFile --keep $dir/IDs.txt --rm-dup exclude-mismatch list --make-bed --out $dir/data_rmdup
plink2 --bfile $dir/data_rmdup --indep-pairwise 50 5 0.2 --out $dir/data_prune
plink2 --bfile $dir/data_rmdup --extract $dir/data_prune.prune.in --make-bed --out $dir/data_pruned
plink2 --bfile $dir/data_pruned --pca 50 --make-bed --out $dir/GWAS_PCA
plink2 --bfile $plinkFile --keep $dir/IDs.txt --rm-dup exclude-mismatch list --make-bed --out $dir/GWAS_input

# Get elbow point to determine number of PCs for covariates
# Calculate the second derivative (n=2) of the eigenvalues and get the maximum (argmax)
# Correct the index to refer back to the original array of eigenvalues (+2)
# Since python vector starts with index 0, so +1 to get correct index
data_file_path="$dir/GWAS_PCA.eigenval"
elbow_point=$(python3 script_kneed.py "$data_file_path")

# Plot PCA scree plot
Rscript plot_PCA.R "$dir" "$elbow_point"

# Dynamically create the covariate list
covarColList="sex"
for i in $(seq 1 $elbow_point); do
  covarColList+=",PC$i"
done

# Dynamically generate the covariate column selection
cols=""
for i in $(seq 3 $((2 + elbow_point))); do
  cols="$cols\$$i,"
done

# Remove the trailing comma
cols="${cols%,}"

# Add the first 2 columns (ID columns) followed by the dynamically generated columns
awk_command="awk 'BEGIN{FS=OFS=\"\t\"} {print \$2, $cols}' $dir/GWAS_PCA.eigenvec > $dir/temp_pca.txt"

# Execute the dynamically generated awk command
eval "$awk_command"

# Merge the phenotype file with the dynamically extracted PCs
awk 'NR==FNR {s[$1] = $2; f[$1] = $3; next} {print $0"\t"s[$1]"\t"f[$1]}' $phenoFile $dir/temp_pca.txt > $dir/pheno_PCs.txt

rm $dir/temp_*

# Run Rscript to fit null logistic mixed model (binary pheno)
Rscript step1_fitNULLGLMM.R     \
    --plinkFile=$dir/GWAS_input  \
    --phenoFile=$dir/pheno_PCs.txt \
    --phenoCol=$phenoCol \
    --covarColList=$covarColList \
    --qCovarColList=$qCovarColList  \
    --sampleIDColinphenoFile=$sampleIDColinphenoFile \
    --traitType=binary        \
    --outputPrefix=$outputPrefix \
    --nThreads=$nThreads    \
    --IsOverwriteVarianceRatioFile=TRUE

for chr in {1..22}; do

    Rscript step2_SPAtests.R        \
        --vcfFile=$vcfFile       \
        --vcfFileIndex=${vcfFile}.csi       \
        --vcfField=DS       \
        --SAIGEOutputFile=$outputPrefix.chr$chr.txt \
        --chrom=$chr       \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile=$outputPrefix.rda \
        --varianceRatioFile=$outputPrefix.varianceRatio.txt   \
        --LOCO=TRUE \
        --is_Firth_beta=TRUE    \
        --pCutoffforFirth=0.01 \
        --is_output_moreDetails=TRUE \
        --is_imputed_data=TRUE
done

# Concatenate all chromosomes
awk 'FNR==1 && NR!=1 { next; } { print; }' ${outputPrefix}.chr{1..22}.txt > ${outputPrefix}.temp.txt

# Sort by p-value
(head -n 1 ${outputPrefix}.temp.txt && tail -n +2 ${outputPrefix}.temp.txt | sort -g -k13,13) > ${outputPrefix}.GWAS.txt

# Remove the temporary file
rm ${outputPrefix}.temp.txt

# FUMA format
awk 'BEGIN {OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "SE", "P"; next} {print $1, $2, $4, $5, $9, $10, $13}' ${outputPrefix}.GWAS.txt > ${outputPrefix}.GWAS_FUMA.txt


############## GWAS FEMALE and MALE ########################

sex_vec=("female" "male")

for sex in "${sex_vec[@]}"; do

    # Set up
    dir="PATH_TO_OUTPUT_DIRECTORY"
    vcfFile="PATH_TO_VCF_FILE"
    plinkFile="PATH_TO_PLINK_FILE"
    phenoFile='phenotype.txt'
    phenoCol='FA'
    outputPrefix="$dir/${cohort}_${sex}"
    sampleIDColinphenoFile='IID'
    nThreads="32"

    # Create directory
    mkdir -p "$dir"

    # Get IDs
    if [ "$sex" == "female" ]; then
        sex_str="Female"
    elif [ "$sex" == "male" ]; then
        sex_str="Male"
    fi

    # Select all controls
    awk -v sex_str="$sex_str" '$2 == sex_str || $3 == 0 {print "0\t"$1}' $phenoFile > $dir/keep_IDs.txt

    # Run PCA
    plink2 --bfile $plinkFile --keep $dir/keep_IDs.txt --rm-dup exclude-mismatch list --make-bed --out $dir/data_rmdup
    plink2 --bfile $dir/data_rmdup --indep-pairwise 50 5 0.2 --out $dir/data_prune
    plink2 --bfile $dir/data_rmdup --extract $dir/data_prune.prune.in --make-bed --out $dir/data_pruned
    plink2 --bfile $dir/data_pruned --pca 50 --make-bed --out $dir/GWAS_PCA
    plink2 --bfile $plinkFile --keep $dir/keep_IDs.txt --rm-dup exclude-mismatch list --make-bed --out $dir/GWAS_input

    # Get elbow point to determine number of PCs for covariates
	data_file_path="$dir/GWAS_PCA.eigenval"
    elbow_point=$(python3 script_kneed.py "$data_file_path")

	# Plot PCA scree plot
	Rscript plot_PCA.R "$dir" "$elbow_point"

	# Dynamically create the covariate list
	covarColList=""
	for i in $(seq 1 $elbow_point); do
	  covarColList+="PC$i,"
	done
	covarColList="${covarColList%,}"

	# Dynamically generate the covariate column selection
	cols=""
	for i in $(seq 3 $((2 + elbow_point))); do
	  cols="$cols\$$i,"
	done
	cols="${cols%,}"

	# Add the first 2 columns (ID columns) followed by the dynamically generated columns
	awk_command="awk 'BEGIN{FS=OFS=\"\t\"} {print \$2, $cols}' $dir/GWAS_PCA.eigenvec > $dir/temp_pca.txt"

	# Execute the dynamically generated awk command
	eval "$awk_command"

	# Merge the phenotype file with the dynamically extracted PCs
	awk 'NR==FNR {s[$1] = $2; f[$1] = $3; next} {print $0"\t"s[$1]"\t"f[$1]}' $phenoFile $dir/temp_pca.txt > $dir/pheno_PCs.txt

    rm $dir/temp_*

    # Run Rscript to fit null logistic mixed model (binary pheno)
    Rscript step1_fitNULLGLMM.R     \
        --plinkFile=$dir/GWAS_input  \
        --phenoFile=$dir/pheno_PCs.txt \
        --phenoCol=$phenoCol \
        --covarColList=$covarColList \
        --sampleIDColinphenoFile=$sampleIDColinphenoFile \
        --traitType=binary        \
        --outputPrefix=$outputPrefix \
        --nThreads=$nThreads    \
        --IsOverwriteVarianceRatioFile=TRUE

    for chr in {1..22}; do
        Rscript step2_SPAtests.R        \
            --vcfFile=$vcfFile       \
            --vcfFileIndex=${vcfFile}.csi       \
            --vcfField=DS       \
            --SAIGEOutputFile=$outputPrefix.chr$chr.txt \
            --chrom=$chr       \
            --minMAF=0 \
            --minMAC=20 \
            --GMMATmodelFile=$outputPrefix.rda \
            --varianceRatioFile=$outputPrefix.varianceRatio.txt   \
            --LOCO=TRUE \
            --is_Firth_beta=TRUE    \
            --pCutoffforFirth=0.01 \
            --is_output_moreDetails=TRUE \
            --is_imputed_data=TRUE
    done

    # Concatenate all chromosomes
    awk 'FNR==1 && NR!=1 { next; } { print; }' ${outputPrefix}.chr{1..22}.txt > ${outputPrefix}.temp.txt

    # Sort by p-value
    (head -n 1 ${outputPrefix}.temp.txt && tail -n +2 ${outputPrefix}.temp.txt | sort -g -k13,13) > ${outputPrefix}.GWAS.txt

    # Remove the temporary file
    rm ${outputPrefix}.temp.txt

    # FUMA format
    awk 'BEGIN {OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "SE", "P"; next} {print $1, $2, $4, $5, $9, $10, $13}' ${outputPrefix}.GWAS.txt > ${outputPrefix}.GWAS_FUMA.txt

done
