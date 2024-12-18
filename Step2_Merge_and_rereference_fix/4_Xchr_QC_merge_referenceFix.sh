#!/bin/bash

mkdir -p QC

#======= GOFA n4230 =======#
# These dataset was QC'd
zcat GOFA_n4230_dose.vcf.gz | grep -E '^#|TYPED' > GOFA_n4230_genotyped.vcf

plink2 --vcf GOFA_n4230_genotyped.vcf --set-all-var-ids @:# --chr X  --snps-only --make-bed --out GOFA_n4230

# Remove individual that failed QC (from autosomal chromosomes QC)
plink2 --bfile GOFA_n4230 --remove fail_IDs/GOFA_n4230.fail.IDs --make-bed --out QC/GOFA_n4230

#======= GOFA n285 =======#
plink2 --pedmap ILLUMINA_PLINK_FILE \
	   --remove Samples_to_remove \
	   --chr X,XY  \
	   --make-pgen  \
	   --out GOFA_n285_sorted  \
	   --set-all-var-ids @:#  \
	   --snps-only --sort-vars

plink2 --pfile GOFA_n285_sorted --make-bed --rm-dup exclude-mismatch list --out GOFA_n285

# Remove individual that failed QC (from autosomal chromosomes QC)
# Filter SNP missingness, MAF and HWE (1e-5 for cases)
plink2 --bfile GOFA_n285 --remove fail_IDs/GOFA_n285.fail.IDs --allow-no-sex --geno 0.01 --maf 0.01 --hwe 1e-5 --make-bed --out QC/GOFA_n285

#======= GOFA n95 =======#
plink2 --pedmed ILLUMINA_PLINK_FILE \
	   --chr X,XY  \
	   --make-pgen  \
	   --out GOFA_n95_sorted  
	   --set-all-var-ids @:#  --snps-only --sort-vars

plink2 --pfile GOFA_n95_sorted --make-bed --rm-dup exclude-mismatch list --out GOFA_n95

# Remove individual that failed QC (from autosomal chromosomes QC)
# Filter SNP missingness, MAF and HWE (1e-5 for cases)
plink2 --bfile GOFA_n95 --remove fail_IDs/GOFA_n95.fail.IDs --allow-no-sex --geno 0.01 --maf 0.01 --hwe 1e-5 --make-bed --out QC/GOFA_n95

#======= KFO =======#
plink2 --pedmap ILLUMINA_PLINK_FILE \
	   --chr X,XY  \
	   --make-pgen  \
	   --out KFO_sorted  \
	   --set-all-var-ids @:#  --snps-only --sort-vars

plink2 --pfile KFO_sorted --make-bed --rm-dup exclude-mismatch list --out KFO

# Update phenotype in fam file
awk 'NR == FNR {
    split($1, id_parts, "_");
    id = id_parts[2];
    pheno[id] = $3;   # Phenotype file: '0' for controls, '1' for cases
    next
    }{
    if ($2 in pheno) 
        $6 = pheno[$2] + 1;  # PLINK: '1' for controls, '2' for cases
    print
}' pheno.txt KFO.fam > KFO_1.fam
tr -cd '[:print:]\n' < KFO_1.fam > KFO_1_clean.fam
rm KFO_1.fam
mv KFO.fam KFO.fam_temp
mv KFO_1_clean.fam KFO.fam

# Remove individual that failed QC (from autosomal chromosomes QC)
# Filter SNP missingness, MAF and HWE (1e-5, PLINK2 applied HWE filter on everyone)
plink2 --bfile KFO --remove fail_IDs/KFO.fail.IDs --allow-no-sex --geno 0.01 --maf 0.01 --hwe 1e-5 'midp' --make-bed --out QC/KFO_temp
# Filter HWE additionally in controls (1e-4, more stringent)
plink2 --bfile QC/KFO_temp --hwe 1e-4 'midp' --keep-if "PHENO1==control" --make-bed --out QC/KFO_ctrl
plink2 --bfile QC/KFO_temp --extract QC/KFO_ctrl.bim --make-bed --out QC/KFO

#======= BCN =======#
plink2 --pedmap ILLUMINA_PLINK_FILE  \
	   --chr X,XY  \
	   --make-pgen  \
	   --out BCN_sorted  \
	   --set-all-var-ids @:#  --snps-only --sort-vars

plink2 --pfile BCN_sorted --make-bed --rm-dup exclude-mismatch list --out BCN

# Remove individual that failed QC (from autosomal chromosomes QC)
# Filter SNP missingness, MAF and HWE (1e-5 for cases)
plink2 --bfile BCN --remove fail_IDs/BCN.fail.IDs --allow-no-sex --geno 0.01 --maf 0.01 --hwe 1e-5 --make-bed --out QC/BCN

#======= GCAT =======#
plink2 --bfile GCATCoreSpain_V1  \
	   --chr X  \
	   --make-bed  \
	   --out GCAT  \
	   --rm-dup exclude-mismatch list  \
	   --set-all-var-ids @:#  --snps-only

# Sample QC was done in original paper
# Filter SNP missingness, MAF and HWE (1e-4 for controls)
plink2 --bfile GCAT --allow-no-sex --geno 0.01 --maf 0.01 --hwe 1e-4 --make-bed --out QC/GCAT

#================ Fix to reference ==================#

mkdir -p QC/fix_ref

datasets=("GOFA_n4230" "GOFA_n285" "GOFA_n95" "BCN" "KFO" "GCAT")

for dataset in "${datasets[@]}"; do

	# Set up
	qcdir="QC/fix_ref/$dataset"
	mkdir $qcdir

	# Calculate frequency
	plink2 --bfile QC/GOFA_n4230 --freq --out QC/${dataset}_freq

	awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
	     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' QC/${dataset}_freq.afreq > QC/${dataset}_freq.frq

	sed -i 's/XY/X/g' QC/${dataset}.bim
	sed -i 's/X/23/g' QC/${dataset}.bim
	sed -i 's/XY/X/g' QC/${dataset}_freq.frq
	sed -i 's/X/23/g' QC/${dataset}_freq.frq

	#=== HRC ===#
	perl HRC-1000G-check-bim.pl -b QC/${dataset}.bim -f QC/${dataset}_freq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o $qcdir/

	# Don't generate vcf and chr23 file
	sed -i '/vcf/d' $qcdir/Run-plink.sh
	sed -i '/--chr/d' $qcdir/Run-plink.sh

	# Change the directory of TEMP files
	sed -i 's/rm TEMP/rm QC\/TEMP/g' $qcdir/Run-plink.sh

	# Change plink to plink2
	sed -i '/--a2-allele/ s/plink/plink2/g' $qcdir/Run-plink.sh
	sed -i '/--real-ref-alleles/ s/plink/plink2/g' $qcdir/Run-plink.sh

	# # Run the bash code
	sed -i '1i#!/bin/bash' $qcdir/Run-plink.sh
	source $qcdir/Run-plink.sh

	sed -i 's/X/23/g' $qcdir/${dataset}-updated.bim

done

#=============== Merge ALL ================#

mkdir -p QC/merged

sort QC/fix_ref/GOFA_n4230/GOFA_n4230-updated.bim QC/fix_ref/GOFA_n285/GOFA_n285-updated.bim QC/fix_ref/GOFA_n95/GOFA_n95-updated.bim QC/fix_ref/KFO/KFO-updated.bim QC/fix_ref/BCN/BCN-updated.bim QC/fix_ref/GCAT/GCAT-updated.bim | cut -f 2 | uniq -c | awk '$1 == 6 {print $2}' > QC/merged/overlap_SNPs.txt

cat > QC/files.txt << EOF
QC/fix_ref/GOFA_n285/GOFA_n285-updated.bed QC/fix_ref/GOFA_n285/GOFA_n285-updated.bim QC/fix_ref/GOFA_n285/GOFA_n285-updated.fam
QC/fix_ref/GOFA_n95/GOFA_n95-updated.bed QC/fix_ref/GOFA_n95/GOFA_n95-updated.bim QC/fix_ref/GOFA_n95/GOFA_n95-updated.fam
QC/fix_ref/KFO/KFO-updated.bed QC/fix_ref/KFO/KFO-updated.bim QC/fix_ref/KFO/KFO-updated.fam
QC/fix_ref/BCN/BCN-updated.bed QC/fix_ref/BCN/BCN-updated.bim QC/fix_ref/BCN/BCN-updated.fam
QC/fix_ref/GCAT/GCAT-updated.bed QC/fix_ref/GCAT/GCAT-updated.bim QC/fix_ref/GCAT/GCAT-updated.fam
EOF

plink --bfile QC/fix_ref/GOFA_n4230/GOFA_n4230-updated --merge-list QC/files.txt --extract QC/merged/overlap_SNPs.txt --make-bed --out QC/merged/ALL_overlap

# Calculate frequency
plink2 --bfile QC/merged/ALL_overlap --freq --out QC/merged/ALL_overlap_freq

awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' QC/merged/ALL_overlap_freq.afreq > QC/merged/ALL_overlap_freq.frq

sed -i 's/X/23/g' QC/merged/ALL_overlap.bim
sed -i 's/X/23/g' QC/merged/ALL_overlap_freq.frq

#=== HRC ===#
perl HRC-1000G-check-bim.pl -b QC/merged/ALL_overlap.bim -f QC/merged/ALL_overlap_freq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o QC/merged/

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' QC/merged/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' QC/merged/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' QC/merged/Run-plink.sh
source QC/merged/Run-plink.sh

# Generate VCF file
bgzip -c QC/merged/ALL_overlap-updated-chr23.vcf > QC/merged/ALL_overlap-updated-chr23.vcf.gz
tabix --csi -p vcf QC/merged/ALL_overlap-updated-chr23.vcf.gz

# Fix ploidy
bcftools +guess-ploidy -g b37 QC/merged/ALL_overlap-updated-chr23.vcf.gz > QC/merged/samples.txt
bcftools +fixploidy QC/merged/ALL_overlap-updated-chr23.vcf -Ov -o QC/merged/ALL_overlap_final.vcf -- -s QC/merged/samples.txt
bgzip -c QC/merged/ALL_overlap_final.vcf > QC/merged/ALL_overlap_final.vcf.gz


#=============== Merge German ================#

mkdir -p QC/merged_german

sort QC/fix_ref/GOFA_n4230/GOFA_n4230-updated.bim QC/fix_ref/GOFA_n285/GOFA_n285-updated.bim QC/fix_ref/GOFA_n95/GOFA_n95-updated.bim QC/fix_ref/KFO/KFO-updated.bim | cut -f 2 | uniq -c | awk '$1 == 4 {print $2}' > QC/merged_german/overlap_SNPs.txt

cat > QC/files_german.txt << EOF
QC/fix_ref/GOFA_n285/GOFA_n285-updated.bed QC/fix_ref/GOFA_n285/GOFA_n285-updated.bim QC/fix_ref/GOFA_n285/GOFA_n285-updated.fam
QC/fix_ref/GOFA_n95/GOFA_n95-updated.bed QC/fix_ref/GOFA_n95/GOFA_n95-updated.bim QC/fix_ref/GOFA_n95/GOFA_n95-updated.fam
QC/fix_ref/KFO/KFO-updated.bed QC/fix_ref/KFO/KFO-updated.bim QC/fix_ref/KFO/KFO-updated.fam
EOF

plink --bfile QC/fix_ref/GOFA_n4230/GOFA_n4230-updated --merge-list QC/files_german.txt --extract QC/merged_german/overlap_SNPs.txt --make-bed --out QC/merged_german/GERMAN_overlap

# Calculate frequency
plink2 --bfile QC/merged_german/GERMAN_overlap --freq --out QC/merged_german/GERMAN_overlap_freq

awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' QC/merged_german/GERMAN_overlap_freq.afreq > QC/merged_german/GERMAN_overlap_freq.frq

sed -i 's/X/23/g' QC/merged_german/GERMAN_overlap.bim
sed -i 's/X/23/g' QC/merged_german/GERMAN_overlap_freq.frq

#=== HRC ===#
perl HRC-1000G-check-bim.pl -b QC/merged_german/GERMAN_overlap.bim -f QC/merged_german/GERMAN_overlap_freq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o QC/merged_german/

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' QC/merged_german/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' QC/merged_german/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' QC/merged_german/Run-plink.sh
source QC/merged_german/Run-plink.sh

# Generate VCF file
bgzip -c QC/merged_german/GERMAN_overlap-updated-chr23.vcf > QC/merged_german/GERMAN_overlap-updated-chr23.vcf.gz
tabix --csi -p vcf QC/merged_german/GERMAN_overlap-updated-chr23.vcf.gz

# Fix ploidy
bcftools +guess-ploidy -g b37 QC/merged_german/GERMAN_overlap-updated-chr23.vcf.gz > QC/merged_german/samples.txt
bcftools +fixploidy QC/merged_german/GERMAN_overlap-updated-chr23.vcf -Ov -o QC/merged_german/GERMAN_overlap_final.vcf -- -s QC/merged_german/samples.txt
bgzip -c QC/merged_german/GERMAN_overlap_final.vcf > QC/merged_german/GERMAN_overlap_final.vcf.gz


#=============== Merge Spain ================#

mkdir -p QC/merged_spain

sort QC/fix_ref/BCN/BCN-updated.bim QC/fix_ref/GCAT/GCAT-updated.bim | cut -f 2 | uniq -c | awk '$1 == 2 {print $2}' > QC/merged_spain/overlap_SNPs.txt

cat > QC/files_spain.txt << EOF
QC/fix_ref/GCAT/GCAT-updated.bed QC/fix_ref/GCAT/GCAT-updated.bim QC/fix_ref/GCAT/GCAT-updated.fam
EOF

plink --bfile QC/fix_ref/BCN/BCN-updated --merge-list QC/files_spain.txt --extract QC/merged_spain/overlap_SNPs.txt --make-bed --out QC/merged_spain/SPAIN_overlap

# Calculate frequency
plink2 --bfile QC/merged_spain/SPAIN_overlap --freq --out QC/merged_spain/SPAIN_overlap_freq

awk 'BEGIN {OFS="\t"; print "CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"}
     NR > 1 {print $1, $2, $4, $3, sprintf("%.5f", $5), $6}' QC/merged_spain/SPAIN_overlap_freq.afreq > QC/merged_spain/SPAIN_overlap_freq.frq

sed -i 's/X/23/g' QC/merged_spain/SPAIN_overlap.bim
sed -i 's/X/23/g' QC/merged_spain/SPAINN_overlap_freq.frq

#=== HRC ===#
perl HRC-1000G-check-bim.pl -b QC/merged_spain/SPAIN_overlap.bim -f QC/merged_spain/SPAIN_overlap_freq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o QC/merged_spain/

# Change plink to plink2
sed -i '/--a2-allele/ s/plink/plink2/g' QC/merged_spain/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' QC/merged_spain/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' QC/merged_spain/Run-plink.sh
source QC/merged_spain/Run-plink.sh

# Generate VCF file
bgzip -c QC/merged_spain/SPAIN_overlap-updated-chr23.vcf > QC/merged_spain/SPAIN_overlap-updated-chr23.vcf.gz
tabix --csi -p vcf QC/merged_spain/SPAIN_overlap-updated-chr23.vcf.gz

# Fix ploidy
bcftools +guess-ploidy -g b37 QC/merged_spain/SPAIN_overlap-updated-chr23.vcf.gz > QC/merged_spain/samples.txt
bcftools +fixploidy QC/merged_spain/SPAIN_overlap-updated-chr23.vcf -Ov -o QC/merged_spain/SPAIN_overlap_final.vcf -- -s QC/merged_spain/samples.txt
bgzip -c QC/merged_spain/SPAIN_overlap_final.vcf > QC/merged_spain/SPAIN_overlap_final.vcf.gz
