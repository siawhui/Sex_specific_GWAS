#!/bin/bash

# Set up
qcdir='PATH_TO_QC_DIRECTORY'
refdir='PATH_TO_DIRECTORY_CONTAINING_REFERENCE_FILE' 
refname='NAME_OF_REFERENCE_FILE' # 1000G phase3 (PLINK format)
initFile='PATH_TO_GENOTYPE_FILE'
name='PREFIX_OF_OUTPUT_FILE'

# Create directory
mkdir -p "$qcdir"

# Extract missing allele
awk '{if ($5 == "." || $6 == ".") print $2}' $initFile.bim > $name.missing_allele

plink2 --bfile $initFile \
       --chr 1-22  \
       --snps-only  \
       --rm-dup exclude-mismatch list  \
       --make-bed   \
       --out $name  \
       --exclude n95.missing_allele  \
       --set-all-var-ids @:#

#=== Ancestry estimation ===#

#Filter reference and study data for non A-T or G-C SNPs
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink2 --bfile  $refdir/$refname \
       --exclude $qcdir/$refname.ac_gt_snps \
       --make-bed \
       --out $qcdir/$refname.no_ac_gt_snps  \
       --set-all-var-ids @:#

plink2 --bfile  $name \
       --exclude $qcdir/$name.ac_gt_snps \
       --rm-dup exclude-all list \
       --make-bed \
       --out $qcdir/$name.no_ac_gt_snps

# Prune study data
plink2 --bfile  $qcdir/$name.no_ac_gt_snps \
       --exclude range  $refdir/high-LD-regions-hg19-GRCh37.txt \
       --indep-pairwise 50 5 0.2 \
       --out $qcdir/$name.prune

plink2 --bfile  $qcdir/$name.no_ac_gt_snps \
       --extract $qcdir/$name.prune.prune.in \
       --make-pgen \
       --sort-vars \
       --out $qcdir/$name.pruned

plink2 --pfile  $qcdir/$name.pruned \
       --make-bed \
       --out $qcdir/$name.pruned

# Filter reference data for the same SNP set as in study
plink2 --bfile  $qcdir/$refname.no_ac_gt_snps \
       --extract $qcdir/$name.prune.prune.in \
       --rm-dup exclude-mismatch \
       --make-bed \
       --out $qcdir/$refname.pruned

# Check and correct chromosome mismatch (alternative for plink(v1.9, not available on guix) --update-chr)
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1) { $1 = a[$2] } 1' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.updateChr.bim

cp $qcdir/$refname.pruned.bed $qcdir/$refname.updateChr.bed
cp $qcdir/$refname.pruned.fam $qcdir/$refname.updateChr.fam

# Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print $2, a[$2]}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/${refname}.toUpdatePos

# Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/$refname.toFlip

# Upate positions and flip alleles
plink --bfile $qcdir/$refname.updateChr \
      --update-map $qcdir/$refname.toUpdatePos \
      --flip $qcdir/$refname.toFlip \
      --make-bed \
      --out $qcdir/$refname.flipped

# Remove mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
    $qcdir/$refname.mismatch

plink2 --bfile $qcdir/$refname.flipped \
       --exclude $qcdir/$refname.mismatch \
       --make-pgen \
       --sort-vars \
       --out $qcdir/$refname.toMerge

plink2 --pfile $qcdir/$refname.toMerge \
       --make-bed \
       --out $qcdir/$refname.toMerge

awk '{print $2}' $qcdir/$refname.toMerge.bim > $qcdir/$refname.toMerge.snps.txt

plink2 --bfile $qcdir/$name.pruned \
       --extract $qcdir/$refname.toMerge.snps.txt \
       --make-bed \
       --out $qcdir/$name.toMerge

# Merge study genotypes and reference data
plink --bfile $qcdir/$name.toMerge  \
      --bmerge $qcdir/$refname.toMerge.bed $qcdir/$refname.toMerge.bim \
         $qcdir/$refname.toMerge.fam  \
      --make-bed \
      --out $qcdir/$name.merge.$refname

# PCA on the merged data
plink2 --bfile $qcdir/$name.merge.$refname \
       --pca \
       --out $qcdir/$name.$refname