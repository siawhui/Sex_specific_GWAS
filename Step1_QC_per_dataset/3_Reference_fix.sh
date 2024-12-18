#!/bin/bash

# Set up
qcdir='PATH_TO_QC_DIRECTORY'
datadir='PATH_TO_DIRECTORY_CONTAINING_GENOTYPE_FILE'
name='PREFIX_OF_GENOTYPE_FILE'

perlFile='PATH_TO_PERL_SCRIPT'
# HRC or 1000G Imputation preparation and checking
# Downloaded from https://www.chg.ox.ac.uk/~wrayner/tools/
# HRC-1000G-check-bim.pl

HRC='PATH_TO_HRC_REFERENCE'
# Unzipped tab delimited HRC reference (HRC.r1-1.GRCh37.wgs.mac5.sites.tab)
# Downloaded from the Haplotype Reference Consortium Website
# https://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/

mkdir -p $qcdir

# Calculate frequency
plink --bfile $datadir/$name --freq --out $datadir/$name.freq

#=== HRC ===#
perl $perlFile -b $datadir/$name.bim -f $datadir/$name.freq.frq -r $HRC -h -o $qcdir/

# Change plink to plink2 (PLINK1.9 not available on guix)
sed -i '/--a2-allele/ s/plink/plink2/g' $qcdir/Run-plink.sh
sed -i '/--real-ref-alleles/ s/plink/plink2/g' $qcdir/Run-plink.sh

# Run the bash code
sed -i '1i#!/bin/bash' $qcdir/Run-plink.sh
source $qcdir/Run-plink.sh