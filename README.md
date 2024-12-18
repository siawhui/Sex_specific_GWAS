# Sex-specific Genome-Wide Association Study for Food Allergy

This repository contains the code and scripts developed for my Master's thesis: "Sex-Specific Genome-Wide Association Study for Food Allergy." The project focuses on identifying sex-specific genetic loci associated with food allergy using genome-wide association study (GWAS) method.

## Table of Contents

1. [Quality control on each dataset](#quality-control-on-each-dataset-step1_qc_per_dataset)
2. [Merge datasets and reference fix](#merge-datasets-and-reference-fix-step2_merge_and_reference_fix)
3. [Imputation](#imputation)
4. [Quality control after genotype imputation](#quality-control-after-genotype-imputation-step3_qc_after_imputation)
5. [Association analysis](#association-analysis-step4_association_analysis)
6. [Meta-analysis](#meta-analysis-step5_meta_analysis)
7. [Manhattan plot and QQ plot](#manhattan-plot-and-quantile-quantile-qq-plot-step6_manhattan_and_qq_plots)
8. [Conditional analysis](#conditional-analysis-step7_conditional_analysis)
9. [Interaction analysis](#interaction-analysis-step8_interaction_analysis)
10. [Heritability estimation and chi-squared test](Step9_Heritability_and_chisqr_test/)](#heritability-estimation-and-chi-squared-test-step9_heritability_and_chisqr_test)


## Quality control on each dataset (`Step1_QC_per_dataset/`)
This folder contains 3 scripts that perform quality control (QC) in each genotype dataset separately. Before starting the QC, convert the genotype PED/MAP files into PLINK BED/BIM/FAM files. 

1. Check Ancestry (`1_Check_ancestry.sh`)
The genotypes of the study population is combined with genotypes of a reference dataset consisting of individuals from known ethnicities. Principal component analysis (PCA) on this combined genotype panel can then be used to detect population structure down to the level of the reference dataset.

2. QC (`2_PlinkQC.Rmd`)
*plinkQC* performs per-individual and per-marker genotype QC and generate QC report. The per-individual QC includes: (i) check heterozygosity and missingness: for the identification of individuals with outlying missing genotype and/or heterozygosity rates, (ii) check relatedness: for the identification of related individuals, (iii) check ancestry: identification of individuals of divergent ancestry. The per-marker QC includes: (i) check SNP missingnes: for the identifying markers with excessive missing genotype rates, (ii) check Hardy-Weinberg equilibrium (HWE): for the identifying markers showing a significant deviation from HWE, (iii) check minor allele frequency (MAF): for the removal of markers with low MAF.

Individuals and markers that fail the quality control can subsequently be removed with plinkQC to generate a new, clean dataset. 

3. Reference fix (`3_Reference_fix.sh`)
Before combining genotype data from different cohorts, the reference fix script ensures that all datasets are harmonized to a common reference genome (Haplotype Reference Consortium, HRC). This is done using the Perl script, `HRC-1000G-check-bim.pl`, downloaded from [https://www.chg.ox.ac.uk/~wrayner/tools/](https://www.chg.ox.ac.uk/~wrayner/tools/).

This step addresses potential issues such as flipped strands, chromosome or position mismatches, and inconsistent allele annotations across cohorts. By aligning all datasets to the same reference, this process minimizes errors during merging and downstream analysis, ensuring the integrity of the combined genotype data.


## Merge datasets and reference fix (`Step2_Merge_and_reference_fix/`)
This folder contains 4 scripts that handle key steps in merging and preparing genotype datasets from different cohorts, ensuring consistency and data quality:

1. Merging Datasets (`1_Merge.sh`)
Cleaned genotype data are combined into three sets: German, Spanish, and a combined dataset including all cohorts.

2. Checking for Duplicated Individuals (`2_plinkQC_relatedness.Rmd`)
This R Markdown script uses *plinkQC* to identify and resolve duplicate individuals who may have been recruited multiple times across or within cohorts.

3. Reference Fix and Conversion to VCF (`3_Reference_fix.sh`)
The Perl script `HRC-1000G-check-bim.pl` is used to ensure harmonization of genotype data to reference genome used for imputation (HRC) by fixing strand mismatches, chromosome/position inconsistencies, and allele alignment issues. The final output is converted to VCF format for imputation.

4. X Chromosome QC and Reference Fix (`4_Xchr_QC_merge_referenceFix.sh`)
This script performs all QC steps, merging, and reference fixing separately for X chromosome data (similar to autosomal chromosomes). Additionally, it includes an extra step to fix the ploidy of the X chromosome data using `bcftools +fixploidy`.

## Imputation
Imputation of each dataset (German, Spanish, and combined dataset) were done on Sanger Imputation Service webserver ([https://imputation.sanger.ac.uk/](https://imputation.sanger.ac.uk/)) following the instructions [here](https://imputation.sanger.ac.uk/?instructions=1#prepareyourdata). 

Further details are provided below to assist with future replication:
- Create a new job on Sanger Imputation Service
  - Reference panel: Haplotype Reference Consortium (r1.1)
  - Pipeline: pre-phase with EAGLE2 and impute
- Receive an email invitation from Globus, indicating that the private endpoint has been created on the Sanger servers
- Connect to "MDC-OPEN" WIFI on MacBook
- Open Globus ad click "Web: Transfer Files"
- The left-hand panel is local path (currently connected to the folder "GlobusConnect" on MacBook desktop). Make sure the zipped VCF file to be imputed in located in the folder.
- For the right-hand panel, search for `#vrimpute#iserver` (can find it in "Shared With You") and open the folder named with the endpoint ID written in the email.
- Upload one **single zipped VCF file** into the folder
- After receiving an email from Globus confirming the transfer was successful, verify the completed upload by clicking the link in the email from Sanger. This action will close the endpoint, and the imputation server will begin processing the data.
- You can track the status of the imputation job [here](https://imputation.sanger.ac.uk/?status=1) by entering the Run ID. If the job is aborted due to an error, you will be notified by Sanger. In that case, fix the issue and create a new job.
- Once the imputation is complete, you will receive a confirmation email from Sanger. Download the imputed data (in VCF format) using the same steps as for uploading. Please note that the data will be deleted within seven days.


## Quality control after genotype imputation (`Step3_QC_after_imputation/`)
In this folder, two scripts are provided for performing QC after imputation on genetic data.

1. QC in autosomal chromosomes (`1_Filter_INFO_and_MAF_autosomal_chr.sh`)
This script filters variants based on MAF, INFO score, and HWE (in controls) for autosomal chromosomes. The filtering thresholds ensure high-quality variants for downstream analysis.

2. QC in X chromosome (`2_Filter_INFO_and_MAF_Xchr.sh`)
This script applies the same QC steps (MAF, INFO score, and HWE filtering) specifically to the X chromosome.


## Association analysis (`Step4_Association_analysis/`)

This folder contains scripts for performing association analyses on imputed genetic data.

### Main Scripts:
1. Association analysis on autosomes using SAIGE (`1_Autosomal_SAIGE.sh`)
Performs association analysis on autosomal chromosomes using the SAIGE method, covering all samples as well as female and male-specific subsets.
2. Association analysis on X chromosome using XCMAX4 (`2_Xchr_XCMAX4.sh`)
Conducts association analysis specifically on the X chromosome using the XCMAX4 method for all samples.
3. Sex-specific association analysis on X chromosome using XWAS (`3_Xchr_sex_specific_XWAS.sh`)
Executes female- and male-specific association analyses on the X chromosome using the XWAS method.

### Functional Scripts:
1. Select optimal PCs (`script_kneed.py`)
Uses the kneed Python package to determine the optimal number of principal components (PCs) by applying the elbow method, which are then used as covariates in the association analysis.
2. Plot PCA scree plot (`plot_PCA.R`)
Generates a PCA scree plot in R to visualize the proportion of variance explained by each principal component and highlight the optimal number of PCs selected.
3. SAIGE step 1: fitting the null logistic mixed model (`step1_fitNULLGLMM.R`)
Performs the first step of the SAIGE to fit the null logistic mixed model. Downloaded from [here](https://github.com/saigegit/SAIGE/tree/main/extdata).
3. SAIGE step 2: performing single-variant association tests (`step2_SPAtests.R`)
Performs the second step, single-variant association tests, implemented SAIGE. Downloaded from [here](https://github.com/saigegit/SAIGE/tree/main/extdata).
4. XCMAX4 function (`XCMAX4.R`)
An R function for running the XCMAX4 analysis. Downloaded from [here](https://github.com/YoupengSU/XCMAX4/tree/main).
5. Run XCMAX4 (`run_XCMAX4.R`)
An R script to execute the XCMAX4 association analysis.


## Meta-analysis (`Step5_Meta_analysis/`)

This folder contains scripts for performing meta-analysis on association results from German and Spanish cohorts.

#### Main Scripts:
1. Meta-analysis on autosomal chromosomes (`1_Autosomes_Xchr_METAL.sh`)
Performs meta-analysis of the German and Spanish association results for autosomal chromosomes (all, female, male samples) and the X chromosome (all samples) using the METAL tool.
2. Meta-analysis on X chromosomes (`2_Xchr_sex_specific_METAL.sh`)
Conducts meta-analysis on sex-specific association results in the X chromosome using METAL.

#### METAL Scripts:
1. `METAL_script_all.txt`: METAL script for meta-analysis of autosomal chromosomes (all samples).
2. `METAL_script_female.txt`: METAL script for meta-analysis of autosomal chromosomes (female-specific analysis).
3. `METAL_script_male.txt`: METAL script for meta-analysis of autosomal chromosomes (male-specific analysis).
4. `METAL_script_XCMAX4.txt`: METAL script for meta-analysis of the X chromosome (all samples).
5. `METAL_script_XWAS_female.txt`: METAL script for meta-analysis of the X chromosome (female-specific analysis).
6. `METAL_script_XWAS_male.txt`: METAL script for meta-analysis of the X chromosome (male-specific analysis).

## Manhattan plot and quantile-quantile (QQ) plot (`Step6_Manhattan_and_QQ_plots/`)

This folder contains a Python script for visualizing the results of the GWAS analysis.

1. Plot Manhattan plot and QQ plot (`1_Manhattan_and QQ_gwaslab.py`)
Uses the `gwaslab` Python function to generate Manhattan and QQ plots from the summary statistics of the GWAS results.


## Conditional analysis (`Step7_Conditional_analysis/`)

This folder contains scripts for performing conditional analysis on the GWAS results.

1. Multiple SNPs conditional analysis (`1_Multiple_SNPs_conditional_analysis_GCTA.sh`)  
Performs stepwise conditional analysis using GCTA (`--cojo` option) to identify independent significant signals.

2. Single SNP conditional analysis (`2_Single_SNP_conditional_analysis_SAIGE.sh`) 
Conducts single SNP conditional analysis using SAIGE to obtain the conditioned p-values of other SNPs.

3. Regional association plot (`3_LocusZoom.Rmd`)
Generates a regional association plot using the `locuszoomr` R package to visualize the results of the conditional analysis or other interested SNPs.


## Interaction analysis (`Step8_Interaction_analysis/`)

This folder contains an R script for performing interaction analysis to examine whether sex-specific differences in effect sizes exist for genetic variants associated with food allergy.

1. Interaction analysis (`1_Interaction_analysis.R`) 
Performs interaction analysis by calculating the Z-scores by subtracting the male beta-value from the female beta-value, then dividing by the square root of the sum of the squared standard errors, and deriving two-tailed p-values from the Z-scores.


## Heritability estimation and chi-squared test (`Step9_Heritability_and_chisqr_test/`)

This folder contains scripts for heritability estimation and performing a chi-squared test for association between sex and case/control status.

1. Heritability estimation (`1_Heritability_GCTA_GREML.sh`)  
Estimates the heritability of food allergy in the German cohort from imputed data using the GREML-LDMS function from GCTA. Please find the details of this method [here](https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata).
2. Chi-squared test (`2_Chi_sqr_test.Rmd`)  
Performs a chi-squared test to determine if there is a significant association between sex and case/control status in the cohort.

Additionally, there is a functional R script `GREML_STEP2_stratify_SNPs_by_LD.R` used in step 2 of the GREML analysis to stratify the SNPs by LD scores of individual SNPs.




