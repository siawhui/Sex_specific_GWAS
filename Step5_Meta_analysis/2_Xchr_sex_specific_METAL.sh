#!/bin/bash

#======== reformat chrX table =======#

datasets=("GERMAN" "SPAIN")
sexes=("female" "male")

for dataset in "${datasets[@]}"; do
   for sex in "${sexes[@]}"; do

      awk '{OFS="\t"} NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P"; next}
        {
           MarkerID = $1":"$2":"$4":"$3;
           CHR = $1;
           BP = $2;
           A1 = $3;
           A2 = $4;
           Beta = $5;
           P = $6;
           print MarkerID, CHR, BP, A1, A2, Beta, P
       }' xwas_${dataset}.02.xstrat_${sex}.logistic > xwas_${dataset}_${sex}.txt

   done
done

#========= run METAL =========#
metal METAL_script_XWAS_female.txt
metal METAL_script_XWAS_male.txt

# Sort
sort -k6 -g METAANALYSIS_XWAS_female1.tbl > METAANALYSIS_XWAS_female_sorted.tbl
sort -k6 -g METAANALYSIS_XWAS_male1.tbl > METAANALYSIS_XWAS_male_sorted.tbl

# FUMA format
awk '{OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "P"; next}
  {
     split($1, arr, ":");
     CHR = arr[1];
     BP = arr[2];
     A2 = toupper($3); #ref
     A1 = toupper($2); #alt
     Beta = $5;
     P = $6;
     print CHR, BP, A2, A1, Beta, P
}' METAANALYSIS_XWAS_female_sorted.tbl > METAANALYSIS_XWAS_female_FUMA.txt

awk '{OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "P"; next}
  {
     split($1, arr, ":");
     CHR = arr[1];
     BP = arr[2];
     A2 = toupper($3); #ref
     A1 = toupper($2); #alt
     Beta = $5;
     P = $6;
     print CHR, BP, A2, A1, Beta, P
}' METAANALYSIS_XWAS_male_sorted.tbl > METAANALYSIS_XWAS_male_FUMA.txt

# GWASLab format
awk '{OFS="\t"} NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P"; next}
  {
     MarkerID = $1":"$2":"$4":"$3;
     CHR = $1;
     BP = $2;
     A1 = $3;
     A2 = $4;
     Beta = $5;
     P = $6;
     print MarkerID, CHR, BP, A1, A2, Beta, P
 }' METAANALYSIS_XWAS_female_FUMA.txt > METAANALYSIS_XWAS_female_GWASLab.txt

awk '{OFS="\t"} NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P"; next}
  {
     MarkerID = $1":"$2":"$4":"$3;
     CHR = $1;
     BP = $2;
     A1 = $3;
     A2 = $4;
     Beta = $5;
     P = $6;
     print MarkerID, CHR, BP, A1, A2, Beta, P
 }' METAANALYSIS_XWAS_male_FUMA.txt > METAANALYSIS_XWAS_male_GWASLab.txt
