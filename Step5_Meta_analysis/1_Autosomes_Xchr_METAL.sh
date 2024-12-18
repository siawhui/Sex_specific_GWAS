#!/bin/bash

mkdir -p results

#======== reformat XCMAX4 result table =======#

awk '{OFS="\t"} NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P", "SE"; next}
  {
     MarkerID = $1":"$2":"$4":"$3;
     CHR = $1;
     BP = $2;
     A1 = $3;
     A2 = $4;
     Beta = $5;
     P = $6;
     SE = $7;
     print MarkerID, CHR, BP, A1, A2, Beta, P, SE
 }' XCMAX4_GERMAN/results.txt > XCMAX4_GERMAN/results_METAL.txt

awk '{OFS="\t"} NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P", "SE"; next}
  {
     MarkerID = $1":"$2":"$4":"$3;
     CHR = $1;
     BP = $2;
     A1 = $3;
     A2 = $4;
     Beta = $5;
     P = $6;
     SE = $7;
     print MarkerID, CHR, BP, A1, A2, Beta, P, SE
 }' XCMAX4_SPAIN/results.txt > XCMAX4_SPAIN/results_METAL.txt

#========= run METAL =========#
metal METAL_script_XCMAX4.txt
metal METAL_script_all.txt
metal METAL_script_female.txt
metal METAL_script_male.txt

# Sort
sort -k6 -g results/METAANALYSIS_XCMAX41.tbl > results/METAANALYSIS_XCMAX4_sorted.tbl
sort -k10 -g results/METAANALYSIS_all1.tbl > results/METAANALYSIS_all_sorted.tbl
sort -k10 -g results/METAANALYSIS_female1.tbl > results/METAANALYSIS_female_sorted.tbl
sort -k10 -g results/METAANALYSIS_male1.tbl > results/METAANALYSIS_male_sorted.tbl

# FUMA format
for file in results/METAANALYSIS_all_sorted.tbl results/METAANALYSIS_female_sorted.tbl results/METAANALYSIS_male_sorted.tbl; do
    output_file="${file%_sorted.tbl}_FUMA.txt"  # Create output file name by appending _FUMA
    awk '{OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "SE", "P"; next}
     {
        split($1, arr, ":");
        CHR = arr[1];
        BP = arr[2];
        A2 = toupper($2);
        A1 = toupper($3);
        Beta = $8;
        SE = $9;
        P = $10;
        print CHR, BP, A2, A1, Beta, SE, P
    }' "$file" > "$output_file"
done

awk '{OFS="\t"} NR==1 {print "CHR", "BP", "A2", "A1", "Beta", "SE", "P"; next}
  {
     split($1, arr, ":");
     CHR = arr[1];
     BP = arr[2];
     A2 = toupper($2); #ref
     A1 = toupper($3); #alt
     Beta = $4;
     SE = $5;
     P = $6;
     print CHR, BP, A2, A1, Beta, SE, P
 }' results/METAANALYSIS_XCMAX4_sorted.tbl > results/METAANALYSIS_XCMAX4_FUMA.txt

# GWASLab format
# Filter out SNPs that are missing in one dataset or have opposite direction
for file in all female male; do
   awk '{OFS="\t"} 
   NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P"; next}
   $11 == "++" || $11 == "--" {
       MarkerID = $1;
       split($1, arr, ":");
       CHR = arr[1];
       BP = arr[2];
       A1 = toupper($2);
       A2 = toupper($3);
       Beta = $8;
       SE = $9;
       P = $10;
       print MarkerID, CHR, BP, A1, A2, Beta, P
   }' results/METAANALYSIS_${file}_sorted.tbl > results/METAANALYSIS_${file}_GWASLab.txt
done

awk '{OFS="\t"} 
NR==1 {print "MARKER", "CHR", "BP", "A1", "A2", "Beta", "P"; next}
$7 == "++" || $7 == "--" {
   MarkerID = $1;
   split($1, arr, ":");
   CHR = arr[1];
   BP = arr[2];
   A2 = toupper($3); #ref
   A1 = toupper($2); #alt
   Beta = $4;
   SE = $5;
   P = $6;
   print MarkerID, CHR, BP, A1, A2, Beta, P
 }' results/METAANALYSIS_XCMAX4_sorted.tbl > results/METAANALYSIS_XCMAX4_GWASLab.txt

sed -i 's/X/23/g' results/METAANALYSIS_XCMAX4_GWASLab.txt
