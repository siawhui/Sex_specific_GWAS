args <- commandArgs(trailingOnly = TRUE)

ld_file <- args[1]
output_name <- args[2]

# ld_file <- "W:/Siaw/Sex_specific_GWAS/conditional_analysis/heritability/GERMAN.score.ld"
# output_name <- "W:/Siaw/Sex_specific_GWAS/conditional_analysis/heritability/GERMAN"

lds_seg = read.table(ld_file,header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5]) # quartiles[4] is mean
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, paste0(output_name, "_snp_group1.txt"), row.names=F, quote=F, col.names=F)
write.table(lb2_snp, paste0(output_name, "_snp_group2.txt"), row.names=F, quote=F, col.names=F)
write.table(lb3_snp, paste0(output_name, "_snp_group3.txt"), row.names=F, quote=F, col.names=F)
write.table(lb4_snp, paste0(output_name, "_snp_group4.txt"), row.names=F, quote=F, col.names=F)