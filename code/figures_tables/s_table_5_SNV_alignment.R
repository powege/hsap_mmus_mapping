rm(list = ls())
graphics.off()

library(data.table)


### SET VARS
in.clinvar <- "/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_by_chr.csv"
in.gwas <- "/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv"
out.plot.table <- "/Figure_2_hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_GWAS_SNV_alignment.csv"
out.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_4_hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_GWAS_SNV_alignment.csv"

### IMPORT

cv_raw <- fread(in.clinvar)
gwas_raw <- fread(in.gwas)

### FORMAT TABLE

cv_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(cv_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(cv_raw$chromosome))
cv <- setDT(cv_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation)]
gwas_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(gwas_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(gwas_raw$chromosome))
gwas <- setDT(gwas_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation)]

cv$n_conserved[cv$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA
gwas$n_conserved[gwas$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA

cv$percentage_aligned <- (cv$n_aligned / cv$n_total)*100
cv$percentage_conserved <- (cv$n_conserved / cv$n_total)*100
gwas$percentage_aligned <- (gwas$n_aligned / gwas$n_total)*100
gwas$percentage_conserved <- (gwas$n_conserved / gwas$n_total)*100

cv$group <- "Mendelian"
gwas$group <- "Complex"
out_table <- rbind(cv, gwas)

### QC
out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Mendelian"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "Mendelian"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Mendelian"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "Mendelian"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "Mendelian"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Mendelian"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Mendelian"]

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Mendelian"] +
  out_table$n_total[out_table$annotation == "NonCDS" & out_table$group == "Mendelian"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Mendelian"]

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Complex"]

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "NonCDS" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Complex"]


### CREATE TABLE

out_table <- out_table[,c("group", "annotation", "n_total", "n_aligned", "n_conserved", "percentage_aligned", "percentage_conserved")]
colnames(out_table) <- c("Disease class",
                         "Annotation",
                         "N SNV",
                         "N SNV aligned",
                         "N SNV conserved",
                         "SNV aligned (%)",
                         "SNV conserved (%)")
temp.table <- xtable(out_table, digits = 1)
temp.table$`SNV aligned (%)` <- round(temp.table$`SNV aligned (%)`, 1)
temp.table$`SNV conserved (%)` <- round(temp.table$`SNV conserved (%)`, 1)


### EXPORT

fwrite(temp.table, out.table)
fwrite(out_table, out.plot.table)


#####

# cv_raw <- fread(in.clinvar)
# cv_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(cv_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
# length(unique(cv_raw$chromosome))
# for (i in 1:22){
#   cv_sub <- cv_raw[chromosome == i]
#   
#   print(cv_sub$n_total[cv_sub$annotation == "Exon - CDS"] +
#     cv_sub$n_total[cv_sub$annotation == "Functional - nonCDS"] ==
#     cv_sub$n_total[cv_sub$annotation == "Functional - all"])
#   
#   print(cv_sub$n_total[cv_sub$annotation == "Intron - distal"] +
#     cv_sub$n_total[cv_sub$annotation == "Unannotated"] +
#     cv_sub$n_total[cv_sub$annotation == "Functional - all"] ==
#     cv_sub$n_total[cv_sub$annotation == "Total"])
#   
#   print(i)
# }

