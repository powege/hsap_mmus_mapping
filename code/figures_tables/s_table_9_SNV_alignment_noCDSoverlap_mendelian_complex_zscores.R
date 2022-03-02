# https://www.statisticshowto.com/probability-and-statistics/hypothesis-testing/z-test/
# https://www.cyclismo.org/tutorial/R/pValues.html

rm(list = ls())
graphics.off()

library(data.table)


### SET VARS
in.clinvar <- "/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_noCDSoverlap_by_chr.csv"
in.gwas <- "/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_noCDSoverlap_by_chr.csv"
out.table <- "/s_table_5_hsap_v_mmus_RegBuild_v101_SNV_alignment_noCDSoverlap_mendelian_complex_zscores.csv"

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

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Complex"]

ann <- unique(out_table$annotation)
a_z_vec <- rep(NA, length(ann))
a_p_vec <- rep(NA, length(ann))
for(i in seq_along(ann)){
  
  a1 = out_table$n_aligned[out_table$annotation == ann[i] & out_table$group == "Mendelian"]
  a2 = out_table$n_aligned[out_table$annotation == ann[i] & out_table$group == "Complex"]
  b1 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Mendelian"]
  b2 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Complex"]
  
  p1 <- a1/b1
  p2 <- a2/b2
  p <- (a1 + a2)/(b1 + b2)
  
  z <- ( (p1 - p2) - 0) / sqrt((p*(1-p)) * ((1/b1)+(1/b2)) )
  pvalue2sided=2*pnorm(-abs(z))
  
  a_z_vec[i] <- z
  a_p_vec[i] <- pvalue2sided
}

c_z_vec <- rep(NA, length(ann))
c_p_vec <- rep(NA, length(ann))
for(i in seq_along(ann)){
  
  a1 = out_table$n_conserved[out_table$annotation == ann[i] & out_table$group == "Mendelian"]
  a2 = out_table$n_conserved[out_table$annotation == ann[i] & out_table$group == "Complex"]
  b1 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Mendelian"]
  b2 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Complex"]
  
  p1 <- a1/b1
  p2 <- a2/b2
  p <- (a1 + a2)/(b1 + b2)
  
  z <- ( (p1 - p2) - 0) / sqrt((p*(1-p)) * ((1/b1)+(1/b2)) )
  pvalue2sided=2*pnorm(-abs(z))
  
  c_z_vec[i] <- z
  c_p_vec[i] <- pvalue2sided
}

out_z <- data.table(Annotation = ann, 
                    `Alignment z-score` = a_z_vec,
                    `Alignment p-value` = a_p_vec,
                    `Conservation z-score` = c_z_vec,
                    `Conservation p-value` = c_p_vec)


### CREATE TABLE
temp.table <- xtable(out_z, digits = 1)


### EXPORT
# fwrite(out_z, out.table)

