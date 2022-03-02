# https://www.statisticshowto.com/probability-and-statistics/hypothesis-testing/z-test/
# https://www.cyclismo.org/tutorial/R/pValues.html

rm(list = ls())
graphics.off()

library(data.table)


### SET VARS
in.gwas <- "/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv"
in.multicell.file <- "/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary_by_chr.csv"
out.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_8_hsap_v_mmus_RegBuild_v101_SNV_alignment_mendelian_total_zscores.csv"

### IMPORT

gwas_raw <- fread(in.gwas)
multi_raw <- fread(in.multicell.file)

### FORMAT TABLE

gwas_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(gwas_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(gwas_raw$chromosome))
gwas <- setDT(gwas_raw)[, .(n_total = sum(n_total), 
                            n_aligned = sum(n_aligned),
                            n_conserved = sum(n_conserved)), by = .(annotation)]
multi_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(multi_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(multi_raw$chromosome))
multi <- setDT(multi_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation)]

gwas$n_conserved[gwas$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA
multi$n_conserved[multi$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA

gwas$group <- "Complex"
multi$group <- "Total"
out_table <- rbind(gwas, multi)

### QC
out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "Complex"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Complex"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Complex"]

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "Total"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "Total"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Total"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "Total"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "Total"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "Total"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "Total"]

ann <- unique(out_table$annotation)
a_z_vec <- rep(NA, length(ann))
a_p_vec <- rep(NA, length(ann))
for(i in seq_along(ann)){
  
  a1 = out_table$n_aligned[out_table$annotation == ann[i] & out_table$group == "Complex"]
  a2 = out_table$n_aligned[out_table$annotation == ann[i] & out_table$group == "Total"]
  b1 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Complex"]
  b2 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Total"]
  
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
  
  a1 = out_table$n_conserved[out_table$annotation == ann[i] & out_table$group == "Complex"]
  a2 = out_table$n_conserved[out_table$annotation == ann[i] & out_table$group == "Total"]
  b1 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Complex"]
  b2 = out_table$n_total[out_table$annotation == ann[i] & out_table$group == "Total"]
  
  p1 <- a1/b1
  p2 <- a2/b2
  p <- (a1 + a2)/(b1 + b2)
  
  z <- ( (p1 - p2) - 0) / sqrt((p*(1-p)) * ((1/b1)+(1/b2)) )
  pvalue2sided=2*pnorm(-abs(z))
  
  c_z_vec[i] <- z
  c_p_vec[i] <- pvalue2sided
}

a_z_vec <- round(a_z_vec, digits = 1)
a_p_vec <- format(a_p_vec, digits = 2)
c_z_vec <- round(c_z_vec, digits = 1)
c_p_vec <- format(c_p_vec, digits = 2)

out_z <- data.table(Annotation = ann, 
                    `Alignment z-score` = a_z_vec,
                    `Alignment p-value` = a_p_vec,
                    `Conservation z-score` = c_z_vec,
                    `Conservation p-value` = c_p_vec)


### CREATE TABLE
temp.table <- xtable(out_z, digits = 1)


### EXPORT
fwrite(out_z, out.table)

