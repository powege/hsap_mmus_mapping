rm(list = ls())
graphics.off()

library(data.table)
library(xtable)

### SET VARS 

in.multicell.file <- "/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary_by_chr.csv"
out.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_3_hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment.csv"

### IMPORT
multi_raw <- fread(in.multicell.file)

### FORMAT TABLE 
multi_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(multi_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(multi_raw$chromosome))
multi <- setDT(multi_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation)]
multi$percentage_aligned <- round( (multi$n_aligned / multi$n_total)*100, 1)
multi$percentage_conserved <- round( (multi$n_conserved / multi$n_total)*100, 1)

### QC
multi$n_total[multi$annotation == "Exon - CDS"] +
  multi$n_total[multi$annotation == "Functional - nonCDS"] ==
  multi$n_total[multi$annotation == "Functional - all"]

multi$n_total[multi$annotation == "Intron - distal"] +
  multi$n_total[multi$annotation == "Unannotated"] +
  multi$n_total[multi$annotation == "Functional - all"] ==
  multi$n_total[multi$annotation == "Total"]

multi$n_total[multi$annotation == "Exon - CDS"] +
  multi$n_total[multi$annotation == "NonCDS"] ==
  multi$n_total[multi$annotation == "Total"]

### CREATE TABLE
colnames(multi) <- c("Annotation",
                     "N human bp",
                     "N bp aligned to mouse",
                     "N bp conserved in mouse",
                     "Bp aligned to mouse (%)",
                     "Bp conserved in mouse (%)")
temp.table <- xtable(multi, digits = 1)


### EXPORT
fwrite(multi, out.table)



