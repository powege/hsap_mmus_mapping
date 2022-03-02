rm(list = ls())
graphics.off()

library(data.table)
library(xtable)

### SET VARS 

functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
h.ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
in.multicell.file <- "/hsap_grch38_v_mmus_grcm38_v101_tissue_feature_alignment_summary_by_chr.csv"
out.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_4_hsap_grch38_v_mmus_grcm38_v101_tissue_specific_enhancer_alignment.csv"

### FUNCTIONS
source(functions.file)

### IMPORT
h_ann <- fread(paste0("gunzip -cq ", h.ann.file))
multi_raw <- fread(in.multicell.file)

### FORMAT TABLE 
h_ann_total <- collapse.overlap(h_ann[,1:3])
h_ann_total <- sum( (h_ann_total$end +1) - h_ann_total$start)

multi_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(multi_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(multi_raw$chromosome))
multi <- setDT(multi_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation, tissue)]
multi$percentage_aligned <- round( (multi$n_aligned / multi$n_total)*100, 1)
multi$percentage_conserved <- round( (multi$n_conserved / multi$n_total)*100, 1)

multi$genomic_coverage_human <- round( (multi$n_total / h_ann_total)*100, 1)

### CREATE TABLE
colnames(multi) <- c("Annotation",
                     "Tissue",
                     "N human bp",
                     "N bp aligned to mouse",
                     "N bp conserved in mouse",
                     "Bp aligned with mouse (%)",
                     "Bp conserved in mouse (%)",
                     "Human genomic coverage (%)")
temp.table <- xtable(multi, digits = 1)


### EXPORT
fwrite(multi, out.table)



