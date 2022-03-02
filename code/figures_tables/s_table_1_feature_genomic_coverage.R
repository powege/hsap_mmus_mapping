rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### SET VARS
functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
h.ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
m.ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_1_multicell_annotation_gennomic_coverage.csv"

### FUNCTIONS
source(functions.file)

### IMPORT
m_ann <- fread(paste0("gunzip -cq ", m.ann.file))
h_ann <- fread(paste0("gunzip -cq ", h.ann.file))

### QC 
h_ann$category[h_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"
m_ann$category[m_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"

# m_ann <- m_ann[percentage_Nmask < 1]
# h_ann <- h_ann[percentage_Nmask < 1]

### FORMAT
m_ann_total <- collapse.overlap(m_ann[,1:3])
m_ann_total <- sum( (m_ann_total$end +1) - m_ann_total$start)
  
h_ann_total <- collapse.overlap(h_ann[,1:3])
h_ann_total <- sum( (h_ann_total$end +1) - h_ann_total$start)

ann <- unique(m_ann$category)
out_list <- list()
for (i in 1:length(ann)){
  m_sub <- collapse.overlap(m_ann[category == ann[i]][,1:3])
  h_sub <- collapse.overlap(h_ann[category == ann[i]][,1:3])
  out_list[[i]] <- data.table(category = ann[i], 
                              n_bp_mouse = sum( (m_sub$end +1) - m_sub$start),
                              n_bp_human = sum( (h_sub$end +1) - h_sub$start))
}
output1 <- do.call("rbind", out_list)

### RUN FOR "Functional - all" "Functional - non-CDS" and "Total"

h_cds <- collapse.overlap(h_ann[category %in% c("Exon - CDS")][,c("chromosome", "start", "end")])
h_fun_all <- collapse.overlap(h_ann[!category %in% c("Intron - distal", "Unannotated")][,c("chromosome", "start", "end")])
h_fun_nonCDS <- bed_antijoin(bed1 = h_fun_all, bed2 = h_cds)
h_ann_2_list <- list(functional = h_fun_all,
                         non_cds =  h_fun_nonCDS)

m_cds <- collapse.overlap(m_ann[category %in% c("Exon - CDS")][,c("chromosome", "start", "end")])
m_fun_all <- collapse.overlap(m_ann[!category %in% c("Intron - distal", "Unannotated")][,c("chromosome", "start", "end")])
m_fun_nonCDS <- bed_antijoin(bed1 = m_fun_all, bed2 = m_cds)
m_ann_2_list <- list(functional = m_fun_all,
                         non_cds =  m_fun_nonCDS)

ann <- c("Functional - all", "Functional - nonCDS")
out_list <- list()
for (i in 1:length(ann)){
  m_sub <- collapse.overlap(m_ann_2_list[[i]][,1:3])
  h_sub <- collapse.overlap(h_ann_2_list[[i]][,1:3])
  out_list[[i]] <- data.table(category = ann[i], 
                              n_bp_mouse = sum( (m_sub$end +1) - m_sub$start),
                              n_bp_human = sum( (h_sub$end +1) - h_sub$start))
}
output2 <- do.call("rbind", out_list)

output <- rbind(output1, output2)
output$total_bp_mouse <- m_ann_total
output$total_bp_human <- h_ann_total
output$genomic_coverage_mouse <- round( (output$n_bp_mouse / m_ann_total)*100, 1)
output$genomic_coverage_human <- round( (output$n_bp_human / h_ann_total)*100, 1)
output$genomic_coverage_mouse <- format(output$genomic_coverage_mouse, nsmall = 1)
output$genomic_coverage_human <- format(output$genomic_coverage_human, nsmall = 1)

### EXPORT
fwrite(output, out.file)


#####
