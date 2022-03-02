rm(list = ls())
graphics.off()

library(data.table)

### IMPORT ==========
h_mtx <- fread("/h.sap_GRC38_v101_whole_genome_features_multicell_overlap.matrix")
m_mtx <- fread("/m.mus_GRC38_v101_whole_genome_features_multicell_overlap.matrix")

### FORMAT ==========
ann <- colnames(h_mtx)

h_mtx <- as.data.table(t(h_mtx))
m_mtx <- as.data.table(t(m_mtx))

h_mtx[] <- lapply(h_mtx, round, digits = 1)
m_mtx[] <- lapply(m_mtx, round, digits = 1)

colnames(h_mtx) <- ann
colnames(m_mtx) <- ann

h_mtx <- cbind(data.table(Species = rep("Human", nrow(h_mtx)), Annotation_A = colnames(h_mtx)), h_mtx)
m_mtx <- cbind(data.table(Species = rep("Mouse", nrow(m_mtx)), Annotation_A = colnames(m_mtx)), m_mtx)
output <- rbind(h_mtx, m_mtx)

### EXPORT ==========
fwrite(output, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_table_2_feature_overlap.csv")
