rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)
library(RColorBrewer)

### SET VARS

in.clinvar <- "/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_noCDSoverlap_by_chr.csv"
in.gwas <- "/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_noCDSoverlap_by_chr.csv"
out.plot.table <- "/Figure_2_hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_GWAS_SNV_alignment_noCDSoverlap.csv"
out.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_figure_3_hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_GWAS_SNV_alignment_noCDSoverlap.jpeg"


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

cv$group <- "ClinVar"
gwas$group <- "GWAS"
out_table <- rbind(cv, gwas)

### QC
out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "ClinVar"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "ClinVar"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "ClinVar"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "ClinVar"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "ClinVar"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "ClinVar"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "ClinVar"]

out_table$n_total[out_table$annotation == "Exon - CDS" & out_table$group == "GWAS"] +
  out_table$n_total[out_table$annotation == "Functional - nonCDS" & out_table$group == "GWAS"] ==
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "GWAS"]

out_table$n_total[out_table$annotation == "Intron - distal" & out_table$group == "GWAS"] +
  out_table$n_total[out_table$annotation == "Unannotated" & out_table$group == "GWAS"] +
  out_table$n_total[out_table$annotation == "Functional - all" & out_table$group == "GWAS"] ==
  out_table$n_total[out_table$annotation == "Total" & out_table$group == "GWAS"]

### FORMAT FIGURES

order <- c("Exon - CDS",
           "Exon - UTR",
           "Exon - other",
           "Promoter",
           "Intron - proximal",
           "Enhancer - proximal",
           "Enhancer - distal",
           "CTCF binding",
           "TAD boundry",
           "Miscellaneous",
           "Intron - distal",
           "Unannotated")
order <- rev(order)
GW_cv_alignment <- cv$percentage_aligned[cv$annotation == "Total"]
GW_gwas_alignment <- gwas$percentage_aligned[gwas$annotation == "Total"]

cv_plot <- cv[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
cv_plot$x_lab <- paste0(cv_plot$annotation, "\n(n=", cv_plot$n_total, ")")
gwas_plot <- gwas[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
gwas_plot$x_lab <- paste0(gwas_plot$annotation, "\n(n=", gwas_plot$n_total, ")")

cv_plot$annotation <- factor(cv_plot$annotation, levels = as.character(order))
cv_plot <- cv_plot[order(annotation)]
lab_order <- cv_plot$x_lab
cv_plot$x_lab <- factor(cv_plot$x_lab, levels = as.character(lab_order))
gwas_plot$annotation <- factor(gwas_plot$annotation, levels = as.character(order))
gwas_plot <- gwas_plot[order(annotation)]
lab_order <- gwas_plot$x_lab
gwas_plot$x_lab <- factor(gwas_plot$x_lab, levels = as.character(lab_order))

cv_plot$percentage_not_conserved <- cv_plot$percentage_aligned - cv_plot$percentage_conserved
cv_plot <- melt(cv_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
cv_plot$variable <- factor(cv_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
cv_plot$col <- ifelse(cv_plot$variable == "percentage_not_conserved", "black", as.character(cv_plot$annotation))
cv_plot$col <- factor(cv_plot$col, levels = c("black", as.character(order)))
gwas_plot$percentage_not_conserved <- gwas_plot$percentage_aligned - gwas_plot$percentage_conserved
gwas_plot <- melt(gwas_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
gwas_plot$variable <- factor(gwas_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
gwas_plot$col <- ifelse(gwas_plot$variable == "percentage_not_conserved", "black", as.character(gwas_plot$annotation))
gwas_plot$col <- factor(gwas_plot$col, levels = c("black", as.character(order)))

cv_color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
cv_color_map <- setNames(cv_color_scale, levels(unique(cv_plot$col)))
gwas_color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
gwas_color_map <- setNames(gwas_color_scale, levels(unique(gwas_plot$col)))
out_table <- rbind(cv_plot, gwas_plot)

### PLOT

p_cv <- ggplot(cv_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% SNV sites)") +
  ggtitle("Mendelian") +
  scale_fill_manual(values = cv_color_map) +
  geom_hline(yintercept=GW_cv_alignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_cv

p_gwas <- ggplot(gwas_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% SNV sites)") +
  ggtitle("Complex") +
  scale_fill_manual(values = gwas_color_map) +
  geom_hline(yintercept=GW_gwas_alignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_gwas

### EXPORT

pout <- grid.arrange(p_cv, p_gwas, nrow = 2, ncol = 1)
ggsave(out.figure, plot = pout, height = 12, width = 7)
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

