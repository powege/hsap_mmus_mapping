rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)
library(RColorBrewer)

### SET VARS 

in.multicell.file <- "/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary_noCDSoverlap_by_chr.csv"
out.multicell.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_figure_2_hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_stacked_noCDSoverlap.jpg"
out.plot.table <- "/s_figure_2_hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_stacked_noCDSoverlap.csv"


### IMPORT
multi_raw <- fread(in.multicell.file)

### FORMAT TABLE 
multi_raw[, c("n_total", "n_aligned", "n_conserved")] <- lapply(multi_raw[, c("n_total", "n_aligned", "n_conserved")], as.numeric)
length(unique(multi_raw$chromosome))
multi <- setDT(multi_raw)[, .(n_total = sum(n_total), 
                              n_aligned = sum(n_aligned),
                              n_conserved = sum(n_conserved)), by = .(annotation)]
multi$percentage_aligned <- round( (multi$n_aligned / multi$n_total)*100, 2)
multi$percentage_conserved <- round( (multi$n_conserved / multi$n_total)*100, 2)
multi$genomic_coverage <- round( (multi$n_total/multi$n_total[multi$annotation == "Total"])*100, 1)

### QC
multi$n_total[multi$annotation == "Exon - CDS"] +
  multi$n_total[multi$annotation == "Functional - nonCDS"] ==
  multi$n_total[multi$annotation == "Functional - all"]

multi$n_total[multi$annotation == "Intron - distal"] +
  multi$n_total[multi$annotation == "Unannotated"] +
  multi$n_total[multi$annotation == "Functional - all"] ==
  multi$n_total[multi$annotation == "Total"]

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
GW_Halignment <- multi$percentage_aligned[multi$annotation == "Total"]

multi_plot <- multi[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
multi_plot$x_lab <- paste0(multi_plot$annotation, "\n(", format(multi_plot$genomic_coverage, nsmall = 1), "%)")

multi_plot$annotation <- factor(multi_plot$annotation, levels = as.character(order))
multi_plot <- multi_plot[order(annotation)]
lab_order <- multi_plot$x_lab
multi_plot$x_lab <- factor(multi_plot$x_lab, levels = as.character(lab_order))

multi_plot$percentage_not_conserved <- multi_plot$percentage_aligned - multi_plot$percentage_conserved
multi_plot <- melt(multi_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
multi_plot$variable <- factor(multi_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
multi_plot$col <- ifelse(multi_plot$variable == "percentage_not_conserved", "black", as.character(multi_plot$annotation))
multi_plot$col <- factor(multi_plot$col, levels = c("black", as.character(order)))

color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
color_map <- setNames(color_scale, levels(unique(multi_plot$col)))

### PLOT

p_multi <- ggplot(multi_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% bp)") +
  # ggtitle("B") +
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_multi

### EXPORT
ggsave(out.multicell.figure, plot = p_multi, height = 6, width = 7)
fwrite(multi_plot, out.plot.table)



