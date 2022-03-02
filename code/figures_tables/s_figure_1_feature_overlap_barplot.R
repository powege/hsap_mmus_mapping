rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


### IMPORT ==========
h_mtx <- fread("/h.sap_GRC38_v101_whole_genome_features_multicell_overlap.matrix")
m_mtx <- fread("/m.mus_GRC38_v101_whole_genome_features_multicell_overlap.matrix")

### FORMAT ==========
h_mtx[] <- lapply(h_mtx, as.numeric)
m_mtx[] <- lapply(m_mtx, as.numeric)

h_mtx <- cbind(data.table(annotation = colnames(h_mtx)), h_mtx)
m_mtx <- cbind(data.table(annotation = colnames(m_mtx)), m_mtx)

h_dt <- melt(h_mtx, id.vars = "annotation",
             measure.vars = colnames(h_mtx)[2:ncol(h_mtx)])
m_dt <- melt(m_mtx, id.vars = "annotation",
             measure.vars = colnames(m_mtx)[2:ncol(m_mtx)])

# set factor order
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
h_dt$annotation <- factor(h_dt$annotation, levels = as.character(order))
h_dt$variable <- factor(h_dt$variable, levels = as.character(order))
m_dt$annotation <- factor(m_dt$annotation, levels = as.character(order))
m_dt$variable <- factor(m_dt$variable, levels = as.character(order))

m_dt <- m_dt[value != 100]
h_dt <- h_dt[value != 100]

### PLOT Human bar

pC <- ggplot() +
  geom_bar(data=h_dt, aes(x=variable, y=value, fill=annotation), colour = "black", stat="identity") +
  xlab("Annotation A") +
  ylab("Overlap with annotation B (%)") +
  ggtitle("Human") +
  # labs(fill = "") +
  # scale_y_continuous(breaks = c(0, 20, 40, 60, 80),
  #                    # trans = "reverse",
  #                    limits = c(0, 85)) +
  coord_flip() +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    # legend.title = element_text(hjust = 0.5),
    # legend.key.size = unit(2, 'lines'),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pC


### PLOT Mouse bar

pD <- ggplot() +
  geom_bar(data=m_dt, aes(x=variable, y=value, fill=annotation), colour = "black", stat="identity") +
  xlab("Annotation A") +
  ylab("Overlap with annotation B (%)") +
  ggtitle("Mouse") +
  # labs(fill = "") +
  # scale_y_continuous(breaks = c(0, 20, 40, 60, 80),
  #                    # trans = "reverse",
  #                    limits = c(0, 85)) +
  coord_flip() +
  scale_fill_brewer(palette="Set3", name = "Annotation B") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    # legend.position = "none",
    # legend.title = element_text(hjust = 0.5),
    # legend.key.size = unit(2, 'lines'),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pD


### EXPORT ==========
poutCD <- grid.arrange(pC, pD, nrow = 1, widths = c(5, 7))
ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Manuscript/Figures_and_tables/s_figure_1_feature_overlap.jpg", plot = poutCD, height = 5, width = 12)

