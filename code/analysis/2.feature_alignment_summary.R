rm(list=ls())
graphics.off()

library(data.table)

### SET VARS 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
h.ann.file <- args[1]
m.ann.file <- args[2]
align.file <- args[3]
chr <- as.character(args[4])
functions.file <- args[5]
out.all.file <- args[6]

### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))
align <- fread(paste0("gunzip -cq ", align.file))

### FORMAT
human_ann$category[human_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"
mouse_ann$category[mouse_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"

# subset by human chromosome
human_ann <- human_ann[chromosome == chr]
align <- align[chromosome_human == chr]

# QC
# human_ann <- human_ann[percentage_Nmask < 1]
# mouse_ann <- mouse_ann[percentage_Nmask < 1]


### RUN FOR ALL ANNNOTATION

# for each human annotation in chromosome:
ann <- unique(human_ann$category)
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  h_ann_colapse <- collapse.overlap(human_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_total[i] <- sum(abs( (h_ann_colapse$end + 1) - h_ann_colapse$start))
  
  # orthologous sequences for human annotation to get total alignment
  ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                               "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = human_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_aligned[i] <- sum(abs( (ann_align$end_A + 1) - ann_align$start_A))

  # total of mouse annotation in orthologous sequences for conservation
  m_ann_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = mouse_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_conserved[i] <- sum(abs ((m_ann_conserved$end + 1) - m_ann_conserved$start))
  
  print(ann[i])
}

out_dt_all <- data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved)
# out_dt_all$p_aligned <- out_dt_all$n_aligned / out_dt_all$n_total
# out_dt_all$p_conserved <- out_dt_all$n_conserved / out_dt_all$n_total
# rm(ann_align, h_ann_colapse, m_ann_conserved, ann, i, n_aligned, n_conserved, n_total)


### RUN FOR "Functional - all" "Functional - non-CDS" "Non-CDS" and "Total"

h_cds <- collapse.overlap(human_ann[category %in% c("Exon - CDS")][,c("chromosome", "start", "end")])
h_fun_all <- collapse.overlap(human_ann[!category %in% c("Intron - distal", "Unannotated")][,c("chromosome", "start", "end")])
h_fun_nonCDS <- bed_antijoin(bed1 = h_fun_all, bed2 = h_cds)
h_total <- collapse.overlap(human_ann[,c("chromosome", "start", "end")])
h_nonCDS <- bed_antijoin(bed1 = h_total, bed2 = h_cds)
human_ann_2_list <- list(functional = h_fun_all,
                    fun_non_cds =  h_fun_nonCDS,
                    non_cds = h_nonCDS,
                    total = human_ann)

m_cds <- collapse.overlap(mouse_ann[category %in% c("Exon - CDS")][,c("chromosome", "start", "end")])
m_fun_all <- collapse.overlap(mouse_ann[!category %in% c("Intron - distal", "Unannotated")][,c("chromosome", "start", "end")])
m_fun_nonCDS <- bed_antijoin(bed1 = m_fun_all, bed2 = m_cds)
m_total <- collapse.overlap(mouse_ann[,c("chromosome", "start", "end")])
m_nonCDS <- bed_antijoin(bed1 = m_total, bed2 = m_cds)
mouse_ann_2_list <- list(functional = m_fun_all,
                         fun_non_cds =  m_fun_nonCDS,
                         non_cds = m_nonCDS,
                         total = mouse_ann)

ann <- c("Functional - all", "Functional - nonCDS", "NonCDS", "Total")
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  human_ann_2 <- human_ann_2_list[[i]]
  mouse_ann_2 <- mouse_ann_2_list[[i]]
  
  h_ann_colapse <- collapse.overlap(human_ann_2[,c("chromosome", "start", "end")])
  n_total[i] <- sum(abs( (h_ann_colapse$end + 1) - h_ann_colapse$start))
  
  # orthologous sequences for human annotation to get total alignment
  ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = h_ann_colapse)
  n_aligned[i] <- sum(abs((ann_align$end_A + 1) - ann_align$start_A))
  sum(abs( (ann_align$end_B + 1) - ann_align$start_B))
  
  # total of mouse annotation in orthologous sequences for conservation
  # m_ann_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
  #                                  bed2 = mouse_ann_2[,c("chromosome", "start", "end")])
  # n_conserved[i] <- sum(abs( (m_ann_conserved$end + 1) - m_ann_conserved$start))
  
  print(ann[i])
}

out_dt_2 <- data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved)
# out_dt_2$p_aligned <- out_dt_2$n_aligned / out_dt_2$n_total
# out_dt_2$p_conserved <- out_dt_2$n_conserved / out_dt_2$n_total


### EXPORT

out_all <- rbind(out_dt_all, out_dt_2)
fwrite(out_all, out.all.file, append = T)


#####
