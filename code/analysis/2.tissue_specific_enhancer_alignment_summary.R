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
out.specific.file <- args[6]

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

### RUN FOR TISSUE SPECIFIC ANNNOTATION

human_enhancer <- human_ann[category == "Enhancer - distal" | category == "Enhancer - proximal"]
mouse_enhancer <- mouse_ann[category == "Enhancer - distal" | category == "Enhancer - proximal"]
human_enhancer_specific <- list(heart = human_enhancer[activity_heart == "ACTIVE" | activity_heart == "POISED"],
                                kidney = human_enhancer[activity_kidney == "ACTIVE" | activity_kidney == "POISED"],
                                spleen = human_enhancer[activity_spleen == "ACTIVE" | activity_spleen == "POISED"])
mouse_enhancer_specific <- list(heart = mouse_enhancer[activity_heart == "ACTIVE" | activity_heart == "POISED"],
                                kidney = mouse_enhancer[activity_kidney == "ACTIVE" | activity_kidney == "POISED"],
                                spleen = mouse_enhancer[activity_spleen == "ACTIVE" | activity_spleen == "POISED"])

out_list <- list()
for (j in 1:length(human_enhancer_specific)){
  
  # for each human annotation in chromosome:
  human_enhancer_sub <- human_enhancer_specific[[j]]
  mouse_enhancer_sub <- mouse_enhancer_specific[[j]]
  
  ann <- unique(human_enhancer_sub$category)
  n_total <- rep(NA, length(ann))
  n_aligned <- rep(NA, length(ann))
  n_conserved <- rep(NA, length(ann))
  for(i in 1:length(ann)){
    
    # collapse annotation to get total bases
    h_enhancer_colapse <- collapse.overlap(human_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_total[i] <- sum(abs( (h_enhancer_colapse$end + 1) - h_enhancer_colapse$start))
    
    # orthologous sequences for human annotation to get total alignment
    ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                      "chromosome_mouse", "start_mouse", "end_mouse")],
                                 bed = human_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_aligned[i] <- sum(abs((ann_align$end_A + 1) - ann_align$start_A))
    
    # total of mouse annotation in orthologous sequences for conservation
    m_enhancer_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                          bed2 = mouse_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_conserved[i] <- sum(abs( (m_enhancer_conserved$end + 1) - m_enhancer_conserved$start))
    
    print(ann[i])
  }
  
  out_list[[j]] <- data.table( tissue = names(human_enhancer_specific)[j],
                               chromosome = chr,
                               annotation = ann,
                               n_total = n_total,
                               n_aligned = n_aligned,
                               n_conserved = n_conserved)
  
  print(j)
}

out_dt_specific <- do.call("rbind", out_list)
# out_dt_specific$p_aligned <- out_dt_specific$n_aligned / out_dt_specific$n_total
# out_dt_specific$p_conserved <- out_dt_specific$n_conserved / out_dt_specific$n_total


### EXPORT
fwrite(out_dt_specific, out.specific.file, append = T)



