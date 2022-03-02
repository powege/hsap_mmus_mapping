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
h.mendelian.file <- args[4]
h.gwas.file <- args[5]
chr <- as.character(args[6])
functions.file <- args[7]
out.mendelian.file <- args[8]
out.gwas.file <- args[9]


### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))
align <- fread(paste0("gunzip -cq ", align.file))
mendel <- fread(h.mendelian.file)
gwas <- fread(h.gwas.file)

### FORMAT 
human_ann$category[human_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"
mouse_ann$category[mouse_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"

# mendel <- mendel[,1:7]
colnames(mendel) <- c("chromosome", "start", "dbSNP", "Reference", "Alternate",
                      "ClinicalSignificance", "OriginSimple")
mendel <- mendel[OriginSimple == "germline"]
mendel <- mendel[ClinicalSignificance %like% "Pathogenic" | 
                   ClinicalSignificance %like% "Likely pathogenic" | 
                   ClinicalSignificance %like% "pathogenic"]
# mendel <- mendel[OriginSimple == "germline" & ClinicalSignificance %like% "Pathogenic"]
mendel <- mendel[chromosome == chr]
mendel <- mendel[,c("chromosome", "start")]
mendel$end <- mendel$start
mendel <- unique(mendel)

colnames(gwas) <- c("chromosome", "start", "Parent_term")
gwas <- gwas[chromosome == chr]
gwas <- gwas[,c("chromosome", "start")]
gwas$end <- gwas$start
gwas <- unique(gwas)

# subset by human chromosome
human_ann <- human_ann[chromosome == chr]
align <- align[chromosome_human == chr]

# QC
# human_ann <- human_ann[percentage_Nmask < 1]
# mouse_ann <- mouse_ann[percentage_Nmask < 1]

### RUN FOR ALL ANNNOTATION

alakazam <- function(human_ann, mouse_ann, align, snv_bed){
# for each human annotation in chromosome:
ann <- unique(human_ann$category)
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  ann_sub <- collapse.overlap(human_ann[category == ann[i]][,c("chromosome", "start", "end")])
  snv_ann <- unique(bed.intersect(bed1 = snv_bed, bed2 = ann_sub))
  n_total[i] <- sum(abs( (snv_ann$end + 1) - snv_ann$start))

  if (n_total[i] != 0){
  # orthologous sequences for human annotation to get total alignment
  snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = snv_ann)
  snv_ann_align_h <- unique(snv_ann_align[,c("chromosome_A", "start_A", "end_A")])
  n_aligned[i] <- sum(abs( (snv_ann_align_h$end_A + 1) - snv_ann_align_h$start_A))

  # total of mouse annotation in orthologous sequences for conservation
  snv_ann_conserved <- unique(bed.intersect(bed1 = snv_ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = mouse_ann[category == ann[i]][,c("chromosome", "start", "end")]))
  n_conserved[i] <- sum(abs( (snv_ann_conserved$end + 1) - snv_ann_conserved$start))
  } else {
    n_aligned[i] <- 0
    n_conserved[i] <- 0
  }
  
  print(ann[i])
}

return(data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved))
rm(snv_ann, snv_ann_align, snv_ann_conserved, ann, i, n_aligned, n_conserved, n_total)
}

mendel_out_1 <- alakazam(human_ann, mouse_ann, align, mendel)
gwas_out_1 <- alakazam(human_ann, mouse_ann, align, gwas)


### RUN FOR "Functional - all" "Non-CDS" and "Total"

# snv_bed <- mendel
kadabra <- function(human_ann, mouse_ann, align, snv_bed){
  
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
    
    snv_ann <- unique(bed.intersect(bed1 = snv_bed, bed2 = human_ann_2[,c("chromosome", "start", "end")]))
    n_total[i] <- sum(abs( (snv_ann$end + 1) - snv_ann$start))
    
    if (n_total[i] != 0){
      # orthologous sequences for human annotation to get total alignment
      snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                            "chromosome_mouse", "start_mouse", "end_mouse")],
                                       bed = snv_ann)
      snv_ann_align_h <- unique(snv_ann_align[,c("chromosome_A", "start_A", "end_A")])
      n_aligned[i] <- sum(abs( (snv_ann_align_h$end_A + 1) - snv_ann_align_h$start_A))
      
      # total of mouse annotation in orthologous sequences for conservation
      # snv_ann_conserved <- unique(bed.intersect(bed1 = snv_ann_align[,c("chromosome_B", "start_B", "end_B")],
      #                                           bed2 = mouse_ann_2[,c("chromosome", "start", "end")]))
      # n_conserved[i] <- nrow(snv_ann_conserved)
    } else {
      n_aligned[i] <- 0
      # n_conserved[i] <- 0
    }
    print(ann[i])
  }
  
  return(data.table( chromosome = chr,
                     annotation = ann,
                     n_total = n_total,
                     n_aligned = n_aligned,
                     n_conserved = n_conserved))
}

mendel_out_2 <- kadabra(human_ann, mouse_ann, align, mendel)
gwas_out_2 <- kadabra(human_ann, mouse_ann, align, gwas)

mendel_out <- rbind(mendel_out_1, mendel_out_2)
gwas_out <- rbind(gwas_out_1, gwas_out_2)

### EXPORT

fwrite(mendel_out, out.mendelian.file, append = T)
fwrite(gwas_out, out.gwas.file, append = T)

#####


