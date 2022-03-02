rm(list=ls())
graphics.off()

library(data.table)

# INPUT:
# human annotation bed 
# mouse annotation bed

# OUTPUT:
# human annotation overlap matrix
# mouse annotation overlap matrix

h.ann.file <- "../../data/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
m.ann.file <- "../../data/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
functions.file <- "../functions/FUNCTIONS.R"
out.mouse.file <- "m.mus_GRC38_v101_whole_genome_features_multicell_overlap.matrix"
out.human.file <- "h.sap_GRC38_v101_whole_genome_features_multicell_overlap.matrix"


### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))

### FORMAT
human_ann$category[human_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"
mouse_ann$category[mouse_ann$category %in% c("Exon - 3'UTR", "Exon - 5'UTR" )] <- "Exon - UTR"

### QC 
# human_ann <- human_ann[percentage_Nmask < 1]
# mouse_ann <- mouse_ann[percentage_Nmask < 1]

### FORMAT

# colapse overlap
annotation <- unique(human_ann$category)
human_list <- list()
mouse_list <- list()
for (i in 1:length(annotation)){
  human_tmp <- collapse.overlap(human_ann[category == annotation[i]][,1:3])
  human_tmp$category <- annotation[i]
  human_list[[i]] <- human_tmp
  
  mouse_tmp <- collapse.overlap(mouse_ann[category == annotation[i]][,1:3])
  mouse_tmp$category <- annotation[i]
  mouse_list[[i]] <- mouse_tmp
  print(i)
}
human_squish <- do.call("rbind", human_list)
mouse_squish <- do.call("rbind", mouse_list)
rm(human_tmp, mouse_tmp, human_list, mouse_list, human_ann, mouse_ann)

human_list <- list()
mouse_list <- list()
for (i in 1:length(annotation)){
  
  human_sub1 <- human_squish[category == annotation[i]]
  mouse_sub1 <- mouse_squish[category == annotation[i]]
  
  human_vec <- rep(NA, length(annotation))
  mouse_vec <- rep(NA, length(annotation))
  
  for (j in 1:length(annotation)){
    
    human_sub2 <- human_squish[category == annotation[j]]
    mouse_sub2 <- mouse_squish[category == annotation[j]]
    
    human_intersect <- bed.intersect(human_sub1[,1:3], human_sub2[,1:3])
    mouse_intersect <- bed.intersect(mouse_sub1[,1:3], mouse_sub2[,1:3])
    
    human_vec[j] <- sum(abs((human_intersect$end + 1) - human_intersect$start)) 
    mouse_vec[j] <- sum(abs((mouse_intersect$end + 1) - mouse_intersect$start)) 
    
    print(j)
  }

  human_list[[i]] <- human_vec
  mouse_list[[i]] <- mouse_vec
  
  print(annotation[i])
}

human_overlap <- do.call("rbind", human_list)
mouse_overlap <- do.call("rbind", mouse_list)

colnames(human_overlap) <- annotation
rownames(human_overlap) <- annotation
colnames(mouse_overlap) <- annotation
rownames(mouse_overlap) <- annotation

human_totals <- rep(NA, nrow(human_overlap))
mouse_totals <- rep(NA, nrow(mouse_overlap))
for(i in 1:nrow(human_overlap)){
  human_totals[i] <- human_overlap[i,i]
  mouse_totals[i] <- mouse_overlap[i,i]
}
names(human_totals) <- annotation
names(mouse_totals) <- annotation

# convert to percentage
for (i in 1:length(human_totals)){
  human_overlap[,i] <- (human_overlap[,i] / human_totals[i]) * 100
  mouse_overlap[,i] <- (mouse_overlap[,i] / mouse_totals[i]) * 100
}


### EXPORT
fwrite(human_overlap, out.human.file)
fwrite(mouse_overlap, out.mouse.file)








