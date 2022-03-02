### SCRIPT that formats Ensembl gene and Regulatory Build annotation, and TAD boundries
### Outputs all seqences regardless of overlap
### inferes Intron and Unannotated regions by chromosome
### Outputs cols: chromosome; start; end; category; strand; id

rm(list=ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
reg.build.heart.gff <- args[1] # .gff file
reg.build.kidney.gff <- args[2] # .gff file
reg.build.spleen.gff <- args[3] # .gff file
gencode.gtf <- args[4] # .gtf file
tad.file <- args[5]
chr.nbp.file <- args[6]
N.mask.file <- args[7]
sm.file <- args[8]
out.file <- args[9] 
chr <- as.character(args[10])

### FUNCTIONS

### FUNCTION that returns the seqence match between bed file and vector
# (ie for each sequence in bed file1, how many integers overlap with vector?)
bed_match <- function(bed, vec){
  
  require(data.table)
  
  bed_dt <- setDT(bed) # set as bed as data.table
  colnames(bed_dt) <- c("start", "end")
  bed_dt[, ind := .I] # add uniqe index to data.table
  setkey(bed_dt) # sets keys // order data by all columns
  
  vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
  vec_dt[, vec2 := vec] # dublicate column
  
  # Fast overlap join:
  ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
  setorder(bed_dt, ind) # reorder by ind to get inital order
  bed_dt[, ind := NULL] # deletes ind colum
  bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
  
  return(bed_dt$n_match) 
}

# sub <- subset(exon_all, exon_all$transcript_id == "ENST00000623180")
### FUNCTION that infers intron - distil POS from exon POS
get_intron_proximal_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome,
                      start = NA,
                      end = NA,
                      transcript_id = NA)
  }
  if (nrow(sub) > 1){
    sub <- sub[order(sub$end),] # order by exon end POS!!!
    start <- c((sub$end[1:(nrow(sub))-1] + 1), (sub$start[2:(nrow(sub))] - 10))
    end <- c((sub$end[1:(nrow(sub))-1] + 10), (sub$start[2:(nrow(sub))] - 1))
    chromosome <- rep(sub$chromosome[1], length(start))
    transcript_id <- rep(sub$transcript_id[1], length(start))
    out <- data.frame(chromosome,
                      start,
                      end,
                      transcript_id)
  }
  return(out)
}

# sub <- subset(exon_all, exon_all$transcript_id == "ENST00000450305")
### FUNCTION that infers intron POS from exon POS
get_intron_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome,
                      start = NA,
                      end = NA)
  }
  if (nrow(sub) > 1){
    # sub <- sub[order(sub$exon_number),] # order by exon number
    sub <- sub[order(sub$end),] # order by exon end POS!!!
    start <- sub$end[1:(nrow(sub))-1] + 1
    end <- sub$start[2:(nrow(sub))] - 1
    chrom <- rep(sub$chromosome[1], length(start))
    out <- data.frame(chromosome = chrom,
                      start = start,
                      end = end)
  }
  # dif <- sub$exon_chrom_end - sub$exon_chrom_start
  # len <- sum(dif) + length(dif)
  # out <- c(sub$ensembl_transcript_id[1], len)
  # names(out) <- c("ensembl_transcript_id", "intron_length")
  return(out)
}

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT 
  
reg_build_heart <-  fread(paste0("gunzip -cq ", reg.build.heart.gff))
reg_build_kidney <-  fread(paste0("gunzip -cq ", reg.build.kidney.gff))
reg_build_spleen <-  fread(paste0("gunzip -cq ", reg.build.spleen.gff))
gene <- fread(paste0("gunzip -cq ", gencode.gtf))
tad_pos <- fread(tad.file)
chr_nbp <- fread(chr.nbp.file, select = 1:3)
N_mask <- fread(paste0("gunzip -cq ", N.mask.file))
soft_mask <- fread(paste0("gunzip -cq ", sm.file))

### FORMAT 

# colnames
colnames(N_mask) <- c("chromosome", "start", "end")
colnames(soft_mask) <- c("chromosome", "start", "end")
colnames(chr_nbp) <- c("chromosome", "start", "end")
colnames(gene) <- c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
colnames(reg_build_heart) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
colnames(reg_build_kidney) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
colnames(reg_build_spleen) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

chr_length <- chr_nbp$end[chr_nbp$chromosome == chr] # chromosome length

### FORMAT gene annotation
# pull out transcript_biotype and transcript_id from attribute
# gene$gene_biotype <- str_match(string=gene$attribute, pattern="gene_biotype\\s+\"([^\"]+)\"")[,2]
gene$transcript_biotype <- str_match(string=gene$attribute, pattern="transcript_biotype\\s+\"([^\"]+)\"")[,2]
gene$transcript_id <- str_match(string=gene$attribute, pattern="transcript_id\\s+\"([^\"]+)\"")[,2]
gene$exon_number <- str_match(string=gene$attribute, pattern="exon_number\\s+\"([^\"]+)\"")[,2]
# subset CDS and UTR
cds <- subset(gene, gene$type == "CDS")
utr3 <- subset(gene, gene$type == "three_prime_utr")
utr5 <- subset(gene, gene$type == "five_prime_utr")
# identify transcript_id in CDS and UTR 
PCtranscript <- unique(c(cds$transcript_id, utr5$transcript_id, utr3$transcript_id))
# subset all exons
exon_all <- subset(gene, gene$type == "exon")
# subset all exon with non-PC transcript
exon_nonPC <- exon_all[which(!exon_all$transcript_id %in% PCtranscript),]

### FORMAT reg build

reg_build <- separate(reg_build_heart, attribute, c("activity_heart", "bound_end","bound_start", "description", "epigenome", "feature_type", "regulatory_feature_stable_id"), sep = ";")
reg_build$regulatory_feature_stable_id <- gsub("regulatory_feature_stable_id=", "", reg_build$regulatory_feature_stable_id)
reg_build$activity_heart <- gsub("activity=", "", reg_build$activity_heart)
reg_build_kidney <- separate(reg_build_kidney, attribute, c("activity_kidney", "bound_end","bound_start", "description", "epigenome", "feature_type", "regulatory_feature_stable_id"), sep = ";")
reg_build$activity_kidney <- gsub("activity=", "", reg_build_kidney$activity_kidney)
reg_build_spleen <- separate(reg_build_spleen, attribute, c("activity_spleen", "bound_end","bound_start", "description", "epigenome", "feature_type", "regulatory_feature_stable_id"), sep = ";")
reg_build$activity_spleen <- gsub("activity=", "", reg_build_spleen$activity_spleen)

### FORMAT annotated categories

# subset promotor POS
promoter_pos <- subset(reg_build, reg_build$feature == "promoter")
promoter_pos$category <- "Promoter"
promoter_pos <- promoter_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
promoter_pos <- unique(promoter_pos)
# subset promotor flanking POS
promoter_flank_pos <- subset(reg_build, reg_build$feature == "promoter_flanking_region")
promoter_flank_pos$category <- "Enhancer - proximal"
promoter_flank_pos <- promoter_flank_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
promoter_flank_pos <- unique(promoter_flank_pos)
# subset enhancer POS
enhancer_pos <- subset(reg_build, reg_build$feature == "enhancer")
enhancer_pos$category <- "Enhancer - distal"
enhancer_pos <- enhancer_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
enhancer_pos <- unique(enhancer_pos)
# subset CTCF binding site POS
CTCF_binding_pos <- subset(reg_build, reg_build$feature == "CTCF_binding_site")
CTCF_binding_pos$category <- "CTCF binding"
CTCF_binding_pos <- CTCF_binding_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
CTCF_binding_pos <- unique(CTCF_binding_pos)
# subset TF binding site POS
TF_binding_pos <- subset(reg_build, reg_build$feature == "TF_binding_site")
TF_binding_pos$category <- "TF binding"
TF_binding_pos <- TF_binding_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
TF_binding_pos <- unique(TF_binding_pos)
# subset open chromatin POS
open_chromatin_pos <- subset(reg_build, reg_build$feature == "open_chromatin_region")
open_chromatin_pos$category <- "Open chromatin"
open_chromatin_pos <- open_chromatin_pos[,c("chromosome", "start", "end", "category", "strand", "activity_heart", "activity_kidney", "activity_spleen", "regulatory_feature_stable_id")]
open_chromatin_pos <- unique(open_chromatin_pos)
# subset exon CDS POS
cds$category <- "Exon - CDS"
exon_CDS_pos <- cds[,c("chromosome", "start", "end", "category", "strand", "transcript_id", "exon_number")]
exon_CDS_pos <- unique(exon_CDS_pos)
# subset exon UTR 5' POS
utr5$category <- "Exon - 5'UTR"
exon_UTR5_pos <- utr5[,c("chromosome", "start", "end", "category", "strand", "transcript_id", "exon_number")]
exon_UTR5_pos <- unique(exon_UTR5_pos)
# subset exon UTR 3' POS
utr3$category <- "Exon - 3'UTR"
exon_UTR3_pos <- utr3[,c("chromosome", "start", "end", "category", "strand", "transcript_id", "exon_number")]
exon_UTR3_pos <- unique(exon_UTR3_pos)
# subset exon non-coding POS
exon_nonPC$category <- "Exon - other"
exon_NC_pos <- exon_nonPC[,c("chromosome", "start", "end", "category", "strand", "transcript_id", "exon_number")]
exon_NC_pos <- unique(exon_NC_pos)
# subset TAD boundries
tad_pos$category <- "TAD boundry"
tad_pos$strand <- "."
tad_pos <- unique(tad_pos)

# rbind
annotation <- rbind.fill(exon_CDS_pos,
                         exon_UTR5_pos,
                         exon_UTR3_pos,
                         exon_NC_pos,
                         promoter_pos,
                         promoter_flank_pos,
                         enhancer_pos,
                         CTCF_binding_pos,
                         TF_binding_pos,
                         open_chromatin_pos,
                         tad_pos)
# remove duplicates
annotation <- unique(annotation)
# colSums(is.na(annotation))

### FORMAT introns 

# subset all exons
exon_all <- subset(exon_all, exon_all$chromosome == chr)
# identify intron - proximal from exons
intron_prox_pos <- ddply(exon_all, "transcript_id", get_intron_proximal_POS)
intron_prox_pos <- intron_prox_pos[complete.cases(intron_prox_pos),]
intron_prox_pos$category <- "Intron - proximal"
intron_prox_pos$strand <- "."
intron_prox_pos <- unique(intron_prox_pos)
# subset annotation by chr
annotation <- subset(annotation, annotation$chromosome == chr)
# rbind 
annotation <- rbind.fill(annotation, intron_prox_pos)
# identify intron - distil from exons
intron_dis_all <- ddply(exon_all, "transcript_id", get_intron_POS)
intron_dis_all <- intron_dis_all[complete.cases(intron_dis_all),]
# exclude regions already annotated
tmp_ann <- annotation[,c("start", "end")]
colnames(tmp_ann) <- c("from", "to")
tmp_int <- intron_dis_all[,c("start", "end")]
colnames(tmp_int) <- c("from", "to")
# vector of intron POS not in other annotations
intron_pos <- setdiff(x = unlist(seq2(from = tmp_int$from, to = tmp_int$to)), 
                      y = unlist(seq2(from = tmp_ann$from, to = tmp_ann$to)))
# intron_pos <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
# format vector to dataframe
intron_pos <- sort(intron_pos, decreasing = F)
intron_pos <- t(sapply(split(intron_pos, findInterval(intron_pos, intron_pos[which(c(1, diff(intron_pos)) > 1)])), range))
intron_pos <- as.data.table(intron_pos)
colnames(intron_pos) <- c("start", "end")
intron_pos$chromosome <- chr
intron_pos$category <- "Intron - distal"
intron_pos <- unique(intron_pos)
intron_pos$strand <- "."
# rbind 
annotation <- rbind.fill(annotation, intron_pos)


### FORMAT unannotated

tmp_ann <- annotation[,c("start", "end")]
colnames(tmp_ann) <- c("from", "to")
# vector of intron POS not in other annotations
unan_pos <- setdiff(x = 1:chr_length, 
                    y = unlist(seq2(from = tmp_ann$from, to = tmp_ann$to)))
# unan_pos <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
# format vector to dataframe
unan_pos <- sort(unan_pos, decreasing = F)
unan_pos <- t(sapply(split(unan_pos, findInterval(unan_pos, unan_pos[which(c(1, diff(unan_pos)) > 1)])), range))
unan_pos <- as.data.table(unan_pos)
colnames(unan_pos) <- c("start", "end")
unan_pos$chromosome <- chr
unan_pos$category <- "Unannotated"
unan_pos <- unique(unan_pos)
unan_pos$strand <- "."
# rbind 
annotation <- rbind.fill(annotation, unan_pos)

# complete cases
annotation <- annotation[!is.na("chromosome") | 
                   !is.na("start") |
                   !is.na("end") |
                   !is.na("category")]
# ensure start >= end
annotation <- as.data.table(annotation)
annotation[ start > end, `:=`( start = end, end = start)]
# remove duplicates 
annotation <- unique(annotation) 
# convert "Open chromation" and "TF binding" to "Miscellaneous"
annotation$category <- as.character(annotation$category)
annotation$category[annotation$category == "Open chromatin" | annotation$category == "TF binding"] <- "Miscellaneous"
# add unique ID 
annotation$category_code <- mapvalues(annotation$category, from=c("Exon - CDS",
                                                                  "Exon - 5'UTR",
                                                                  "Exon - 3'UTR",
                                                                  "Exon - other",
                                                                  "Promoter",
                                                                  "Enhancer - proximal",
                                                                  "Enhancer - distal",
                                                                  "CTCF binding",
                                                                  "Miscellaneous",
                                                                  "Intron - proximal",
                                                                  "Intron - distal",
                                                                  "Unannotated",
                                                                  "TAD boundry"), # number code annotations (plyr)
                               to=c("A","B","C","D","E","F","G","H","I","J","K","L","M"))
tmp_list <- list()
for (cat in 1:length(unique(annotation$category))){
  tmp_sub <- subset(annotation, annotation$category == unique(annotation$category)[cat])
  tmp_sub$category_id <- paste(tmp_sub$category_code, tmp_sub$chromosome, 1:length(tmp_sub$category), sep = ".")
  tmp_list[[cat]] <- tmp_sub
}
annotation <- do.call("rbind", tmp_list)
annotation$category_code <- NULL

### CALCULATE FRACTION OF N MASKED BASES

annotation$n_Nmask <- bed_match(bed = annotation[,c("start", "end")],
                                vec = unlist(seq2(from = N_mask$start[N_mask$chromosome == chr], to = N_mask$end[N_mask$chromosome == chr])))
annotation$length <- (annotation$end +1) - annotation$start
annotation$percentage_Nmask <- (annotation$n_Nmask / annotation$length) * 100 
# hist(annotation$percentage_Nmask)

### CALCULATE FRACTION OF SOFT MASKED BASES

annotation$n_sm <- bed_match(bed = annotation[,c("start", "end")],
                                vec = unlist(seq2(from = soft_mask$start[soft_mask$chromosome == chr], to = soft_mask$end[soft_mask$chromosome == chr])))
annotation$length <- (annotation$end +1) - annotation$start
annotation$percentage_sm <- (annotation$n_sm / annotation$length) * 100 
# hist(annotation$percentage_sm)


annotation <- annotation[,c( "chromosome", 
                             "start", 
                             "end", 
                             "category", 
                             "category_id",                 
                             "regulatory_feature_stable_id",
                             "transcript_id",             
                             "exon_number",                 
                             "activity_heart",           
                             "activity_kidney",     
                             "activity_spleen",             
                             "percentage_Nmask",
                             "percentage_sm",
                             "strand")]

### EXPORT 
fwrite(annotation, out.file, compress = "gzip", append = TRUE)


### RESOURCES 
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


#####
