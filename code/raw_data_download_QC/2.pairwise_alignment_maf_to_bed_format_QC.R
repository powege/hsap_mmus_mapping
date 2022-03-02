# SCRIPT that converts Multiple Allignment Format (.MAF) to bed file of orthologous bases
# excludes indels and scaffolds 
# includes autosomes and X chromosome.

rm(list = ls())
graphics.off()

library(stringr)
library(stringi)
library(plyr)
library(dplyr)
library(utils)
library(data.table)
library(zoo)
library(reshape2)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

### set args variables
input.path <- args[1]
output.path <- args[2]
dir <- args[3]
in.prefix <- args[4]
out.prefix <- args[5]
chr.id <- as.integer(args[6])


### FUNCTIONS

# Function that returns sequence chromosome as string from .maf input
turnip_CHR <- function(ref1){
  
  vars <- sub("^[^.]*\\.", "", ref1)
  vars <- strsplit(vars, "\\s+")
  chr <- vars[[1]][1]
  refseq <- vars[[1]][6]
  refseq <- strsplit(refseq, "")[[1]]
  chr_vec <- rep(NA, length(refseq))
  chr_vec[refseq != "-"] <- chr
  
  return(chr_vec)
}


# Function that retuns sequence coordinates as string from .maf input
turnip_POS <- function(ref1){
  
  vars <- sub("^[^.]*\\.", "", ref1)
  vars <- strsplit(vars, "\\s+")
  
  # identify strand direction
  strand <- vars[[1]][4]
  
  # if + strand:
  if (strand == "+"){
    start <- as.integer(vars[[1]][2]) + 1
    refseq <- vars[[1]][6]
    refseq <- strsplit(refseq, "")[[1]]
    refbp <- rep(NA_integer_, length(refseq))
    refbp[refseq != "-"] <- seq(from = start, by = 1, length.out = sum(refseq != "-"))
  }
  
  # if - strand: 
  if (strand == "-"){
    chr_length <- as.integer(vars[[1]][5])
    minus <- as.integer(vars[[1]][2])
    start <- chr_length - minus
    refseq <- vars[[1]][6]
    refseq <- strsplit(refseq, "")[[1]]
    refbp <- rep(NA_integer_, length(refseq))
    refbp[refseq != "-"] <- seq(from = start, by = -1, length.out = sum(refseq != "-"))
  }
  
  return(refbp)
}

# Function that retuns sequence bases as string from .maf input
turnip_REF <- function(ref1){
  
  refseq <- str_extract(ref1, "[\\w\\-]+(?=\\s*)$")
  refseq <- strsplit(refseq, "")[[1]]
  refseq <- replace(refseq, refseq == "-", NA)
  
  return(refseq)
}


# Function that counts the nember of .maf files for a given chromosome
file_count <- function(input.path, dir, in.prefix, chr){
  files <- length(list.files(path = paste0(input.path, dir), pattern = paste0(in.prefix, chr, "_")))
  return(files)
}

# Function that converts long to short by ind
# sub <- subset(df_long, df_long$ID == "1111")
long_to_short <- function(sub){
  data.table(
    CHR_A = sub$spec1.chr[1],
    POS_A_START = sub$spec1.pos[1],
    POS_A_END = sub$spec1.pos[nrow(sub)],
    CHR_B = sub$spec2.chr[1],
    POS_B_START = sub$spec2.pos[1],
    POS_B_END = sub$spec2.pos[nrow(sub)]
  )
}

# Function that converts .maf files to to tab delimited file with each base and genomic coordinates
# output cols: (POS; REF; SPECIES)
MAF_convert <- function(chr, files.chr, input.path, output.path, dir, in.prefix, out.prefix){
  
  # j=30
  for (j in 1:files.chr){
    
    # read
    x <- read.table(paste0(input.path, dir, in.prefix, chr, "_", j, ".maf"), sep = "\t")
    
    # format
    x <- as.character(x$V1)
    
    ind.spec1 <- seq(3, length(x), by=4)
    spec1 <- x[ind.spec1]
    ind.spec2 <- seq(4, length(x), by=4)
    spec2 <- x[ind.spec2]
    
    spec1.pos <- unlist(lapply(spec1, turnip_POS))
    spec2.pos <- unlist(lapply(spec2, turnip_POS))
    
    spec1.chr <- unlist(lapply(spec1, turnip_CHR))
    spec2.chr <- unlist(lapply(spec2, turnip_CHR))
    
    spec1.ref <- unlist(lapply(spec1, turnip_REF))
    spec2.ref <- unlist(lapply(spec2, turnip_REF))
    
    spec1.ref <- toupper(spec1.ref) # remove softmask
    spec2.ref <- toupper(spec2.ref)
    bases <- c("A", "T", "C", "G")
    spec1.ref[!spec1.ref %in% bases] <- NA
    spec2.ref[!spec2.ref %in% bases] <- NA
    
    df_long <- data.frame(spec1.chr, spec1.pos, spec1.ref, spec2.chr, spec2.pos, spec2.ref)
    df_long <- subset(df_long, df_long$spec2.chr %in% c(1:19, "X")) # subset automomes and X alignment (remove patches)
    df_long <- df_long[complete.cases(df_long),] # filter indels
    
    # convert long to short format
    INT_A_split <- cumsum(c(TRUE, abs(diff(df_long$spec1.pos))!=1))
    INT_B_split <- cumsum(c(TRUE, abs(diff(df_long$spec2.pos))!=1))
    CAT1_A_split <- cumsum(c(1, diff(as.numeric(as.factor(df_long$spec1.chr))) != 0))
    CAT1_B_split <- cumsum(c(1, diff(as.numeric(as.factor(df_long$spec2.chr))) != 0))
    df_long_ind <- data.table(INT_A_split,
                              INT_B_split,
                              CAT1_A_split,
                              CAT1_B_split)
    df_long$ID <- apply(df_long_ind, 1, paste0, collapse = "") # unique ID for changes
    df_long$ID <- as.factor(df_long$ID)
    dt_short <- ddply(df_long, "ID", long_to_short)
    df_long$ID <- NULL
    dt_short$ID <- NULL
    
    df_long <- df_long[,c("spec1.pos", "spec1.ref", "spec2.ref")]
    colnames(dt_short) <- c("chromosome_A", "start_A", "end_A",
                            "chromosome_B", "start_B", "end_B")
    # write tables
    fwrite(dt_short, paste0(output.path, out.prefix, chr, ".bed.gz"), append = T, compress = "gzip", col.names = F)
    fwrite(df_long, paste0(output.path, out.prefix, chr, ".long.gz"), append = T, compress = "gzip", col.names = F)
    
    print(j)
  }
}


### SCRIPT

# set chr
chr <- c(1:22, "X")[chr.id]

# count files
files.chr <- file_count(input.path = input.path, 
                        dir = dir, 
                        chr = chr, 
                        in.prefix = in.prefix)

### Convert
# chr = chr
# files.chr = files.chr
# input.path = input.path
# output.path = output.path
# dir = dir
# in.prefix = in.prefix
# out.prefix = out.prefix
MAF_convert(chr = chr,
            files.chr = files.chr,
            input.path = input.path,
            output.path = output.path,
            dir = dir,
            in.prefix = in.prefix,
            out.prefix = out.prefix)


##########
