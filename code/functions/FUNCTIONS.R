### SCRIPT that contains functions used for the workflow

### FUNCTIONS included:
#   orthologous.seq()
#   orthologous.pos() !!!
#   collapse.overlap()
#   bed.intersect()
#   bed.to.long()
#   long.to.bed()
#   vec.to.bed()
#   bed.match.count() !!!


### FUNCTION that returns the orthologous sequence given an alignment and bed file. 
## INPUT:
# alignment - file with columns: "chromosome_A", "start_A", "end_A", "chromosome_B", "start_B", "end_B"
# bed - file with columns: "chromosome_A", "start_A", "end_A"
## OUTPUT:
# file with columns: "chromosome_A", "start_A", "end_A", "chromosome_B", "start_B", "end_B"

### change output to: "chromosome_A", "start_A", "end_A", 
#                     "chromosome_A_match", "start_A_match", "end_A_match",
#                     "chromosome_B_match", "start_B_match", "end_B_match"

orthologous.seq <- function(alignment, bed){
  
  require(data.table)
  
  # format
  alignment <- as.data.table(alignment)
  bed <- as.data.table(bed)
  colnames(alignment) <- c("chromosome_A", "start_A", "end_A", "chromosome_B", "start_B", "end_B")
  colnames(bed) <- c("chromosome_A", "start_A", "end_A")
  alignment[, c(1,4)] <- lapply(alignment[, c(1,4)], as.character)
  alignment[, c(2:3,5:6)] <- lapply(alignment[, c(2:3,5:6)], as.integer)
  bed[, 1] <- lapply(bed[, 1], as.character)
  bed[, c(2:3)] <- lapply(bed[, c(2:3)], as.integer)
  
  # alignment to long
  align_long <- alignment[, .(chromosome_A = chromosome_A, 
                              pos_A = start_A:end_A, 
                              chromosome_B = chromosome_B, 
                              pos_B = start_B:end_B), 
                          by=.(1:nrow(alignment))][, nrow := NULL] # bed file to vector
  
  # bed to long
  bed_long_list <- list() # for loop to avoid Error: negative length vectors are not allowed
  for (j in 1:nrow(bed)){
    bed_sub <- bed[j,]
    bed_long_list[[j]] <- bed_sub[, .(chromosome_A = chromosome_A, 
                                      pos_A = start_A:end_A), 
                                  by=.(1:nrow(bed_sub))][, nrow := NULL] # bed file to vector
  }
  bed_long <- do.call("rbind", bed_long_list)
  
  # subset alignmentlong by bed long
  long_sub <- align_long[bed_long, on = c("chromosome_A", "pos_A")]
  long_sub <- long_sub[complete.cases(long_sub),]
  long_sub <- unique(long_sub)
  
  long_sub[,`:=`(pos1id=rleid(cumsum(abs(diff(c(0,pos_A)))>1)), # add grouping variable
                 pos2id=rleid(cumsum(abs(diff(c(0,pos_B)))>1)),
                 cat1id=rleid(chromosome_A),
                 cat2id=rleid(chromosome_B))][
                   , `:=`(grp=.GRP), by = c("pos1id","pos2id",
                                            "cat1id","cat2id")         
                   ][, `:=`(pos1id = NULL, pos2id = NULL, cat1id = NULL, cat2id = NULL)]
  
  # tmp <- long_sub[pos_A > 45908064 & pos_A < 45908079]
  
  # long to short
  output <- long_sub[order(chromosome_A, pos_A), # order by chromosome and pos
                     .(start_A = min(pos_A),
                       end_A = max(pos_A),
                       start_B = pos_B[1],
                       end_B = pos_B[length(pos_B)]),
                     by = .(grp = grp,
                            chromosome_A = chromosome_A,
                            chromosome_B = chromosome_B) # split by grp
                     ][, grp := NULL]
  
  # if unequal lengths, return warning
  ind <- which(abs(output$end_A - output$start_A) != abs(output$end_B - output$start_B))
  if (length(ind) != 0){ print("TROUBLESHOOT: unequal lengths returned in alignment!") }
  
  output <- output[,c("chromosome_A", "start_A", "end_A", "chromosome_B", "start_B", "end_B")]
  output <- output[complete.cases(output)]
  
  return(output)
}

### FUNCTION that collapses overlapping sequences in a bed file
## INPUT:
# bed - bed file with columns: "chromosome", "start", "end"
## OUTPUT:
# bed file with columns: "chromosome", "start", "end"
collapse.overlap <- function(bed) {
  
  require(GenomicRanges)
  require(data.table)
  
  bed <- as.data.table(bed)
  colnames(bed) <- c("CHROM","START","STOP")
  bed$CHROM <- as.character(bed$CHROM)
  bed[, c(2:3)] <- lapply(bed[, c(2:3)], as.integer)
  setkey(bed)
  bedKey <- c("CHROM","START","STOP")
  setnames(bed,colnames(bed),bedKey)
  
  if(!identical(key(bed),bedKey)) setkeyv(bed,bedKey)
  grBed <- makeGRangesFromDataFrame(bed,
                                    seqnames.field = "CHROM",start.field="START",end.field="STOP")
  grBed <- reduce(grBed)
  grBed <- data.table(
    CHROM=as.character(seqnames(grBed)),
    START=start(grBed),
    STOP=end(grBed),
    key = c("CHROM","START","STOP"))
  
  colnames(grBed) <- c("chromosome", "start", "end")
  
  return(grBed)
}

### FUNCTION that returns intersecting sequennces between two bed files
## INPUT:
# bed1 -- bed file with columns: "chromosome", "start", "end"
# bed2 -- bed file with columns: "chromosome", "start", "end"
## OUTPUT:
# bed file with columns: "chromosome", "start", "end"
bed.intersect <- function(bed1, bed2) {
  
  require(data.table)
  
  rowShift <- function(x, shiftLen = 1L) {
    #Note this function was described in this thread:
    #http://stackoverflow.com/questions/14689424/use-a-value-from-the-previous-row-in-an-r-data-table-calculation
    r <- (1L + shiftLen):(length(x) + shiftLen)
    r[r<1] <- NA
    return(x[r])
  }
  
  collapse.overlap <- function(bed) {
    
    require(GenomicRanges)
    require(data.table)
    
    bed <- as.data.table(bed)
    colnames(bed) <- c("CHROM","START","STOP")
    bed$CHROM <- as.character(bed$CHROM)
    bed[, c(2:3)] <- lapply(bed[, c(2:3)], as.integer)
    setkey(bed)
    bedKey <- c("CHROM","START","STOP")
    setnames(bed,colnames(bed),bedKey)
    
    if(!identical(key(bed),bedKey)) setkeyv(bed,bedKey)
    grBed <- makeGRangesFromDataFrame(bed,
                                      seqnames.field = "CHROM",start.field="START",end.field="STOP")
    grBed <- reduce(grBed)
    grBed <- data.table(
      CHROM=as.character(seqnames(grBed)),
      START=start(grBed),
      STOP=end(grBed),
      key = c("CHROM","START","STOP"))
    
    colnames(grBed) <- c("chromosome", "start", "end")
    
    return(grBed)
  }
  
  colnames(bed1) <- c("CHROM","START","STOP")
  colnames(bed2) <- c("CHROM","START","STOP")
  bed1$CHROM <- as.character(bed1$CHROM)
  bed2$CHROM <- as.character(bed2$CHROM)
  bed1[, c(2:3)] <- lapply(bed1[, c(2:3)], as.integer)
  bed2[, c(2:3)] <- lapply(bed2[, c(2:3)], as.integer)
  bed1[ START > STOP, `:=`( START = STOP, STOP = START)]   # ensure start <= end
  bed2[ START > STOP, `:=`( START = STOP, STOP = START)]   # ensure start <= end
  setkey(bed1)
  setkey(bed2)
  bedKey <- c("CHROM","START","STOP")
  
  
  if(nrow(bed1)>nrow(bed2)) {
    bed <- foverlaps(bed1, bed2, nomatch = 0)
  } else {
    bed <- foverlaps(bed2, bed1, nomatch = 0)
  }
  bed[, START := pmax(START, i.START)]
  bed[, STOP := pmin(STOP, i.STOP)]
  bed[, `:=`(i.START = NULL, i.STOP = NULL)]
  if(!identical(key(bed),bedKey)) setkeyv(bed,bedKey)
  if(any(bed[, STOP+1 >= rowShift(START), by=CHROM][,V1], na.rm = T)) {
    bed <- collapse.overlap(bed)
  }
  
  colnames(bed) <- c("chromosome", "start", "end")
  return(bed)
}


### FUNCTION that converts a bed file to a data.table with all pos (short to long)
### INPUT
# bed - bed file with columns: chromosome, start, end
### OUTPUT
# data.table with columns: chromosome, pos)
bed.to.long <- function(bed){
  require(data.table)
  bed <- as.data.table(bed)
  colnames(bed) <- c("chromosome", "start", "end")
  output <- bed[, .(chromosome = chromosome, pos = start:end), by=.(1:nrow(bed))][, nrow := NULL] # bed file to vector
  return(output)
}

### FUNCTION that converts a long file to a bed file (long to short)
### INPUT
# long - long file with columns: chromosome, pos
### OUTPUT
# bed file with columns: chromosome, start, end)
long.to.bed <- function(long){
  require(data.table)
  long <- as.data.table(long)
  colnames(long) <- c("chromosome", "pos")
  output <- long[order(pos), 
                   .(start = min(pos),
                     end = max(pos)),
                   by = .(grp = rleid(c(0, cumsum(diff(pos) > 1))), 
                          chromosome = chromosome)
                   ][, grp := NULL]
  return(output)
}

### FUNCTION that converts vector of integers to DT of start and end coordinates
vec.to.bed <- function(vec){
  require(data.table)
  vec_dt <- data.table(V1 = vec)
  output <- vec_dt[order(V1),
                   .(chromosome = NA,
                     start = min(V1),
                     end = max(V1)),
                   by = .(grp = rleid(c(0, cumsum(diff(V1) > 1))))
                   ][, grp := NULL]
  return(output)
}

# ### FUNCTION that returns the seqence match between bed files
# # (ie for each sequence in bed file1, how many integers overlap with bed file2?)
# bed.match.count <- function(bed1, bed2){
# 
#   require(data.table)
# 
#   bed1 <- setDT(bed1) # set as bed as data.table
#   bed2 <- setDT(bed2) # set as bed as data.table
#   
#   colnames(bed1) <- c("chromosome", "start", "end")
#   colnames(bed2) <- c("chromosome", "start", "end")
#   
#   #####
#   
#   # bed_dt[, ind := .I] # add uniqe index to data.table
#   # setkey(bed_dt) # sets keys // order data by all columns
#   # 
#   # vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
#   # vec_dt[, vec2 := vec] # dublicate column
#   # 
#   # # Fast overlap join:
#   # ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
#   #                  type = "within", nomatch = 0L)
#   # counts <- ans1[, .N, keyby = ind] # count by ind
#   # # merge to inital data
#   # bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
#   # setorder(bed_dt, ind) # reorder by ind to get inital order
#   # bed_dt[, ind := NULL] # deletes ind colum
#   # bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
#   
#   #####
# 
#   return(bed_dt$n_match)
# }

### FUNCTION that randomply samples positions from a bed file
# bed = bed file input
# nPOS = number of coordinates to sample
bed.sample.pos <- function(bed, nPOS){
  require(data.table)
  
  colnames(bed) <- c("chromosome", "start", "end")
  # calculate size of range
  bed[, size := 1 + end-start]
  # Randomly sample bed file rows, proportional to the length of each range
  simulated.sites <- bed[sample(.N, size=nPOS, replace=TRUE, prob=bed$size)]
  # Randomly sample uniformly within each chosen range !!! taking too long as start:end vectory too large
  # simulated.sites[, position := sample(start:end, size=1), by=1:nrow(simulated.sites)] # old 
  simulated.sites$position <- mapply(function(x, y) sample(seq(x, y), 1), simulated.sites$start, simulated.sites$end)
  # Remove extra columns and format as needed
  simulated.sites[, start  := position]
  simulated.sites[, end := position]
  simulated.sites[, c("size", "position") := NULL]
  
  return(simulated.sites)
}

### FUNCTION that  provides coordinates for sequences that contain all the integers 
# in bed1 that do not overlap with any integers from the sequences in bed2. 
bed_antijoin <- function(bed1, bed2){
  
  require(data.table)

  colnames(bed1) <- c("chromosome", "start", "end")
  colnames(bed2) <- c("chromosome", "start", "end")
  
  # setkey(bed2, chromosome, start, end)
  # ds <- foverlaps(bed1, bed2,  type="any")
  # ds <- ds[,.(chromosome, 
  #       start = fcase(is.na(start) | i.start <= start, i.start,
  #                     i.end >= end, end + 1),
  #       end = fcase(is.na(end) | i.end >= end, i.end,
  #                   i.start <= start, start - 1)
  # )]
  # return(ds)
  
  # Alternative method
  require(GenomicRanges)

  ds2 <- setdiff(makeGRangesFromDataFrame(bed1), makeGRangesFromDataFrame(bed2))
  ds2 <- as.data.table(ds2)
  ds2 <- ds2[,c("seqnames", "start", "end")]
  colnames(ds2) <- c("chromosome", "start", "end")
  return(ds2)
}




