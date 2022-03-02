### SCRIPT that QCs ClinVar SNVs

rm(list = ls())
graphics.off()

library(data.table)

### ARGUMENTS
raw.file <- "variant_summary.txt.gz"
qced.file <- "ClinVar_GRCh38_germline_pathogenic_benign_snps.vcf"

### IMPORT
CV <- fread(paste0("gunzip -cq ", raw.file))

### FORMAT
CV <- CV[Assembly == "GRCh38"] # Subset GRCh38
CV <- CV[Type == "single nucleotide variant"] # Subset SNVs
CV <- CV[Chromosome %in% c(1:22, "X")] # Subset chromosomes
CV <- CV[ClinicalSignificance %like% "Benign" | # Subset benign and pathogenic variants
         ClinicalSignificance %like% "Likely benign" |
         ClinicalSignificance %like% "Pathogenic" | 
         ClinicalSignificance %like% "Likely pathogenic"]
CV <- CV[ReviewStatus == "criteria provided, multiple submitters, no conflicts" | # Subset by review status 
         ReviewStatus == "criteria provided, single submitter" |
         ReviewStatus == "reviewed by expert panel"]
CV <- CV[OriginSimple %like% "germline"] # Subset by origin
output <- CV[,c("Chromosome", "PositionVCF", "RS# (dbSNP)", "ReferenceAlleleVCF", "AlternateAlleleVCF",
            "ClinicalSignificance", "OriginSimple", "Type", "GeneSymbol")]
# remove

### EXPORT
fwrite(output, qced.file, sep = "\t", col.names = F)






