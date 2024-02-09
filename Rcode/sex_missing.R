library(readr)

# This R script performs quality checks on genetic data using PLINK sex check and missing genotype rate analysis.
# It checks for sample sex consistency and proportion of missing genotypes.
# Valid samples are identified and saved separately for further analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
  path2 <- args[2]
} else {
  path <- readline(prompt = "Enter the path to the file: ")
  path2 <- readline(prompt = "Enter the path to the file cvs: ")
}

# Load file
plink_sexcheck <- read.table(paste( path, "/plink.sexcheck", sep = "", collapse = NULL) , header=T)
plink_missing <- read.table(paste( path,"/plink.imiss", sep = "", collapse = NULL), header=T)
# Check if the status is not PROBLEM
sex_mismatch <- subset(plink_sexcheck , STATUS=="OK")

# Check if Proportion of missing SNPs is less of 10%
high_missing <- subset(plink_missing , F_MISS < 0.10)
# Save the analysis 
write.table(sex_mismatch[,c("FID", "IID")], paste( path,"/sex.valid", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 
write.table(high_missing [,c("FID", "IID")], paste( path,"/miss.valid", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 