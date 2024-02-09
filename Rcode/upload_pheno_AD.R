library(dplyr)

# R script for merging and updating phenotypes based on "AD","LMCI" and "EMCI" status
# This script merges two input files (ADNIMERGE.csv and .fam file),
# identifies rows with Alzheimer's Disease (AD) status,
# checks for sample IDs in both files, and updates phenotypes accordingly and rename the RID column.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    path <- args[1]
    path2 <- args[2]
    filename <- args[3]
  } else {
    path <- readline(prompt = "Enter the path to the file fam: ")
    path2 <- readline(prompt = "Enter the path to the file cvs: ")
    filename <- readline(prompt = "Enter the file name: ")
  }

# Read the file
cvs <- read.csv(paste(path2, "/ADNIMERGE.csv", sep = "", collapse = NULL) , header=T)
fam <- read.table(paste(path, "/", filename,".fam", sep = "", collapse = NULL) , header=F)
# Merge the file by sampleID
m <- merge(fam, cvs[,1:13], by.x = c("V2"), by.y=c("PTID"), all.x=T)
# Take the rows with AD
name <- row.names(m[ which( m$DX_bl == "AD"| m$DX_bl == "LMCI"| m$DX_bl == "EMCI"  ),])
d <- m[name,]
# Check if the elements in the fam column are present in the second file 
risultato <- fam$V2 %in% d$V2
# Update the phenotypes
fam$V6[risultato] <- 2  # case
fam$V6[!risultato] <- 1  # control
fam$V1 <- as.integer(sub("^.*_.*_(\\d+)$", "\\1", fam$V2))

# Update Sex
sex <- row.names(m[ which( m$PTGENDER == "Female"),])
d <- m[name,]
risultato <- fam$V2 %in% d$V2
fam$V5[risultato] <- 2  # Female
fam$V5[!risultato] <- 1 #Male
write.table(fam, paste(path,"/", filename,".fam", sep = "", collapse = NULL), row.names=F, col.names=F, sep="\t", quote=F) 
