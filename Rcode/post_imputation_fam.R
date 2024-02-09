library(dplyr)

# This R script reads multiple input files including ADNIMERGE.csv, .fam file, and TEMP6.fam.
# It updates the SEX column in the .fam file based on the information in TEMP6.fam.
# Additionally, it modifies the names of RID and FID columns in the .fam file.
# Samples that are labeled as "Not Hisp/Latino" in ADNIMERGE.csv are removed from the analysis.
# Phenotypes are updated based on specific criteria from ADNIMERGE.csv, and the .fam file is accordingly modified.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    path <- args[1]
    path2 <- args[2]
    filename <- args[3]
    path3<-args[4]
  } else {
    path <- readline(prompt = "Enter the path to the file fam: ")
    path2 <- readline(prompt = "Enter the path to the file cvs: ")
    filename <- readline(prompt = "Enter the file name: ")
    path3 <- readline(prompt = "Enter the path to the file TEMP6: ")
  }

# Read the file
cvs <- read.csv(paste(path2, "/ADNIMERGE.csv", sep = "", collapse = NULL), header=T)
fam <- read.table(paste(path, "/", filename,".fam", sep = "", collapse = NULL), header=F)
file <- read.table(paste(path3,"/TEMP6.fam", sep = "", collapse = NULL), header=F)
# Update the SEX
fam$V5 <- file$V5
# Change the RID and FID names 
parti <- strsplit(as.character(fam$V2), "_")
fam$V2 <- sapply(parti, function(x) paste(x[-1], collapse = "_"))
fam$V1 <- as.integer(sub("^.*_.*_(\\d+)$", "\\1", fam$V2))

# Remove the Samples that are "Not Hisp/Latino"
ccc <- cvs[,c("RID","PTID","PTETHCAT")]
name<- ccc[which(ccc$PTETHCAT!="Not Hisp/Latino"),]
array1 <- fam$V1
uni <- unique(name)
valori_non_comuni <- array1[(array1 %in% uni$RID)]
write.table(valori_non_comuni[order(valori_non_comuni) ],paste(path, "/hispanic.txt", sep = "", collapse = NULL), row.names=F, col.names=F, sep="\t", quote=F) 

# Update the Phenotypes 
m <- merge(fam, cvs[,1:13], by.x = c("V2"), by.y=c("PTID"), all.x=T)
name <- row.names(m[ which( m$DX_bl == "AD"| m$DX_bl == "LMCI"| m$DX_bl == "EMCI"),])
d <- m[name,]
risultato <- fam$V2 %in% d$V2
fam$V6[risultato] <- 2  # case
fam$V6[!risultato] <- 1  # control
write.table(fam[order(fam$V1),], paste(path,"/", filename,".fam", sep = "", collapse = NULL), row.names=F, col.names=F, sep="\t", quote=F) 