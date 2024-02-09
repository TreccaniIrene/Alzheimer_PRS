library(dplyr)

# This script reads a file named "unique_file.txt", filters out data within a specific range of base pair values on chromosome 19, 
# adds information for two SNPs related to APOE-e4 and APOE-e2, and saves the processed data into multiple output files.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    path <- args[1]
    path2 <- args[2]
  } else {
    path <- readline(prompt = "Enter the path: ")
    path2 <- readline(prompt = "Enter the path: ")
  }

# Read the file
base <- read.table(paste(path, "/unique_file.txt", sep = "", collapse = NULL) , header=T)
rmchr19 <- base [which(base['BP'] >= 44400000 & base['BP']<= 46500000),]
rmvar <- base[!(base['BP'] >= 44400000 & base['BP']<= 46500000),]
df <- rmvar[,c('SNP','A1','A2','BETA','SE','P','CHR','BP')]

# Add APOE-e4 and APOE-e2
new_row <- data.frame(
  SNP='rs429358',
  A1= 'C',
  A2='T',
  BETA='-1.350',
  SE='0.027',
  P='0',
  CHR="19",
  BP='45411941'
)

# Add APOE-e4 and APOE-e2
new_row2 <- data.frame(
  SNP='rs7412',
  A1= 'T',
  A2='C',
  BETA='-0.387',
  SE='0.04',
  P='0.000000000000000000000123',
  CHR="19",
  BP='45412079 '
)

# Add the new rows to existing database
df <- rbind(df, new_row)
df <- rbind(df, new_row2)

# Save the results
write.table(rmchr19, paste(path2,"/snp_removed.txt", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 
write.table(df[,c('SNP','A1','A2','BETA','SE')], paste(path,"/snp_ref_prscs.txt", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 
write.table(df, paste(path,"/snp_ref_prscs_PRSice.txt", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 