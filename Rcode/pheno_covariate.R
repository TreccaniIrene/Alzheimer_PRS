
# This R script reads a .fam file and extracts specific columns to create a new covariate file.
# It also reads a PCA eigenvec file to extract principal components (PCs) and merge them with the covariate information.
# The covariate file is then saved as "covariate.txt" with the first six columns.
# Additionally, the script creates a new phenotype file by extracting specific columns from the .fam file and saves it as "recodedpheno.txt".

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
  path2 <- args[2]
  filename <- args[3]
} else {
  path <- readline(prompt = "Enter the path to the directory: ")
  path2 <- readline(prompt = "Enter the path to the file: ")
  filename<- readline(prompt = "Enter the file: ")
}

file <- read.table(paste( path,"/",filename,".fam", sep = "", collapse = NULL), header=F)
# Create a new file with information of covariate
covariate <- file[,c(1,2,5)]
colnames(covariate) <- c("FID","IID", "SEX")
pcs <- read.table(paste(path2, "/PCA.eigenvec", sep = "", collapse = NULL ), header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov[, 1:6],paste(path,"/covariate.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = F)

# Create a new file with information of phenotype
pheno<-file[,c(1,2,6)]
colnames(pheno)<- c("RID","FID","PHENO")
write.table(pheno, paste(path,"/recodedpheno.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = T)