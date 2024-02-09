library(plinkQC)

# Evaluates and depicts results of plink –pca (via run_check_ancestry or externally conducted pca) 
# on merged genotypes from individuals to be QCed and individuals of reference population of known genotypes.
# Currently, check ancestry only supports automatic selection of individuals of European descent. 
# It uses information from principal components 1 and 2 returned by plink –pca to find the center of the European reference samples 
# (mean(PC1_europeanRef), mean(PC2_europeanRef). It computes the maximum Euclidean distance (maxDist) of the European reference samples from this centre. 
# All study samples whose Euclidean distance from the centre falls outside the circle described by the radius 
# r=europeanTh* maxDist are considered non-European and their IDs are returned as failing the ancestry check. 
# check_ancestry creates a scatter plot of PC1 versus PC2 colour-coded for samples of the reference populations and the study population. 

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  path <- args[1]
  name <- args[2]
} else {
  path <- readline(prompt = "Enter the path to the PCA directory: ")
  name <- readline(prompt = "Enter the name of the file: ")
  
}

indir <- system.file("extdata", package="plinkQC")
refname <- 'all.phase3'
prefixMergedDataset <- paste (name,'.all.phase3', sep="", collapse = NULL)

exclude_ancestry <- evaluate_check_ancestry(indir=path, 
                                            name=name,
                                            qcdir=path,
                                            prefixMergedDataset=prefixMergedDataset,
                                            refSamplesFile=paste(system.file("extdata", package="plinkQC"), "/Genomes1000_ID2Pop.txt", sep=""),
                                            refColorsFile=paste(system.file("extdata", package="plinkQC"), "/Genomes1000_PopColors.txt", sep=""),
                                            refPopulation = c("CEU", "TSI","FIN", "GBR"),
                                            europeanTh = 2,
                                            interactive=T)

# Save the people not belonging to European population
write.table(exclude_ancestry[["fail_ancestry"]], paste( path,"/remove.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = F)