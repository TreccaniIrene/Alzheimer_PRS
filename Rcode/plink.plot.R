library(ggplot2)
library(dplyr)
library(magrittr)

# This script reads in PRS data and phenotype information, renames the columns, merges the files,
# and generates a density plot of polygenic scores stratified by phenotype.
# The resulting plot is saved as "Plink_0.5.png".

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
 # path2 <- args[2]
} else {
  path <- readline(prompt = "Enter the path to the directory: ")
 # path2 <- readline(prompt = "Enter the path to the file cvs: ")
}

# Read in the files
prs <- read.table("/data/users/itreccani/ADNI/PRS/prs.plink.0.5.profile", header = T)
phen <- read.table("/data/users/itreccani/ADNI/PRS/recodedpheno.txt", header = T)
                    
# Rename the phen
colnames(phen) <- c("FID","IID", "PHENO")
phen$PHENO <- as.factor(phen$PHENO)
levels(phen$PHENO) <- c("healthy", "diseased")

# Merge the files
dat <- merge(prs, phen,by = c("IID","FID"))
med_salary_df <- dat %>%
  group_by(PHENO.y) %>%
  summarize(median=median(SCORE))

# Start plotting
ggplot(dat, aes(x = SCORE, color = PHENO.y)) +
  geom_density() +
  theme_classic() +
  labs(x = "Polygenic Score", y = "Density") +
  geom_vline(data = med_salary_df, aes(xintercept = median, color = PHENO.y), linewidth = 0.5)

# Save the plot
ggsave( paste( path,"/Plink_0.5.png", sep = "", collapse = NULL), height = 7, width = 7)