library(ggplot2)
library(dplyr)
library(magrittr)

# This script performs an analysis of polygenic scores (PRS) using PRSice-2 output files and recoded phenotype data.
# It generates density plots for PRS1, PRS2, and the best PRS, distinguishing between healthy and diseased individuals.
# The resulting plots are saved as PNG files for further examination.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
  path2 <- args[2]
} else {
  path <- readline(prompt = "Enter the path to the directory: ")
  path2 <- readline(prompt = "Enter the path to the file cvs: ")
}

# Read in the files
prsice <- read.table(paste( path,"/PRS1.all_score", sep = "", collapse = NULL), header=TRUE)
phen <- read.table(paste( path,"/recodedpheno.txt", sep = "", collapse = NULL), header = T)

# Consider only the Pt at 0.01 and 0.5 and rename the columns
prs<- prsice[,c("FID","IID","Pt_0.01","Pt_0.5")]
colnames(prs) <- c("FID","IID", "PRS1","PRS2")
write.table(prs ,paste( path,"/PRSice-0.01-0.5.txt", sep = "", collapse = NULL), row.names=F, col.names=T, sep="\t", quote=F) 
                  
# Rename the phen
colnames(phen) <- c("FID","IID", "PHENO")
phen$PHENO <- as.factor(phen$PHENO)
levels(phen$PHENO) <- c("healthy", "diseased")

# Perform the Density plot for PRS1
# Merge the files
dat <- merge(prs, phen,by = c("IID","FID"))
med_salary_df <- dat %>%
  group_by(PHENO) %>%
  summarize(median=median(PRS1))
# Start plotting
ggplot(dat, aes(x = PRS1, color = PHENO)) +
  geom_density() +
  theme_classic() +
  labs(x = "Polygenic Score", y = "Density") +
  geom_vline(data = med_salary_df, aes(xintercept = median, color = PHENO), linewidth = 0.5)
# Save the plot
ggsave( paste( path,"/PRSice-PRS1.png", sep = "", collapse = NULL), height = 7, width = 7)


# Perform the Density plot for PRS2
# Merge the files
med_salary_df <- dat %>%
  group_by(PHENO) %>%
  summarize(median=median(PRS2))
# Start plotting
ggplot(dat, aes(x = PRS2, color = PHENO)) +
  geom_density() +
  theme_classic() +
  labs(x = "Polygenic Score", y = "Density") +
  geom_vline(data = med_salary_df, aes(xintercept = median, color = PHENO), linewidth = 0.5)
# Save the plot
ggsave( paste( path,"/PRSice-PRS2.png", sep = "", collapse = NULL), height = 7, width = 7)



# Perform the Density plot for PRSbest
prsbest <- read.table(paste( path,"/PRS1.best", sep = "", collapse = NULL), header=TRUE)
# Rename the phen
colnames(phen) <- c("FID","IID", "PHENO")
phen$PHENO <- as.factor(phen$PHENO)
levels(phen$PHENO) <- c("healthy", "diseased")
# Merge the files
dat <- merge(prsbest, phen,by = c("IID","FID"))
med_salary_df <- dat %>%
  group_by(PHENO) %>%
  summarize(median=median(PRS))
# Start plotting
ggplot(dat, aes(x = PRS, color = PHENO)) +
  geom_density() +
  theme_classic() +
  labs(x = "Polygenic Score", y = "Density") +
  geom_vline(data = med_salary_df, aes(xintercept = median, color = PHENO), linewidth = 0.5)

# Save the plot
ggsave( paste( path,"/PRSice-PRSbest.png", sep = "", collapse = NULL), height = 7, width = 7)