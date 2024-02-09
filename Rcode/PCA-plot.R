# Load necessary libraries
library(ggplot2) # For ggplot functionality
library(plotly)  # For interactive plots

# This code handles the reading of .psam and .eigenvec files for PCA analysis.
# If command line arguments are provided, it sets the file paths and the name of the .eigenvec file.
# Otherwise, it prompts the user to manually input the file paths.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  path <- args[1]
  path2 <- args[2]
  filename <- args[3]
} else {
  path2 <- readline(prompt = "Enter the path to the file .psam: ")
  path <- readline(prompt = "Enter the path to the file .eigenvec: ")
  filename <- readline(prompt = "Enter the file name .eigenvec: ")
}

# Read data from files
all_phase <- read.table(paste(path2, "/all_phase3.psam", sep = "", collapse = NULL), header = FALSE)
pca <- read.table(paste(path, "/", filename, ".eigenvec", sep = "", collapse = NULL), header = FALSE)

# Merge data
m <- merge(pca, all_phase, by.x = "V2", by.y = "V1", all = TRUE)
df <- m[3:4]

# Rename columns
colnames(m)[colnames(m) == "V5.y"] <- "Population"
colnames(m)[colnames(m) == "V3.x"] <- "PC1"
colnames(m)[colnames(m) == "V4.x"] <- "PC2"
colnames(m)[colnames(m) == "V5.x"] <- "PC3"

# Set row names
row.names(m) <- m$V2
name <- row.names(m[which(is.na(m$Population == "NA")), ])
m$Population[which(m$V2 %in% name[])] <- "ADNI"

# Perform PCA
pca_res <- prcomp(df, scale = TRUE)

# Create 2D plot using ggplot
p <- ggplot(data = m, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Population))

# Save 2D plot as JPEG
ggsave(paste(path, "/PCA_plot.jpg", sep = "", collapse = NULL),
       plot = p, width = 6, height = 4, units = "in", dpi = 300)

# Create 3D plot using plotly and save as HTML
plot_ly(x = m$PC1, y = m$PC2, z = m$PC3, type = "scatter3d", mode = "markers", color = m$Population) %>%
  htmlwidgets::saveWidget(paste(path, "/PCA3D_plot.html", sep = "", collapse = NULL))
