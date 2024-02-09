library(ggplot2)

# This script reads in the results of a PRS analysis and generates a bar plot to visualize 
# the relationship between different p-value thresholds and the corresponding PRS R-squared values.
# The plot is customized using ggplot2 and saved as "clumping.bar.png".

# Read in the files
prs.result <- read.table( "/data/users/itreccani/ADNI/PRS/plink.prs.result.txt", header = T) 
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[ !is.na( prs.result$print.p ) &
                      prs.result$print.p == 0] <- format(prs.result$P[!is.na(prs.result$print.p) &
                                                                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(aes(label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T)  +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)
  ) +
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size =
                                  18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot
ggsave( "clumping.bar.png", height = 7, width = 7)
