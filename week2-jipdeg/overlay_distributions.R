#============================================
### Metadata and descriptions ###
# Script name "overlay_distributions.R"
# Author: D. Lapin
# Date: 2025-02-10
# Description: the script generates two normal distributions,
# plots them as density plots and saves result into a file.

#============================================
### load libraries ###
library(ggplot2) # visualization of the data

#============================================
### set working directory to the week 2 folder ###

#============================================
### custom functions (if any) ###

#============================================
### data generation or input ###
# generate two distributions and store them in one dataframe
df <- data.frame(Distribution = c(rep("A", 100),
                                  rep("B", 100)),
                 Value = c(rnorm(100, 0, 1),
                           rnorm(100, 2, 1)))

#============================================
### data visualization, cleaning, preprocessing ###
# plot the distributions as density plots
# save plot as a file
g <- ggplot(data = df, aes(x = Value))
g + geom_density(aes(color = Distribution, fill = Distribution), alpha = 0.2) +
  theme_classic()

#============================================
### export data ###
ggsave("density_overlay.png")
