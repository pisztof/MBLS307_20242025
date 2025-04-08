#14-02-2025
data <- read.table(file = "data/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", header = TRUE, sep = "\t")

#load the required libraries
library(ggplot2)
library(tidyverse)
library(ggrepel)

#create volcano plot
data %>% ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP))) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 2, linetype = "solid", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "red") +
  labs(title = "Volcano plot of PAD4 vs YFP",
       x = "log2 abundance PAD4 vs YFP",
       y = "-log10(p-value)") +
  geom_text_repel(data = data[(data$pval_PAD4_vs_YFP < 0.05) & abs(data$log2_PAD4_vs_YFP) > 10,], 
            aes(label = Protein.IDs),
            size = 1.5,
            box.padding = 0.1) +
  theme_classic() +
  ggsave("output/volcano_plot.png", width = 8, height = 6, dpi = 300)

#	Prepare data for boxplot
longer_data <- data %>% pivot_longer(cols = starts_with("Norm_abundance"), names_to = "Sample", values_to = "abundance")
longer_data <- longer_data %>% mutate(Sample = ifelse(grepl("PAD4", Sample), "PAD4_group", "YFP_group"))

# Create boxplot
ggplot(longer_data %>% filter(grepl("EDS1", Fasta.headers)), aes(x=Sample, y=abundance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic() +
  labs(title=paste(longer_data$Protein.IDs,"; log2 ratio",longer_data$log2_PAD4_vs_YFP,"; p-value",longer_data$pval_PAD4_vs_YFP),
       y="Normalised Abundance",
       x="Sample") +
  ggsave("output/boxplot.png", width = 8, height = 6, dpi = 300)