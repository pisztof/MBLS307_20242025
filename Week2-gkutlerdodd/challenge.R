library(tidyverse)
library(ggplot2)

#Process data
data <- read.table("Data/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", sep="\t", header=TRUE)
glimpse(data)
data <- data %>% mutate(
  logpval = -log10(pval_PAD4_vs_YFP)
)

#Volcano plot
ggplot(data, aes(x=log2_PAD4_vs_YFP, y=logpval)) +
  geom_point() +
  geom_vline(xintercept=c(-2,2), color="red") +
  geom_hline(yintercept=2, color="red") +
  labs(x="log2 abundance PAD4 vs YFP", y="-log10(p-value)", title="Volcano plot PAD4 vs YFP") +
  annotate("text", -13.5, 4.5, label="AT4G30920.1",size=2.5) +
  annotate("text", 11, 8.6, label="AT3G48090.1;AT3G48090.2", size=2.5) +
  annotate("text", 16, 6.2, label="AT3G52430.1", size=2.5) +
  annotate("text", 12, 5, label="AT3G48080.1", size=2.5) +
  theme_classic()

ggsave("Plots/volcanoplot.png")

#Tidy the data
data_pivot <- data %>% pivot_longer(
  cols = starts_with("Norm_abundance"),
  names_to = "Protein",
  names_pattern = "[^_]*_[^_]*_([^_]*)",
  values_to = "Norm_abundance")

#Abundance boxplot for AT4G30920.1
ggplot(subset(data_pivot, Protein.IDs == "AT4G30920.1"), aes(x=Protein, y=Norm_abundance)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x="Sample", y="Normalized abundance", title="Normalized abundance of AT4G30920.1, log2 ratio -14.4, p=1.47e-05") +
  theme_classic()

ggsave("Plots/abundance_boxplot.png")

