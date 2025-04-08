#Load library
library(tidyverse)

# Load data
data <- read.table("data/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", header=TRUE, sep="\t")

# Check the data
str(data)

# Subset most differentially abundant proteins
proteins <- subset(data, log2_PAD4_vs_YFP > 9)

# Volcano plot
pdf("output/challengevolcanoplot.pdf")
data %>% ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP))) +
  geom_point() +
  geom_text(aes(label = Protein.IDs), data = proteins, vjust = "outward", nudge_y = 0.1, size = 2) +
  geom_vline(xintercept = c(-2, 2), color = "cyan") +
  geom_hline(yintercept = 2, color = "cyan") +
  labs(title="Volcano plot", x="log2 Abundance PAD4 vs YFP", y="-log10(p-value)") +
  theme_minimal()
dev.off()

#Choose a random protein
protein <- data[111,]
print(protein)

# Transform data for boxplot
boxplotdata <- protein %>% pivot_longer(cols=starts_with("Norm_abundance"), names_to="Sample", values_to="Abundance")
boxplotdata <- boxplotdata %>% mutate(Sample = ifelse(grepl("PAD4", Sample), "PAD4_group", "YFP_group"))
print(boxplotdata)

# Boxplot
pdf("output/challengeboxplot.pdf")
boxplotdata %>% ggplot(aes(x=Sample, y=Abundance)) + 
                        geom_boxplot(outlier.shape = NA) +
                        geom_jitter(width = 0.2, color = "cyan") +
                        labs(title = "Normalized abundance of AT1G14810.1, log2 ratio  0.598, p=0.112", 
                        x = "Sample", 
                        y = "Normalized abundance") +
                        theme_minimal()
dev.off()


