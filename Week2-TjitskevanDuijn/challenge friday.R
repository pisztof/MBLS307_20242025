sun <- read.table(file = 'data/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt', header = TRUE, sep = "\t",
                      stringsAsFactors = TRUE)
library(tidyverse)

str(sun)
# 'data.frame':	1877 obs. of  13 variables:
# $ Fasta.headers           : Factor w/ 1877 levels "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__,..: 1 2 3 4 5 6 7 8 9 10 ...
# $ Protein.IDs             : Factor w/ 1877 levels "AT1G01080.1;AT1G01080.2",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ Majority.protein.IDs    : Factor w/ 1877 levels "AT1G01080.1;AT1G01080.2",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ Norm_abundance_PAD4_rep1: num  22.6 22.4 26.2 24.5 23.1 ...
# $ Norm_abundance_PAD4_rep2: num  23.2 22.5 26.3 24.2 20.8 ...
# $ Norm_abundance_PAD4_rep3: num  22.4 22.4 26.3 24.5 23.5 ...
# $ Norm_abundance_PAD4_rep4: num  23.3 22.3 26.2 23.8 25.2 ...
# $ Norm_abundance_YFP_rep1 : num  19.4 23 24.9 21 22.3 ...
# $ Norm_abundance_YFP_rep2 : num  19.1 22.7 25 20.5 19.8 ...
# $ Norm_abundance_YFP_rep3 : num  20.7 20.1 25.7 19.5 21 ...
# $ Norm_abundance_YFP_rep4 : num  20.7 23.1 25.5 21.2 20.1 ...
# $ log2_PAD4_vs_YFP        : num  2.882 0.212 0.971 3.727 2.313 ...
# $ pval_PAD4_vs_YFP        : num  0.000845 0.77466 0.001873 0.000113 0.074875 ...

#take the negative logarithms of the pvalues and store as new variable
sun$log_pval <- -log10(sun$pval_PAD4_vs_YFP)

# make subset with genes i want to label
significant_genes <- sun[sun$log_pval > 5.5 & abs(sun$log2_PAD4_vs_YFP) > 6, ]

#volplot has the log2 fold change on the x axis, and the pvalue on the y axis
volcano_plot <- ggplot(sun, aes(x = log2_PAD4_vs_YFP, y = log_pval), label = '') +
  geom_point(size = 2) +
  labs(title = "Volcano Plot PAD4 vs YFP", x = "log2 abundance PAD4 vs YFP", y = "-log10(P-value)") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_text(data = significant_genes, aes(label = Protein.IDs), 
            vjust = -0.5, hjust = 1, size = 3, color = "red")
  theme_clear()

# Display the plot
print(volcano_plot)


#make a subset of the data and tidy it
#I choose the gene with the highest pvalue
sun_ordered <- sun[order(sun$log_pval, decreasing = TRUE),]

gene <- sun_ordered[1,] 
gene_long <- gene %>%
  pivot_longer(
    cols = starts_with("Norm_abundance"),  
    names_to = c("condition", "rep"),
    names_pattern = "Norm_abundance_(.*)_rep(\\d)", 
    values_to = "Norm_abundance"
  )

p <- ggplot(gene_long) +
  geom_boxplot(aes(x = condition, y = Norm_abundance)) +
  labs(y = "Normalized abundance", x = "Sample")+
  geom_jitter(aes(y = Norm_abundance, x = condition), colour = "black", height = 0, 
              width = 0.1, alpha = 0.5)+
  ggtitle(label = gene$Protein.IDs, subtitle = paste("log 2 ratio is", gene$log2_PAD4_vs_YFP, "and the pvalue is", gene$pval_PAD4_vs_YFP))
#not very pretty as the data is very clustered. 

ggsave("output/plot.pdf", plot = p, width = 7, height = 5)

#taking a random one now.
gene <- sun_ordered[500,]
#running the gene_long and the ggplot code above
ggsave("output/plot2.pdf", plot = p, width = 7, height = 5)
