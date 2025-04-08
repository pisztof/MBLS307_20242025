# Week 2 challenge
# Volcano plot and boxplot of relative protein
# quantification data from a proteomics experiment.

# > str(proteomics_data)
# 'data.frame':	1877 obs. of  13 variables:
#   $ Fasta.headers           : chr  "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" "AT1G01320.2 | Symbols:  | Tetratricopeptide repeat (TPR)-like superfamily protein | chr1:121582-130099 REVERSE "| __truncated__ "AT1G01610.1 | Symbols: ATGPAT4, GPAT4 | glycerol-3-phosphate acyltransferase 4 | chr1:221950-224255 REVERSE LENGTH=503" ...
# $ Protein.IDs             : chr  "AT1G01080.1;AT1G01080.2" "AT1G01090.1" "AT1G01320.2;AT1G01320.1" "AT1G01610.1" ...
# $ Majority.protein.IDs    : chr  "AT1G01080.1;AT1G01080.2" "AT1G01090.1" "AT1G01320.2;AT1G01320.1" "AT1G01610.1" ...
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

library(tidyverse)

# Fri Feb 14 15:55:28 2025 ------------------------------

proteomics_data <- read.table("data/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", 
                              header = TRUE, sep = "\t")
diff_abundant <- subset(proteomics_data, proteomics_data$log2_PAD4_vs_YFP > 10)

png("output/challenge_volcano_plot.png")
proteomics_data %>% 
  ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP))) +
    geom_point() +
    geom_text(aes(x = log2_PAD4_vs_YFP, 
                  y = -log10(pval_PAD4_vs_YFP), 
                  label = Protein.IDs), 
              data = diff_abundant,
              vjust = "outward", nudge_y = 0.1, size = 2) +
    geom_vline(xintercept = c(-2, 2), color = "red") +
    geom_hline(yintercept = 2, color = "red") +
    labs(title = "Volcano plot PAD4 vs YFP", 
         x = "log2 abundance PAD4 vs YFP", 
         y = "-log10(p-value)") +
    theme_classic()
dev.off()



protein <- proteomics_data[101,]
plot_data <- protein %>% 
  pivot_longer(cols = starts_with("Norm"),
               names_to = c(NA, NA, "Sample","Replicate"),
               names_sep = "_",
               values_to = "Normalized.abundance")

png("output/challenge_boxplot.png")
plot_data %>% 
  ggplot(aes(x = Sample, y = Normalized.abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  labs(title = "Normalized abundance of AT1G13190.1, log2 ratio -2.28, p=0.0188", 
       y = "Normalized abundance") +
  theme_classic()
dev.off()
