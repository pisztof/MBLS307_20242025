library(readr)
library(dplyr)
library(ggplot2)

sun_table <- read.table("Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

#volcano plot
sun_table <- sun_table %>% mutate(log10_pval=-log10(pval_PAD4_vs_YFP))
sun_table %>% ggplot(aes(x=log2_PAD4_vs_YFP, y=log10_pval)) + 
  geom_point() + 
  ggtitle("Volcano plot PAD4 vs YFP") + 
  xlab("log2 abundance PAD4 vs YFP") + 
  ylab("-log10(p-value") + 
  theme_classic() + 
  geom_vline(xintercept = c(-2, 2), color="red") + 
  geom_hline(yintercept = 2, color="red") +
  geom_text(data=sun_table[(sun_table$log2_PAD4_vs_YFP>10),],
            aes(label=Protein.IDs))
#everything before the comma applied for columns and everything after, to rows

png("volcano_plot.png")
sun_table %>% ggplot(aes(x=log2_PAD4_vs_YFP, y=log10_pval)) + 
  geom_point() + 
  ggtitle("Volcano plot PAD4 vs YFP") + 
  xlab("log2 abundance PAD4 vs YFP") + 
  ylab("-log10(p-value") + 
  theme_classic() + 
  geom_vline(xintercept = c(-2, 2), color="red") + 
  geom_hline(yintercept = 2, color="red") +
  geom_text(data=sun_table[(sun_table$log2_PAD4_vs_YFP>10),],
            aes(label=Protein.IDs))
dev.off()


#boxplot
library(tidyverse)

#you normally need to make the table more symmetric to be able to make a plot
ip_abundance <- sun_table %>% filter(grepl("AT1G01080", Protein.IDs)) %>% select(contains("Norm_abundance"))
ip_abundance_longer <- ip_abundance %>% pivot_longer(cols = contains("Norm_abundance"),
                                                     names_to = "sample", values_to = "abundance") 

#as long as you have clear criteria separating each type of data, you will be able to make a plot
ip_abundance_longer <- ip_abundance_longer %>% mutate(sample_group = if_else(grepl("PAD4", sample), "PAD4", "YFP"))

png("boxplot.png")
ip_abundance_longer %>% ggplot(aes(y=abundance, x=sample_group)) + geom_boxplot() + theme_classic() + ggtitle("Normalized abundance of AT1G01080, log2 ratio 2.88, p=0.0009") + 
  xlab("Sample") + 
  ylab("Normalized abundance") + geom_jitter()
dev.off()
