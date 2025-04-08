
install.packages("dplyr")
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)

##Volcano plot of abundance of PAD4 vs YFP
challengetable<- read.table("Week2_Raphael/fridaychallengedata.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

challengetable <- challengetable %>% mutate(log10_pval=-log10(pval_PAD4_vs_YFP))
challengetable%>%ggplot(aes(x=log2_PAD4_vs_YFP, y=log10_pval)) + 
  theme_classic() + 
  geom_point() + 
  ggtitle("Volcano Plot PAD4 vs YFP") + 
  xlab("log2 Abundance PAD4 vs YFP") + 
  ylab("-log10(p)") + 
  geom_text(data=challengetable[(challengetable$log2_PAD4_vs_YFP>9.5),],
            aes(label=Protein.IDs)) +
  geom_hline(yintercept = 2, color="red") +
  geom_vline(xintercept = c(-2, 2), color="red") 
png("Week2_Raphael/volcano_plot_raphael_final.png")



challengetable%>%ggplot(aes(x=log2_PAD4_vs_YFP, y=log10_pval)) + 
  geom_point() + 
  theme_classic() + 
  ggtitle("Volcano plot PAD4 vs YFP") + 
  xlab("log2 abundance PAD4 vs YFP") + 
  ylab("-log10(p)") + 
  geom_text(data=challengetable[(challengetable$log2_PAD4_vs_YFP>10),],
            aes(label=Protein.IDs)) +
  geom_hline(yintercept = 2, color="red") +
  geom_vline(xintercept = c(-2, 2), color="red") 

dev.off()


## Box plot displaying normal abundance of AT1G01080

long_format_abundance <- challengetable %>% 
  pivot_longer(cols = contains("Norm_abundance"), names_to = "sample", values_to = "abundance") %>% 
  filter(!is.na(abundance)) %>% # removing NA values because I had trouble displaying a proper plot
  mutate(sample_group = if_else(grepl("PAD4", sample), "PAD4", "YFP"))

png("boxplotfinal.png", width = 800, height = 600)

ggplot(long_format_abundance, aes(x = sample_group, y = abundance)) + 
  geom_boxplot() +
  theme_classic() + 
  ggtitle("Normalized abundance of AT1G01080, log2 ratio 2.88, p=0.0009") + 
  xlab("Sample") + 
  ylab("Normalized abundance") 

dev.off() 

