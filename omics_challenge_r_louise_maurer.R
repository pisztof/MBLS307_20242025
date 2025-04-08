plant <- read.table(file="C:/Users/mloui/Downloads/omics_r_challenge.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
str(plant)
summary(plant)
#plant2 <- read.table(file="C:/Users/mloui/Downloads/plant2.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
#str(plant2)
#identical(plant,plant2)

#install tidyverse package
install.packages("tidyverse")
install.packages("ggrepel")
library(tidyverse)
library(ggrepel)

#making a volcano plot with ggplot2

#plant$plantlabel <- NA
#plant$plantlabel[plant$log2_PAD4_vs_YFP > 10] <- plant$Protein.IDs[plant$log2_PAD4_vs_YFP > 10]

#add a column with labels for the most differentially abundant proteins
plant <- plant %>% mutate(plantlabel = if_else(log2_PAD4_vs_YFP > 10, Protein.IDs, NA))

#make the volcano plot, add labels and lines
ggplot(data=plant, aes(x=log2_PAD4_vs_YFP, y=-log10(pval_PAD4_vs_YFP), label= plantlabel)) + 
  geom_point(shape=16) + 
  theme_classic() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ggtitle("Volcano plot PAD4 vs YFP") +
  xlab("log2 abundance PAD4 vs YFP") +
  ylab("-log10(p-value)") +
  geom_text_repel()

head(plant)
#subset the dataframe
plant_bp <- plant[1,1:11]
plant_bp

#pivot longer to tidy the data, assign it to subset
subset <- plant_bp %>% 
  pivot_longer(
    cols = starts_with("Norm_abundance"),
    names_to = "replicate",
    values_to = "Norm_ab"
  )

#add a column with PAD4 or YFP
subset <- subset %>% mutate(protein = if_else(grepl("YFP",subset$replicate), "YFP", "PAD4"))

ggplot(data=subset, aes(x=protein, y=Norm_ab)) +
  geom_boxplot() +
  theme_classic() +
  xlab("Sample") +
  ylab("Normalized abundance") +
  
  geom_jitter()



  
  


