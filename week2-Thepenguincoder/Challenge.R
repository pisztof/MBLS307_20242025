#Packages
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)

#Importing the data
challenge <- Challenge <- read.delim("week2-Thepenguincoder/data/Challenge.txt")

#Create a new column with the transformed p values
challenge$logpval <- -log10(challenge$pval_PAD4_vs_YFP)

#Creates a new column that specifies for each point whether or not its significant
challenge$updownsig <- "not significant"
challenge$updownsig[challenge$log2_PAD4_vs_YFP > 2 & 
                      challenge$logpval > -log10(0.05/length(challenge$Protein.IDs))]<- "upregulated"
challenge$updownsig[challenge$log2_PAD4_vs_YFP < -2 & 
                      challenge$logpval > -log10(0.05/length(challenge$Protein.IDs))]<- "downregulated"
challenge$updownsig<-factor(challenge$updownsig, levels=c("upregulated","downregulated","not significant"))

#Makes a small dataframe with the most up or downregulated genes, used later for annotation
top_genes <- challenge[(challenge$log2_PAD4_vs_YFP > 10|challenge$log2_PAD4_vs_YFP < -10) & challenge$logpval > -log10(0.05/length(challenge$Protein.IDs))
                       ,]

#Creating the plot and exporting it
volcanoplot <- ggplot(data = challenge, aes(x=log2_PAD4_vs_YFP, y=logpval, col=updownsig))+
  geom_point()+
  labs(color = 'Significance', #legend_title, 
       x = expression("log"[2]*" abundance PAD4 vs YPF"), y = expression("-log"[10]*"p-value")) +
  geom_vline(xintercept = c(-2, 2), col = "red") +
  geom_hline(yintercept = -log(0.05/length(challenge$Protein.IDs), base=10), col = "red")+
  scale_color_manual(values = c("red", "blue", "grey"))+
  geom_text(data = top_genes, aes(label = Majority.protein.IDs), vjust = -1, size = 2)+
  ggtitle('Volcano plot PAD4 vs YFP')

volcanoplot
ggsave(filename = "volcanoplot.pdf", plot = volcanoplot, device = "pdf", 
       path = "week2-Thepenguincoder/data", width = 250, height = 150, units = "mm", 
       dpi = 500, limitsize = TRUE)



#Select a protein from the dataset
protein <- challenge[30,]

#Pivot the dataframe to make later analysis easier
protein_long <-protein %>% pivot_longer(
  cols = contains("rep"), 
  names_to = "protein", 
  values_to = "abundance"
)

#Make a factor from the two proteins
protein_long$protein[grep("PAD4",protein_long$protein)] <- "PAD4"
protein_long$protein[grep("YFP",protein_long$protein)] <- "YFP"
protein_long$protein <- factor(protein_long$protein)

#Creates and exports the bloxplot
boxplot <- ggplot(aes(y=abundance, x=protein), data=protein_long)+
  geom_boxplot()+
  geom_jitter(width=0.1)+
  xlab("Protein") +
  ylab("Abundance") +
  theme_classic()+
  ggtitle(paste('Normalized abundance of', toString(protein_long$Majority.protein.IDs[1]),
                "\nlog2 ratio:", toString(protein_long$log2_PAD4_vs_YFP[1]),
                "p=", toString(protein_long$pval_PAD4_vs_YFP[1])))

boxplot
ggsave(filename = "challenge_boxplot.pdf", plot = boxplot, device = "pdf", 
       path = "week2-Thepenguincoder/data", width = 250, height = 150, units = "mm", 
       dpi = 500, limitsize = TRUE)



