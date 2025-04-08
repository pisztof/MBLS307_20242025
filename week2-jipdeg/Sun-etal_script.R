#import data
data <- read.table('week2-jipdeg/Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt',header = TRUE, sep='\t')
data$pval.neglog <- -log10(data$pval_PAD4_vs_YFP)

#make a volcano plot with annotated points if PAD4 or YFP has 2^10 more abundance than the other in a protein sample. 
library(tidyverse)
ggplot(data, aes(log2_PAD4_vs_YFP, pval.neglog))+
  geom_point()+
  geom_point(data=data %>% filter(abs(log2_PAD4_vs_YFP)>10), aes(log2_PAD4_vs_YFP, pval.neglog),colour='red')+
  geom_text(aes(label=ifelse(abs(log2_PAD4_vs_YFP)>10,as.character(Protein.IDs),'')),vjust=-1)+
  xlab('log2 abundance PAD4 vs YFP')+
  ylab('-log10(p-value)')

#make boxplots for the different points that were annotated
abund.dif <- subset(data, abs(log2_PAD4_vs_YFP)>10)

#extracting all the rows into their own dataframes to plot with
abund.dif.proteins <- list()
for (i in 1:7){
  abund.dif.proteins[[i]] <- abund.dif[i,] %>%
    pivot_longer( cols = starts_with("Norm_abundance_")) %>%
    mutate( name = ifelse(grepl("PAD4", name), 'PAD4', 'YFP'))
  
}

#Plotting all the dataframes and saving the plots
for(i in 1:7){
  df = abund.dif.proteins[[i]]
  ggplot(df, aes(x = name, y = value)) +
    geom_boxplot() +
    ylab('Normalized abundance') +
    xlab('Sample') +
    ggtitle(paste('Normalized abundance of ', df$Protein.IDs[1], 
                  ', log2 ratio: ', df$log2_PAD4_vs_YFP[1], 
                  ', p-value: ', df$pval_PAD4_vs_YFP[1], sep=''))
  ggsave(paste(df$Protein.IDs[1], '_norm_abundance.png'), device = 'png')
}