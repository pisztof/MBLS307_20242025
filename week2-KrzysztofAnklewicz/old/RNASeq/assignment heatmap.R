library("tidyverse")
library(pheatmap)
setwd("./week2-KrzysztofAnklewicz/RNASeq")

raw_data <- read_csv("../data/BB-RNAseq-input-datasets/Bjornson-et-al-2021-Supp-Table1.csv")
raw_data <- raw_data %>% column_to_rownames(var = "AGI")
filtered <- raw_data[, grepl("Col_nlp20", colnames(raw_data))]
head(filtered)

# Heat map of DEGs only from the 30 minute sample.
overview_res <- read_tsv("./results-deseq2/summary_dgea_deseq2_nlp20.tsv")
overview_res <- overview_res %>% filter(abs(log2FoldChange) >=3)
filtered2 <- filtered %>% filter(row.names(filtered) %in% overview_res$GeneID)
pheatmap(t(scale(t(filtered2))), 
         show_rownames = FALSE, 
         main = paste(nrow(filtered2), "DAGs in Col-0 plants treated with nlp20 (a priori)"),
         angle_col = 0,
         labels_col = c("5 min", "10 min", "30 min", "90 min", "180 min"),
         clustering_distance_rows = "correlation", 
         cluster_cols = FALSE
         )

# Second method - find DAGs from the table, not from the previous question; 
# looking for high absolute means or variability accross time.
filtered3 <- filtered
filtered3["mean"] <- filtered3 %>% rowwise() %>% summarise(Mean = mean(abs(c_across(1:5)))) 
filtered3["stdev"] <- filtered3 %>% rowwise() %>% summarise(stdev = sd(c_across(1:5)))

# scatterplot to get an idea of good cutoff values
ggplot(filtered3, aes(x=stdev, y=mean)) + geom_point() + ggtitle("Scatterplot of means and standard deviations")

filtered3 <- filtered3 %>%  filter(mean>1.5 | stdev>1.5)
head(filtered3)

pheatmap(t(scale(t(filtered3[1:5]))), 
         show_rownames = FALSE, 
         main = paste(nrow(filtered3), "DAGs in Col-0 plants treated with nlp20 (second method)"),
         angle_col = 0,
         labels_col = c("5 min", "10 min", "30 min", "90 min", "180 min"),
         clustering_distance_rows = "correlation", 
         cluster_cols = FALSE
         )

# Third method - filter for genes significantly altered at any timepoint.
filtered4 <- filtered[apply(filtered, 1, function(row) any(abs(row) >= 3)), ]

pheatmap(t(scale(t(filtered4))), 
         show_rownames = FALSE, 
         main = paste(nrow(filtered4), "DAGs in Col-0 plants treated with nlp20 (third method)"),
         angle_col = 0,
         labels_col = c("5 min", "10 min", "30 min", "90 min", "180 min"),
         clustering_distance_rows = "correlation", 
         cluster_cols = FALSE,
)
