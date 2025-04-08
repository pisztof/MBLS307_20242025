#Week 2 Challenge 
#volcano plot 
df_prot <- read.table("Sun_et_al.txt", header = TRUE, sep = "\t")
str(df_prot)

# $ Fasta.headers           : chr  "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" "AT1G01320.2 | Symbols:  | Tetratricopeptide repeat (TPR)-like superfamily protein | chr1:121582-130099 REVERSE "| __truncated__ "AT1G01610.1 | Symbols: ATGPAT4, GPAT4 | glycerol-3-phosphate acyltransferase 4 | chr1:221950-224255 REVERSE LENGTH=503" ...
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

#on the volcano plot 
  #x is log2(Fold Change) --> RATIO BETWEEN GENE EXPRESSION FOR THE 2 GROUPS (GroupB/GroupA) (abundance)
  #y is the -log10(p-value) between the 2 
#aka the last two columns


library(ggplot2)
#create the plot
ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP) ), data = df_prot) + geom_point()

#now we try to add labels for the y and x axis 
ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP) ), data = df_prot) + geom_point() +
  xlab("log2(PAD4 vs YPF Abundance)")+ 
  ylab("-log10(p-value)") 
?geom_text()

#trying to add labels for the most differentially abundance proteins
#add a column of transformed (-log10/p) to the data frame 
df_prot$minus_log10_p <- -log10(df_prot$pval_PAD4_vs_YFP)
# identify the proteins that need to be labelled
library(tidyverse)

#add row of proteins that need to be labelled to the data frame
df_prot <- df_prot %>% mutate(label_prot = if_else(minus_log10_p > 5 & log2_PAD4_vs_YFP > 10, Protein.IDs, NA))
#label_prot <- df_prot[df_prot$minus_log10_p > 10 & df_prot$log2_PAD4_vs_YFP > 5, 2 ]
 

ggplot(aes(x = log2_PAD4_vs_YFP, y = -log10(pval_PAD4_vs_YFP) ), data = df_prot) + geom_point() +
  xlab("log2(PAD4 vs YPF Abundance)")+ 
  ylab("-log10(p-value)") +
  geom_text(aes(label = label_prot), size =3, vjust = 0, nudge_y = 0.5) +
  geom_hline(yintercept = 5, color = "red")+
  geom_vline(xintercept = 10, color = "red")+
  theme_bw()

#try the boxplot +

#first, we need to arrange the data from rows to columns (and tidy it up)

boxplot_data <- df_prot %>%
  pivot_longer(
    cols = matches("Norm_abundance_(PAD4|YFP)_rep\\d+"),  # Select both PAD4 and YFP columns
    names_to = c("Protein", "Replicate"),
    names_pattern = "Norm_abundance_(PAD4|YFP)_rep(\\d+)",  # Extract both protein name and replicate number
    values_to = "Normalized_Abundance"
  ) %>%
  mutate(Replicate = as.numeric(Replicate))  # Convert replicate number to numeric

str(boxplot_data)
# tibble [15,016 × 10] (S3: tbl_df/tbl/data.frame)
# $ Fasta.headers       : chr [1:15016] "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ "AT1G01080.1 | Symbols:  | RNA-binding (RRM/RBD/RNP motifs) family protein | chr1:45503-46789 REVERSE LENGTH=293"| __truncated__ ...
# $ Protein.IDs         : chr [1:15016] "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" ...
# $ Majority.protein.IDs: chr [1:15016] "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" "AT1G01080.1;AT1G01080.2" ...
# $ log2_PAD4_vs_YFP    : num [1:15016] 2.88 2.88 2.88 2.88 2.88 ...
# $ pval_PAD4_vs_YFP    : num [1:15016] 0.000845 0.000845 0.000845 0.000845 0.000845 ...
# $ minus_log10_p       : num [1:15016] 3.07 3.07 3.07 3.07 3.07 ...
# $ label_prot          : chr [1:15016] NA NA NA NA ...
# $ Protein             : chr [1:15016] "PAD4" "PAD4" "PAD4" "PAD4" ...
# $ Replicate           : num [1:15016] 1 2 3 4 1 2 3 4 1 2 ...
# $ Normalized_Abundance: num [1:15016] 22.6 23.2 22.4 23.3 19.4 ...


#select one specific protein to show in the boxplot (the first one)
sel_prot <- "AT1G01080.1;AT1G01080.2"
boxplot_data_final <- boxplot_data %>% 
  filter(str_detect(Protein.IDs, sel_prot))


str(boxplot_data_final)
# tibble [8 × 10] (S3: tbl_df/tbl/data.frame)
# $ Fasta.headers       : chr [1:8] "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" "AT1G01090.1 | Symbols: PDH-E1 ALPHA | pyruvate dehydrogenase E1 alpha | chr1:47705-49166 REVERSE LENGTH=428" ...
# $ Protein.IDs         : chr [1:8] "AT1G01090.1" "AT1G01090.1" "AT1G01090.1" "AT1G01090.1" ...
# $ Majority.protein.IDs: chr [1:8] "AT1G01090.1" "AT1G01090.1" "AT1G01090.1" "AT1G01090.1" ...
# $ log2_PAD4_vs_YFP    : num [1:8] 0.212 0.212 0.212 0.212 0.212 ...
# $ pval_PAD4_vs_YFP    : num [1:8] 0.775 0.775 0.775 0.775 0.775 ...
# $ minus_log10_p       : num [1:8] 0.111 0.111 0.111 0.111 0.111 ...
# $ label_prot          : chr [1:8] NA NA NA NA ...
# $ Protein             : chr [1:8] "PAD4" "PAD4" "PAD4" "PAD4" ...
# $ Replicate           : num [1:8] 1 2 3 4 1 2 3 4
# $ Normalized_Abundance: num [1:8] 22.4 22.5 22.4 22.3 23 ...


ggplot(boxplot_data_final, aes(x = Protein, y = Normalized_Abundance, fill = Protein)) +
  geom_boxplot() +  
  geom_jitter() +  
  scale_fill_manual(values = c("PAD4" = "skyblue", "YFP" = "salmon")) +  
  labs(
    title = "Normalized Abundance of AT1G01080 in PAD4 vs. YFP",
    subtitle = "log2 ratio = 2.88, p = 0.0009",
    x = "Sample (PAD4 vs. YFP)",
    y = "Normalized Abundance"
  ) 
  



  



