

### load libraries
library(tidyverse)
library(ggplot2)

# to install tidyproteomics, first install these dependencies:
BiocManager::install(c("limma","qvalue","fgsea"))
library(limma)

# then install tidyproteomics from GitHub
devtools::install_github("jeffsocal/tidyproteomics")
library(tidyproteomics)

##############################################################################################################################
##############################################################################################################################
### read input MaxQuant data using the import() function from tidyproteomics
# the code assumes that the input data are in the folder 'data' in your working directory.
data_prot <- "./data/reform-proteinGroups.txt" %>% 
  import('MaxQuant', 'proteins')
data_prot %>% summary("sample")

### assign sample type to each sample
data_prot <- data_prot %>% 
  reassign(sample == 'mcf10a_gpa_1', .replace = 'MCF10A') %>%
  reassign(sample == 'mcf10a_gpa_2', .replace = 'MCF10A') %>%
  reassign(sample == 'mcf10a_gpa_3', .replace = 'MCF10A') %>%
  reassign(sample == 'mcf7_gpa_1', .replace = 'MCF7') %>%
  reassign(sample == 'mcf7_gpa_2', .replace = 'MCF7') %>%
  reassign(sample == 'mcf7_gpa_3', .replace = 'MCF7') %>%
  reassign(sample == 'mda231_gpa_1', .replace = 'MDA231') %>%
  reassign(sample == 'mda231_gpa_2', .replace = 'MDA231') %>%
  reassign(sample == 'mda231_gpa_3', .replace = 'MDA231')

data_prot %>% summary("sample")

##############################################################################################################################
##############################################################################################################################
### To facilitate the data analysis, we will the Gene Ontology annotation information for the proteins.
# We can obtain the GO annotation from the Uniprot database which is often used in proteomics research https://www.uniprot.org
# Go to UniProt KnowledgeBase https://www.uniprot.org/uniprotkb?query=* and select the on left hand side "Human" from the list of popular organisms.
# The human reference proteome has ~205k proteins (March 2025).
# Click on Download -> Select TSV format -> Customize columns by keeping only `Entry Name`, `Protein names`, `Gene Ontology (biological process)`
# You can also select other types of GO annotations, but for the sake of simplicity, we will focus on Biological Processes
# Download the file, unzip it and put in the data folder in your working directory. The download can take a few minutes.
# Check the file name before reading it with the read_tsv function

data_go <- read_tsv("./data/uniprotkb_AND_model_organism_9606.tsv")

# check how the GO annotation file looks like
data_go

# make the protein GO annotation table tidier
# if you have selected GO category other than biological process, you need to change the code a bit
data_go <- data_go %>%
  select(protein = Entry,
         biological_process = `Gene Ontology (biological process)`) %>%
  # separate the GO terms so we get 1/row
  separate_rows(biological_process, sep="\\;\\s") %>%
  # remove the [GO:accession]
  mutate(biological_process = sub("\\s\\[.+", "", biological_process)) %>%
  # pivot to the needed format
  pivot_longer(biological_process,
               names_to = 'term',
               values_to = 'annotation')

# check that the table is indeed tidier
data_go

# append the GO annotation information to the imported MaxQuant data using the tidyproteomics function annotate()
data_prot_go <- data_prot %>% annotate(data_go)

# now one can search GO annotation with protein IDs;
# for instance the protein IDs P68431 and P62805 are associated with the GO term "telomere organization"
data_prot_go$annotations %>% filter(protein %in% c('P68431', 'P62805'))

# we can also find what proteins are associated with a particular GO term, e.g., "apoptotic process"
data_prot_go %>% 
  subset(biological_process == 'apoptotic process')


##############################################################################################################################
##############################################################################################################################
### Now we will go over the basic steps in the differential abundance testing in quantitative proteomics
# The process includes three main steps: imputation, normalization and differential abundance testing.
# Additional steps such as PCA and GO term enrichment analyses can also be performed.


### Naïve data differential abundance testing without imputation and normalization
# PCA plot with the raw log2-transformed LFQ values
# PC1 - 45% variance, PC2 - 31%; samples are grouped according to the cell line type, but some samples deviate from the rest in the group 
data_prot_go %>% plot_pca()

# Make a heatmap
# One can see many missing values
data_prot_go %>% plot_heatmap()

# If we perform differential abundance testing using the t-test for MCF10A and MDA231 cell lines,
# very few proteins will be detected as differentially abundant.
# Missing values in the dataset will also lead to inaccurate p-values and log2FC estimates - see the horizontal line of dots
data_prot_go %>% 
  expression(MCF10A/MDA231, .method = stats::t.test) %>%
  plot_volcano(MCF10A/MDA231,
               significance_column = 'adj_p_value',
               significance_max = 0.01,
               log2fc_min = 2,
               color_positive = 'orange',
               color_negative = 'purple',
               show_lines = FALSE,
               show_fc_scale = FALSE,
               show_pannels = FALSE,
               labels_column = 'protein') +
  labs(title = "Volcano plot", subtitle = "MCF10A vs. MDA231")
ggsave("no-imputation-volcano.pdf")


### Data analysis with imputation and normalization
# The typical differential abundance testing includes: imputation and normalization of the LFQ abundance values
# The imputation replaces missing (zero) values with low, 'random' values drawn based on the measured LFQ abundances of detected but very low abundant proteins.
# Several methods for the imputation were proposed. We will use the recommended random forest-based method which takes into account LFQ values for all proteins
# in the samples and abundance of the given protein across samples (matrix type).
# Imputation can take several minutes. After the imputation is done, the plot_heatmap() function creates the heatmap.
# Comparison of the heatmaps with and without imputation shows that the missing (zero) values are now replaced with the small non-zero values.
# The imputed dataset can be used for the statistical analysis where the dispersion can now be properly estimated for each protein and each sample group
data_prot_go %>% 
  impute(.function = impute.randomforest, method = 'matrix') %>%
  plot_heatmap()


# After the imputation, we usually run the protein abundance normalization.
# The normalization and the differential abundance testing can be performed 
# using the methods implemented in the limma package.
# Limma stands for Linear Models for Microarray and RNA-Seq Data. Limma implements empirical Bayes (EB) moderation to variance estimates -
# it borrows the variance information across features (e.g., proteins), leading to more reliable statistical inference for the low abundant proteins.
# The same idea was used in the DESeq2 package when we analyzed RNAseq data.
# The code below will execute the imputation (because we didn't save its results above), normalization, and differential abundance testing for the lines MCF10A and MDA231
imp_norm_daa <- data_prot_go %>% 
  impute(.function = impute.randomforest, method = 'matrix') %>% 
  normalize(.method = 'limma') %>%
  expression(MCF10A/MDA231, .method = "limma")

# now let's check what has happened with the variance and LFQ distributions across the samples after the imputation and normalization
imp_norm_daa %>% plot_normalization() # note that the mean and variance became very similar between sample groups
imp_norm_daa %>% plot_variation_cv()
imp_norm_daa %>% plot_pca # the samples groups are still well separated but now the variation between samples within the sample groups is reduced
imp_norm_daa %>% plot_heatmap() # the heatmap, as we have seen before no longer displays the missing, zero values.

# How does the volcano plot look like after the imputation and normalization?
# The long horizontal lines have now disappeared and the number of differentially abundant proteins has increased drastically,
# showing that these steps helped us to increase the power the differential gene expression analysis.
imp_norm_daa %>%
  plot_volcano(MCF10A/MDA231,
               significance_column = 'adj_p_value',
               significance_max = 0.01,
               log2fc_min = 2,
               color_positive = 'orange',
               color_negative = 'purple',
               show_lines = FALSE,
               show_fc_scale = FALSE,
               show_pannels = FALSE,
               labels_column = 'protein') +
  labs(title = "Volcano plot", subtitle = "MCF10A vs. MDA231")
ggsave("imputation-volcano.pdf")

# We will finish the differential abundance testing with GO term enrichment analysis using the Biological process as a category and Fisher's exact test.
# The GO term annotation comes from the file we have downloaded from UniProt.
enrich_results <- imp_norm_daa %>%
  enrichment(MCF10A/MDA231, 
             .terms = 'biological_process',
             .method = 'fishers_exact')

# extract results of the GO term enrichment analysis in the tabular format
enrich_results <- enrich_results$analysis$`MCF10A/MDA231`$enrichment$biological_process$data

# plot the results of the GO term enrichment analysis for the significant terms
g <- ggplot(data = enrich_results %>% filter(adj_p_value <= 0.05),
            aes(x = enrichment, y = fct_reorder(annotation, enrichment)))
g + geom_point(aes(fill = -log10(adj_p_value), size = size)) +
  theme_bw() +
  labs(x = "Enrichment level, log2 scale", y = "Annotation")
ggsave("GO_enrichment.pdf")

############################################################################################
### THE END ################################################################################
############################################################################################
