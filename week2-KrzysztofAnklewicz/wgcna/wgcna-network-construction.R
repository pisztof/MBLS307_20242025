### weighted gene co-expression network analysis
### for RNAseq data from Bjornson et al., 2021 (Nature Plants)

### install and load libraries
# BiocManager::install("impute") # dependency for WGCNA
# BiocManager::install("preprocessCore") # dependency for WGCNA
library("WGCNA") #install.packages("WGCNA"); first install dependencies above
library("tidyverse")
library("readxl")

### load and select the input expression data
# load log2FC data for PAMP treatments, values are relative to time point 0 min
input_lfc <- read_excel("input-data/Bjornson-etal-2021-STable1.xlsx",
                        sheet = "Bjornson_Supp_Table1",
                        na = c("NA", ""))

# load FDR-corrected p-values for the lfc values above
input_fdr <- read_excel("input-data/Bjornson-etal-2021-STable2.xlsx",
                        sheet = "Bjornson_Supp_Table2",
                        na = c("NA", ""))

# check format of the input data
head(input_lfc)
head(input_fdr)
str(input_lfc)
str(input_fdr)

# limit the dataset to the wild type Col-0 samples across all PAMP and mock treatments
# and genes differentially expressed in at least one of the samples
# (|log2FC| >= 3 at least in one of the samlples, p_adj <= 0.01)

genes_with_padj <- input_fdr %>%
  select(contains("Col_") | contains("AGI")) %>%
  filter(if_any(everything(), ~ . <= 0.01)) %>%
  select("AGI") %>%
  pull() %>%
  unique()

input_subset <- input_lfc %>% select(contains("Col_") | contains("AGI")) %>%
  filter(if_any(where(is.numeric), ~ abs(.) > 3)) %>%
  filter(AGI %in% genes_with_padj) %>%
  column_to_rownames(var = "AGI")

### calculate correlation coefficients for genes based on the lfc values
# WGCNA package expects the input data to have geneIDs as column names and sample as rownames
# we will therefore transpose the input table first
input_subset <- t(input_subset)
str(input_subset)

# allow multithreading on your laptop to speed up the calculations
allowWGCNAThreads()

# In WGCNA, we build a coexpression (correlation) based similarity matrix.
# After that, it is converted into the adjency matrix, informing the network layout.
# The conversion is made by raising the correlation coefficient into a power of something.
# We will select a couple of soft-thresholding powers and check how the obtained
# initial network fits the scale-free topology model (power law).
# We will choose the soft threshold using diagnostic plots.

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 30, by = 4))
powers

# Call the network topology analysis function
# The function will calculate the similarity coefficients between the genes
# and then raise them to the powers above. Finally, the function will test
# fit of the observed values in the adjacency to what is expected for the scale-free topology
# (power law distribution)
sft <- pickSoftThreshold(
  input_subset, # <= Input data
  #blockSize = 30, # used if your laptop cannot handle calculating the large number of pairwise similaritys
  powerVector = powers,
  verbose = 5
)

# plots for scales
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# How good is the fit to the scale-free topology?
# pick a power value providing a combination of a fit to the scale-free topology
# and not too many connections between the nodes.
picked_power <- 24

# Now, we will use the function blockwiseModules to calculate the correlation matrix,
# convert it to adjency matrix with the power coefficient from above,
# and form modules of tightly connected nodes/genes.

temp_cor <- cor # set aside the default correlation function from R. We will put it back at the end of module detection step.

cor <- WGCNA::cor # Force to use WGCNA cor function (fix a namespace conflict issue)

netwk <- blockwiseModules(input_subset, # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power, # <= power here
                          networkType = "signed", # both positive and negative r values are relevant
                          
                          # == Tree and Block Options ==
                          deepSplit = 4,
                          pamRespectsDendro = F,
                          #detectCutHeight = 0.95,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          # TOM = Topological Overlap Matrix (TOM)
                          # it is used to cluster genes not based on the coexpression only,
                          # but also using the adjancency matrix based on the coexpression values.
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

# check how many modules were detected
table(netwk$colors)

# Return cor function to original namespace
cor <- temp_cor

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors.
# If you want to try different splitting parameters for the module detection,
# you can change the options deepSplit (4 is the max) and detectCurHeight
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


# extract the module assignment per gene
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

head(module_df)

# write results into a file
write_tsv(module_df, "gene_modules.tsv")

### Now, we will try to interpret results of the module detection step:
# Does the module detection reflect what we see in the expression profiles?
# For this, we will calculate eigenvalues per module and plot correlation
# of the first one with genes in the clusters on the heatmap
# against the samples.
# Get Module Eigengenes per cluster ("characteristic gene")
MEs0 <- moduleEigengenes(input_subset, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add sample names
MEs0$sample = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-sample) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=sample, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Sample Relationships", y = "Modules", fill="corr")

# what are the differences and similarity between the clusters?

### At the next step, we will visualize profiles of the genes in every cluster.
# In the code below, we first test extract gene expression values (lfc's) for the modules of interest
# pick out modules of interest here
modules_of_interest = c("grey", "turquoise", "blue", "brown")

# Pull out list of genes in that modules
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
input_subset[1:5,1:10]


subexpr = t(input_subset)
subexpr = subexpr[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  scale_color_manual(values = c("blue" = "blue", "brown" = "brown", "grey" = "black", "turquoise" = "turquoise")) +
  facet_grid(rows = vars(module)) +
  labs(x = "sample",
       y = "log2FC relative to 0 min")


### Recalculate TOM for modules of interest
# Strictly speaking, this is not required, because we have chosen all modules
# but let's keep the workflow general
TOM = TOMsimilarityFromExpr(input_subset,
                            power = picked_power,
                            networkType = "signed",
                            TOMType = "signed")
hist(TOM)

### Now, we will extract the similarity matrix from the TOM object into a dataframe
# that can be used as an input for Cytoscape
# Add gene names to row and columns
row.names(TOM) = colnames(input_subset)
colnames(TOM) = colnames(input_subset)

# extract the similarity matrix into a form suitable for Cytoscape
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, similarity = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

### Quality control of the extracted annotated values from the similarity matrix
head(edge_list)

par(mar = c(5,5,5,5))
hist(edge_list$similarity)

# let's double-check whether the within module similarity is higher than between modules.
g <- ggplot(edge_list %>% mutate(module_comb = paste0(module1, "_", module2)),
            aes(x = module_comb, y = similarity))

g + geom_boxplot(outlier.color = NA) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  labs(title = "similarity between nodes in the modules",
       x = "Module combination")

# For Cytoscape, we will reduce connectivity in the network to allow its quicker
# visualization.
edge_list <- edge_list %>% filter(abs(similarity) > 0.4)

### save the network as a TSV file
write_tsv(edge_list,
          file = "edgelist.tsv")

### open the network in Cytoscape

###### enriching the network with additional sources of information
### add information about the known protein-protein interactions (PPIs)
# load PPI data from BioGrid
# search what BioGrid is https://thebiogrid.org and how it differs from StringDB
library("vroom")
ppi <- vroom("./input-data/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-4.4.243.tab3.txt")

# select gene ids on the network to limit the PPI document
gene_ids <- unique(c(edge_list$gene1, edge_list$gene2))

# select relevant PPIs
ppi <- ppi %>% filter(`Systematic Name Interactor A` %in% gene_ids & `Systematic Name Interactor B` %in% gene_ids) %>%
  filter(`Systematic Name Interactor A` != `Systematic Name Interactor B`) %>%
  select(`Systematic Name Interactor A`, `Systematic Name Interactor B`) %>%
  distinct()

# to facilitate adding the PPI information to the network, we will create
# new columns with gene pairs merged in one string; since PPIs are not directional,
# we will create two pairs, 1 and 2, but altering the relative positions
# of the interactors A and B
ppi <- ppi %>% mutate(Gene_pair1 = paste(`Systematic Name Interactor A`, `Systematic Name Interactor B`),
                      Gene_pair2 = paste(`Systematic Name Interactor B`, `Systematic Name Interactor A`)) %>%
  pivot_longer(cols = c(Gene_pair1, Gene_pair2), values_to = "Gene_pair", names_to = "Orientation") %>%
  select(Gene_pair) %>%
  pull()

# create a similar temporary merge column in the network table,
# add a new column with PPI if the PPI is found in BioGrid,
# remove the temporary helper column.
edge_list_annotated <- edge_list %>% mutate(Gene_pair = paste(`gene1`, `gene2`)) %>%
  mutate(PPI = if_else(Gene_pair %in% ppi, "PPI", "")) %>%
  select(!Gene_pair)

### add information about the log2FC to the network
# we will add log2FC in Col-0 after the PAMP nlp20 treatment as annotations to the nodes
edge_list_annotated <- left_join(edge_list_annotated, input_lfc[, c("AGI", "Col_nlp20_180")],
                                 by = c("gene1" = "AGI")) %>%
  rename(gene1_Col_nlp20_180_lfc = Col_nlp20_180)

edge_list_annotated <- left_join(edge_list_annotated, input_lfc[, c("AGI", "Col_nlp20_180")],
                                 by = c("gene2" = "AGI")) %>%
  rename(gene2_Col_nlp20_180_lfc = Col_nlp20_180)

# Export Network file to be read into Cytoscape
write_tsv(edge_list_annotated,
          file = "edgelist-annotated.tsv")

### Visualize and analyze the network in Cytoscape