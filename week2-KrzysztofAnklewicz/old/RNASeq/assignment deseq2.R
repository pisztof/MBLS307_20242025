library("tidyverse")
library("DESeq2")

### Cleaning and formatting the data

# read raw read counts and the design matrix
raw_counts <- read_tsv("../data/BB-RNAseq-input-datasets/PRJEB25079_RawCounts_reformatted.tsv")
design <- read_tsv("../data/BB-RNAseq-input-datasets/design-matrix.tsv")


# reformat the read count table for DEseq2 so the gene names are not in the column but in row names
raw_counts <- raw_counts %>% column_to_rownames(var = "GeneID")
view(raw_counts)


keep_sample <- c("Col_mock_030", "Col_mock_0301", "Col_mock_0302", "Col_mock_0303", "Col_nlp20_030", "Col_nlp20_0301", "Col_nlp20_0302", "Col_nlp20_0303")
keep_factor <- c("Treatment")
DESeq2_design <- design %>% column_to_rownames(var = "SampleID")

DESeq2_design <- DESeq2_design[keep_sample, ] %>% select(all_of(keep_factor))
dds <- DESeqDataSetFromMatrix(countData = raw_counts[, keep_sample],
                              colData = DESeq2_design,
                              design = ~ Treatment)
dim(dds) # dimensions
head(rownames(dds)) # gene names
head(colnames(dds)) # sample names

hist(log2(counts(dds)), main = "Distribution of the read counts across the dataset before removing low expressed genes", xlab = "log2 raw read count")

keep <- rowSums(counts(dds) >= 5) >= 3 #genes with more than 4 counts in more than 2 samples
dds <- dds[keep,]
transformed_counts <- vst(dds, blind = FALSE)


### DGEA

padj_cutoff <- 0.05
lfc_cutoff <- 1
dds2 <- DESeq(dds, test = "Wald") 
res <- results(dds2, contrast = c("Treatment", "nlp20", "mock"), 
               pAdjustMethod = "BH", alpha = padj_cutoff) # base for the comparison at the end
overview_res <- data.frame(GeneID = row.names(res), log2FoldChange = res$log2FoldChange, p_adj = res$padj)
# write deseq2 results
write_tsv(overview_res, "./results-deseq2/summary_dgea_deseq2_nlp20.tsv")


# Count DEGs
count_up_deseq2 <- overview_res %>% filter(log2FoldChange >= lfc_cutoff & p_adj <= padj_cutoff) %>% nrow()
count_down_deseq2 <- overview_res %>% filter(log2FoldChange <= -lfc_cutoff & p_adj <= padj_cutoff) %>% nrow()
count_total_deseq2 <- overview_res %>% nrow()

print(paste(
  count_up_deseq2, "upregulated genes,", 
  count_down_deseq2, "downregulated genes,",
  count_total_deseq2, "total genes."
))




### PCA plot

for (exp_factor in c("Treatment")){
  pcaData <- plotPCA(transformed_counts, intgroup = exp_factor, ntop = nrow(transformed_counts), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  g <- ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(group)))
  g <- g + geom_point(size = 4, alpha = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_classic()
}
g

# make plot
g <- ggplot(overview_res, aes(x = log2FoldChange, y = -log10(p_adj)))
g + geom_point() +
  geom_hline(yintercept=-log10(0.05), color = "red") + 
  geom_vline(xintercept=1, color = "red") +
  geom_vline(xintercept=-1, color = "red") +
  # give title to the plot, change names of axes
  labs(title = paste(count_total_deseq2, "genes tested,", count_up_deseq2, "up- and", count_down_deseq2,
                     "downregulated; |LFC|>=1, p_adj <= 0.05"),
       y = "-log10(p-value)", x = "log2 FC mock vs. nlp20, 30 min") +
  theme_classic()

