library("tidyverse")
library("DESeq2")

# read raw read counts and the design matrix
raw_counts <- read_tsv("../data/BB-RNAseq-input-datasets/PRJEB25079_RawCounts_reformatted.tsv")
design <- read_tsv("../data/BB-RNAseq-input-datasets/design-matrix.tsv")

# reformat the read count table for DEseq2 so the gene names are not in the column but in row names
raw_counts <- raw_counts %>% column_to_rownames(var = "GeneID")

keep_sample <- c("Col_mock_030", "Col_mock_0301", "Col_mock_0302", "Col_mock_0303", 
                 "Col_mock_180", "Col_mock_1801", "Col_mock_1802", "Col_mock_1803",
                 "Col_nlp20_030", "Col_nlp20_0301", "Col_nlp20_0302", "Col_nlp20_0303", 
                 "Col_nlp20_180", "Col_nlp20_1801", "Col_nlp20_1803")                 # there is no Col_nlp20_1802
keep_factor <- c("Treatment", "Timepoint")
DESeq2_design <- design %>% column_to_rownames(var = "SampleID")

DESeq2_design <- DESeq2_design[keep_sample, ] %>% select(all_of(keep_factor))
dds <- DESeqDataSetFromMatrix(countData = raw_counts[, keep_sample],
                              colData = DESeq2_design,
                              design = ~ Treatment + Timepoint)

keep <- rowSums(counts(dds) >= 5) >= 3 #genes with more than 5 counts in more than 2 samples

dds <- dds[keep,]
transformed_counts <- vst(dds, blind = FALSE)


### PCA plot

for (exp_factor in c("Treatment")){
  pcaData <- plotPCA(transformed_counts, intgroup = keep_factor, 
                     ntop = nrow(transformed_counts), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  g <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = group))
  g <- g + geom_point(size = 4, alpha = 1) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_classic() + 
    ggtitle("PCA plot, nlp20 vs. mock treatment") 
}
g


### DGEA

res <- results(dds, contrast = c("Treatment", "nlp20", "mock")) # base for the comparison at the end
overview_res <- data.frame(GeneID = row.names(res), log2FoldChange = res$log2FoldChange, p_adj = res$padj)
# thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1
count_up_deseq2 <- overview_res %>% filter(log2FoldChange >= lfc_cutoff & p_adj <= padj_cutoff) %>% nrow()
count_down_deseq2 <- overview_res %>% filter(log2FoldChange <= -lfc_cutoff & p_adj <= padj_cutoff) %>% nrow()
count_total_deseq2 <- overview_res %>% nrow()
