library("clusterProfiler")
suppressPackageStartupMessages(library("org.At.tair.db"))
library("biomartr")
library("biomaRt")
library("tidyverse")
library("ggpubr")

# get data
overview_res <- read_tsv("./results-deseq2/summary_dgea_deseq2_nlp20.tsv")
upregulated <- overview_res %>% dplyr::filter(log2FoldChange >= 1 & p_adj <= 0.05) %>% dplyr::select(GeneID) %>% pull()
downregulated <- overview_res %>% dplyr::filter(log2FoldChange <= -1 & p_adj <= 0.05) %>% dplyr::select(GeneID) %>% pull()
all_expressed <- overview_res %>% dplyr::select(GeneID) %>% pull()


### KEGG analysis 

# KEGG analysis - upregulated
kegg_upreg <- enrichKEGG(gene = upregulated,
                            universe = all_expressed,
                            organism = "ath",
                            keyType = "kegg",
                            minGSSize = 10,
                            maxGSSize = 500,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            use_internal_data = FALSE)
plot_kegg_upreg <- dotplot(kegg_upreg,
        color = "qvalue",
        showCategory = 10,
        size = "Count",
        title = "KEGG term enrichment analysis for upregulated genes")

# KEGG analysis - downregulated
kegg_downreg <- enrichKEGG(gene = downregulated,
                            universe = all_expressed,
                            organism = "ath",
                            keyType = "kegg",
                            minGSSize = 10,
                            maxGSSize = 500,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            use_internal_data = FALSE)
plot_kegg_downreg <- dotplot(kegg_downreg,
        color = "qvalue",
        showCategory = 10,
        size = "Count",
        title = "KEGG term enrichment analysis for downregulated genes")

# Plot KEGG dotplots
ggarrange(plot_kegg_upreg, plot_kegg_downreg, ncol=1, nrow=2, common.legend=TRUE, legend="right")


### GO analysis 

# Get ENTREZ
upregulated_entrez <- base::unique(select(org.At.tair.db, keys = upregulated,
                                          column = c('ENTREZID', 'SYMBOL'), keytype = 'TAIR')[,'ENTREZID'])
upregulated_entrez <- upregulated_entrez[!is.na(upregulated_entrez)]

downregulated_entrez <- base::unique(select(org.At.tair.db, keys = downregulated,
                                            column = c('ENTREZID', 'SYMBOL'), keytype = 'TAIR')[,'ENTREZID'])  # ignore returned 1:many mapping between keys and columns
downregulated_entrez <- downregulated_entrez[!is.na(downregulated_entrez)]

all_expressed_entrez <- base::unique(select(org.At.tair.db, keys = all_expressed,
                                            column = c('ENTREZID', 'SYMBOL'), keytype = 'TAIR')[,'ENTREZID'])  # ignore returned 1:many mapping between keys and columns
all_expressed_entrez <- all_expressed_entrez[!is.na(all_expressed_entrez)]

# GO analysis - upregulated

go_upreg <- data.frame(enrichGO(gene = upregulated_entrez,
                                    universe = all_expressed_entrez,
                                    OrgDb = org.At.tair.db,  
                                    keyType = "ENTREZID",
                                    ont = "BP",              
                                    pAdjustMethod = "BH",
                                    qvalueCutoff = 0.05,
                                    readable = TRUE,
                                    pool = FALSE))

go_upreg_ordered <- go_upreg[order(-go_upreg$Count), ]
go_upreg_top15 <- head(go_upreg_ordered, 15)

plot_go_upreg <- ggplot(go_upreg_top15, aes(x=Count, y=reorder(Description, p.adjust), fill=p.adjust)) +
  theme_classic() +
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y='', title = 'GO term enrichment analysis for upregulated genes')

# GO analysis - downregulated

go_downreg <- data.frame(enrichGO(gene = downregulated_entrez,
                                    universe = all_expressed_entrez, 
                                    OrgDb = org.At.tair.db,  
                                    keyType = "ENTREZID",
                                    ont = "BP",              
                                    pAdjustMethod = "BH",
                                    qvalueCutoff = 0.05,
                                    readable = TRUE, 
                                    pool = FALSE))

go_downreg_ordered <- go_downreg[order(-go_downreg$Count), ]
go_downreg_top15 <- head(go_downreg_ordered, 15)

plot_go_downreg <- ggplot(go_downreg_top15, aes(x=Count, y=reorder(Description, p.adjust), fill=p.adjust)) +
  theme_classic() +
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y='', title = 'GO term enrichment analysis for downregulated genes')

# Plot GO dotplots
ggarrange(plot_go_upreg, plot_go_downreg, ncol=1, nrow=2, align="v", common.legend=TRUE, legend="right")

