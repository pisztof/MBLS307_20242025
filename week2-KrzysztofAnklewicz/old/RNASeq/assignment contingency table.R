library("tidyverse")
library("UpSetR")

# Contingency table flg22 vs. nlp20
overview_res_nlp20 <- read_tsv("./results-deseq2/summary_dgea_deseq2_nlp20.tsv")
up_nlp20 <- overview_res_nlp20 %>% dplyr::filter(log2FoldChange >= 1 & p_adj <= 0.05) %>% dplyr::select(GeneID) %>% pull()
overview_res_nlp20 <- overview_res_nlp20 %>% dplyr::select(GeneID) %>% pull()

overview_res_flg22 <- read_tsv("./results-deseq2/summary_dgea_deseq2_flg22.tsv")
up_flg22 <- overview_res %>% dplyr::filter(log2FoldChange >= 1 & p_adj <= 0.05) %>% dplyr::select(GeneID) %>% pull()
overview_res_flg22 <- overview_res_flg22 %>% dplyr::select(GeneID) %>% pull()

# calculating the number of genes common to both analyses to use as universe size
universe <- list(nlp20 = overview_res_nlp20, flg22 = overview_res_flg22) 
upset(fromList(universe), order.by = "freq") # 20485 common genes

sets <-list(up_flg22 = up_flg22, up_nlp20 = up_nlp20)
upset(fromList(sets), order.by = "freq") # Overlap: 1297, flg22: 818, nlp20: 39
20485 - 1297 - 818 - 39 # = 18331 <--- "background" genes, not induced by either

cont_table <- data.frame(among_flg22_up = c(1297, 818), not_among_flg22_up = c(39, 18331),
                         row.names = c("among_nlp20_up", "not_among_nlp20_up"))
print(cont_table)
fisher.test(cont_table)

# using the total number of aradopsis genes
27562 - 1297 - 818 - 39 # = 25408
cont_table <- data.frame(among_flg22_up = c(1297, 818), not_among_flg22_up = c(39, 25408),
                         row.names = c("among_nlp20_up", "not_among_nlp20_up"))
print(cont_table)
fisher.test(cont_table)

# Contingency table nlp20 vs. WRKY targets

targets <- read_tsv("../data/BB-RNAseq-input-datasets/WRKY-targets-Birkenbihl-et-al-2018.tsv")
sets2 <-list(wrky_target = targets$GeneID, up_nlp20 = up_nlp20)
upset(fromList(sets2), order.by = "freq")

# Overlap: 662, nlp20: 1453, WRKY: 1288
27562 - 662 - 1453 - 1288 # <--- background: 24159
cont_table2 <- data.frame(among_nlp20_up = c(662, 1453), not_among_nlp20_up = c(1288, 24159),
                         row.names = c("among_wrky_up", "not_among_wrky_up"))
print(cont_table2)
fisher.test(cont_table2)
