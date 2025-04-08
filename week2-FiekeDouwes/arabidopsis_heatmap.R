# Heatmap exercise from intermediate trajectory
# Fri Feb 14 15:04:41 2025 ------------------------------


### make a data frame with gene expression values. The gene IDs should be in row names.
df <- data.frame(row.names = c("geneA", "geneB", "geneC"),
                 sample1 = runif(3, min = -3, max = 3),
                 sample2 = runif(3, min = -3, max = 3),
                 sample3 = runif(3, min = -3, max = 3),
                 sample4 = runif(3, min = -3, max = 3))

### make a simple heatmap
# load library
library(pheatmap) 

# make the heatmap
pheatmap(df)

# re-draw the heatmap without clustering genes
pheatmap(df,
         cluster_cols = FALSE)

ten_genes <- data.frame(row.names = c("gene1","gene2","gene3","gene4","gene5",
                                      "gene6","gene7","gene8","gene9","gene10"),
                        condition1 = runif(10, min = -3, max = 3),
                        condition2 = runif(10, min = -3, max = 3),
                        condition3 = runif(10, min = -3, max = 3),
                        condition4 = runif(10, min = -3, max = 3),
                        condition5 = runif(10, min = -3, max = 3))
png("output/heatmap_exercise.png")
pheatmap(ten_genes, border_color = NA, cluster_cols = FALSE)
dev.off()
