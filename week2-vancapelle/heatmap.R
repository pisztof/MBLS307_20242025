# 14-2-2025
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

df2 <- data.frame(row.names = c("geneA", "geneB", "geneC", "geneD", "geneE"),
                 sample1 = runif(5, min = -3, max = 3),
                 sample2 = runif(5, min = -3, max = 3),
                 sample3 = runif(5, min = -3, max = 3),
                 sample4 = runif(5, min = -3, max = 3),
                 sample5 = runif(5, min = -3, max = 3))
pheatmap(df2,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA)
png("output/heatmap.png")
pheatmap(df2,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA)
dev.off()
