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
pheatmap(df, cluster_cols = FALSE)



df2 <- data.frame(row.names = paste('gene',letters[1:10],sep=''), matrix(runif(50,min=-3,max=3), nrow=10,ncol=5))
colnames(df2) <- paste('condition', 1:ncol(df2), sep='')

pheatmap(df2, cluster_cols=FALSE, border_color = NA, filename = 'rheatmap.png')

