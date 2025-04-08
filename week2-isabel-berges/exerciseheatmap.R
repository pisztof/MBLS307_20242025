install.packages("pheatmap")
library(pheatmap)
df <- data.frame(row.names = c("geneA", "geneB", "geneC"),
                 sample1 = runif(3, min = -3, max = 3),
                 sample2 = runif(3, min = -3, max = 3),
                 sample3 = runif(3, min = -3, max = 3),
                 sample4 = runif(3, min = -3, max = 3))

pheatmap(df)

#change color scheme for legend
pheatmap(df, cluster_cols = FALSE)

#exercise
?pheatmap
#with the runif() function you generate 10 random values between the min and the max
#the conditions are the columns and the genes are the rows
df_10 <- data.frame(row.names = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG", "geneH", "geneI", "geneJ"),
                    condition1 = runif(10, min = -3, max = 3),
                    condition2 = runif(10, min = -3, max = 3),
                    condition3 = runif(10, min = -3, max = 3),
                    condition4 = runif(10, min = -3, max = 3),
                    condition5 = runif(10, min = -3, max = 3))

pheatmap(df_10, border_color = NA)

#save as png
png("ten_genes_heatmap.png")
pheatmap(df_10, border_color = NA)
dev.off()
