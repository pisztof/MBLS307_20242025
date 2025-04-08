#heatmappps
install.packages("pheatmap")
library(pheatmap)
df <- data.frame(row.names = c("geneA", "geneB", "geneC"),
                 sample1 = runif(3, min = -3, max = 3),
                 sample2 = runif(3, min = -3, max = 3),
                 sample3 = runif(3, min = -3, max = 3),
                 sample4 = runif(3, min = -3, max = 3))
pheatmap(df)
# re-draw the heatmap without clustering genes
pheatmap(df,
         cluster_cols = FALSE)

#excercise (making a heatmap)
# Using the above code as a basis, make a heatmap for the expression of 10 genes across 5 conditions
df_2 <- data.frame(row.names = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG", "geneH", "geneI", "geneL"),
                   condition1 = runif(10, min =-3, max =3),
                   condition2 = runif(10, min =-3, max =3),
                   condition3 = runif(10, min =-3, max =3),
                   condition4 = runif(10, min =-3, max =3),
                   condition5 = runif(10, min =-3, max =3))
## number of values to be drawn by runif() should be the same as the number of rows (genes)
?pheatmap
pheatmap(df_2,
         cluster_cols = FALSE, border_color = NA)                   
