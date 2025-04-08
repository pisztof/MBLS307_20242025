library(pheatmap)

df <- data.frame(row.names = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"),
                 condition1 = runif(10, min = -3, max = 3),
                 condition2 = runif(10, min = -3, max = 3),
                 condition3 = runif(10, min = -3, max = 3),
                 condition5 = runif(10, min = -3, max = 3),
                 condition5 = runif(10, min = -3, max = 3))

pheatmap(df, border_color = NA, cluster_cols=FALSE, filename="Plots/heatmap.png")
?pheatmap