install.packages("pheatmap")
library(pheatmap)
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

#exercisegenes
# Load the necessary library
library(pheatmap)

# Generate data for 10 genes across 5 conditions
set.seed(123)
df <- data.frame(
  row.names = paste0("gene", LETTERS[1:10]),
  sample1 = runif(10, min = -3, max = 3),
  sample2 = runif(10, min = -3, max = 3),
  sample3 = runif(10, min = -3, max = 3),
  sample4 = runif(10, min = -3, max = 3),
  sample5 = runif(10, min = -3, max = 3)
)

# Save the heatmap as a PNG file
png("heatmap2.png")

# Create the heatmap without the grey borders and clustering columns
pheatmap(df, 
         cluster_cols = FALSE, 
         border_color = NA)  # Remove grey border lines
dev.off()


