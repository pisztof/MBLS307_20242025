days <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
daysf <- factor(days, levels=c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
daysf

days2 <- c(1,5,2,6,3,7,4,2,3,4)
days2f <- factor(days2)
levels(days2f) <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
days2f

install.packages('pheatmap', dependencies = TRUE)
library(pheatmap)


samples <- data.frame(row.names = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG", "geneH", "geneI", "geneJ"),
                 condition1 = runif(10, min = -3, max = 3),
                 condition2 = runif(10, min = -3, max = 3),
                 condition3 = runif(10, min = -3, max = 3),
                 condition4 = runif(10, min = -3, max = 3),
                 condition5 = runif(10, min = -3, max = 3))

pheatmap(samples,
         cluster_cols = FALSE,
         border_color = NA)
