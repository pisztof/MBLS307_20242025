install.packages("ggplot2")
library(ggplot2)
### exercise - basic use of ggplot2 for data visualization
install.packages("tidyverse")
# load libraries
library("tidyverse")

# we use the mpg dataset available by default in R
##car class is in the column "class",
##size of the engine (L) is in the column "displ".
ggplot(mpg, aes(x = class, y = displ)) +
  geom_boxplot() +
  geom_jitter() + ## induces a slight random spread of the points, helpful  otherwise would be clustered
  theme_classic()
    
ggplot(mpg, aes(x = class, y = displ)) +
  geom_boxplot() +
  geom_point() + #shows lesser data points
  theme_classic() #theme is set to classic
