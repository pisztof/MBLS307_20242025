install.packages("tidyverse")
library(tidyverse)
mpg_df <- mpg
#Using the ggplot function
#make a boxplot displaying 
#engine size per car class
#Add individual datapoints on top 
#of the boxplot and
#change the theme to “classic”.

mpg_df %>% ggplot(aes(x=class, y=displ)) + geom_point() + geom_boxplot()+ theme_classic()
#i can not manage to make the scatterplot and the boxplot to be plotted independently of each other on the same graph

ggplot(mpg, aes(x = class, y = displ)) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()
#i was missing geom_jitter
