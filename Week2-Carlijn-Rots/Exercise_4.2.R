# Exercise
# The mpg dataframe (see example in Chapter 2 “Making beautiful plots”) 
# has information about the engine size (column “displ”) 
# and the class of the car (“class”). Using the ggplot function, 
# make a boxplot displaying engine size per car class.
# Add individual datapoints on top of the boxplot and 
# change the theme to “classic”. 
# The expected plot is shown in Figure 3.

install.packages("tidyverse")
library(tidyverse)
mpg_df <- mpg
str(mpg_df)
ggplot(mpg_df) +
  geom_boxplot(aes(y = displ, x = class)) +
  geom_jitter(aes(y = displ, x = class), colour = "black", height = 0.1, width = 0.2) + 
  xlab("Class of vehicle") +
  ylab("Engine size (L)") +
  ggtitle(label = "Boxplot on ex 4.2") +
  theme_bw() +
  theme(panel.grid = element_blank())
