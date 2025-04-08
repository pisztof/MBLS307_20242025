#Exercise about pipes

#Rewrite the code to calculate the average engine size in the dataset 
#Auto-MPG mean(mpg$displ). For that, use pipes and tidyverse functions 
#select() and pull() to get the engine size values.

# Fri Feb 14 11:23:48 2025 ------------------------------

library(tidyverse)
mpg_df <- mpg
mean(mpg_df$displ)
mpg_df %>% pull(displ) %>% mean()
mpg_df %>% select(displ) %>% mean()

