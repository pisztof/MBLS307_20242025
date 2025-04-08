library(ggplot2)
library(tidyverse)

#Piping exercise
mpg %>% pull(displ) %>% mean()

#Data tidying
#One approach for table 2 would be to have two lists, cases and population, and iterate through each row and store the variable in the "count" column into the corresponding list based on what is in the "type" column. Then the rate could be calculated through vector operations, and this new vector could be repeated (A,B,C) -> (A,A,B,B,C,C) and appended as a new column
#For table 3, we could get the rate column, split based on "/", convert the 2 resulting strings to numeric, and then perform the calculation. This could then be saved as a new column (all of these steps could be performed at once in mutate)
