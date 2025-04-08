#Exercises from chapter 5 of R for Data Science
#https://r4ds.hadley.nz/data-tidy
# %>% is sometimes shown as |>

# Fri Feb 14 11:35:45 2025 ------------------------------

library(tidyverse)

#5.2.1
#1
#table1: observations are specific countries in specific years.
# First two columns are the country and year of the observation, and last two
# Last two columns are cases and size of population in country in that year.
#table2: cases and population are merged into 1 column 
# and split into 2 observations.
#table3: cases and population are merged into one column as "rate"

#2
#table1: extract cases & population columns and perform mathematical operation.
# Place back as rate column.
#table2: extract cases & population using groupby(type)
#table3: rate already given

# dataset %>% 
#   pivot_longer(
#   cols = col:col or starts_with("str") etc,
#   names_to = "name of new column" 
#     or c("names","of","new","columns") if multiple vars included in column names
#     or c(".value","num") if values and variable names included in column names (eg. dob_child1 and name_child1 give (".value","child"))
#   names_sep = "sep" in case multiple vars in column names (eg. sp_m_014)
#   values_to = "value"
#   values_drop_na = TRUE if there are NA's
 
# dataset %>% 
#   pivot_wider(
#   id_cols = columns with identifiers. Otherwise, tibble still isn't shortened.
#   names_from = column_with_names
#   values_from = column_with_values
#   )