install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(tidyr)
mpg_df <- mpg

table1 %>% mutate(rate=cases/population*10000)
table1 %>% group_by(year) %>% summarize(total_cases=sum(cases))

ggplot(table1, aes(x=year, y=cases)) + geom_line(aes(group = country), color="grey50")+geom_point(aes(color=country, shape=country)) + scale_x_continuous(breaks = c(1999, 2000))

#1
#the tables show the number of documented cases of TB
table1
#table 1 shows country, year, cases per 10000 and total population
table2
#table 2 shows country, year, type and count
table3
#table 3 shows country, year and rate of cases

#2 calc rate for table2 and table3
#a.extract TB cases/country/year
#b.extract matching popu/country/year
#c.divide cases by popu, x10000
#d.store back

#table2 was in a pivot_long mode and we made it wider
table2_wide <- table2 %>% pivot_wider(names_from = type, values_from = count)
table2_wide <- table2_wide %>% mutate(rate=(cases/population)*10000)
table2_wide

#table 3 has the rate but uncalculated and as a string
#first we make 2 separate columns with the data of the cases and the population, using the "/" to distinguish between what data should go yo what column
table3_int <- table3 %>% separate(rate, into = c("cases", "population"), sep="/", convert = TRUE)
table3_int #now we have 2 separate columns + the data is in integer format

table3_int <- table3_int %>% mutate(rate = (cases/population)*10000)
table3_int
