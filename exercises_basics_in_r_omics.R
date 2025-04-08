#exercise 2
log(12.43)
log10(12.43)
log2(12.43)
sqrt(12.43)
exp(12.43)
area_circle <- pi*20^2
area_circle
(14*0.51)^1/3
weight <- c(69,62,57,59,59,64,56,66,67,66)
mean(weight)
var(weight)
sd(weight)
range(weight)
length(weight)
first_five <- weight[1:5]
first_five
height <- c(112,102,83,84,99,90,77,112,133,112)
summary(height)
some_child <- height[c(2,3,9,10)]
some_child
shorter_child <- height[height <= 99]
shorter_child
bmi <- weight/((height/100)^2)
bmi
seq(from = 0, to = 1, by = 0.1)
seq1 <- seq(from = 0, to = 1, by = 0.1)
seq1
seq2 <- rev(seq(from = 1, to = 10, by = 0.5))
seq2
seq3 <- rep(1:3, times = 3)
seq3
seq4 <- rep(c("a","c", "e", "g"), times =3)
seq4
seq5 <- rep(c("a","c", "e", "g"), each =3)
seq5
seq6 <- rep(1:3, each =3, times = 2)
seq6
seq7 <- rep(1:5,times = c(5,4,3,2,1))
seq7
seq8 <- rep(c(7,2,8,1), times = c(4,3,1,5))
seq8
height_sorted <- sort(height)
height_sorted
height_sorted2 <- sort(height, decreasing =TRUE)
height_sorted2
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")
height_ord <- order(height)
height_ord
names_sort <- child_name[height_ord]
names_sort
#louise is shortest, ronald is tallest
weight_ord <- rev(order(weight))
weight_ord
weight_rev <- child_name[weight_ord]
weight_rev
#alfred is heaviest, louise is lightest
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mean(mydata)
mean(mydata, na.rm = TRUE)
help(mean)
rm(seq1)
#chapter 3
flowers <- read.table(file="C:/Users/mloui/Downloads/all_data/data/flower.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
str(flowers)
# 'data.frame':	96 obs. of  8 variables:
#   $ treat    : Factor w/ 2 levels "notip","tip": 2 2 2 2 2 2 2 2 2 2 ...
# $ nitrogen : Factor w/ 3 levels "high","low","medium": 3 3 3 3 3 3 3 3 3 3 ...
# $ block    : int  1 1 1 1 1 1 1 1 2 2 ...
# $ height   : num  7.5 10.7 11.2 10.4 10.4 9.8 6.9 9.4 10.4 12.3 ...
# $ weight   : num  7.62 12.14 12.76 8.78 13.58 ...
# $ leafarea : num  11.7 14.1 7.1 11.9 14.5 12.2 13.2 14 10.5 16.1 ...
# $ shootarea: num  31.9 46 66.7 20.3 26.9 72.7 43.1 28.5 57.8 36.9 ...
# $ flowers  : int  1 10 10 1 4 9 7 6 5 8 ...
names(flowers)
flowers$height
f_height <- flowers$height
mean(f_height)
summary(f_height)
mean(flowers$height)
summary(flowers$height)
flowers[1,4]
flowers$height[1]
flowers[1:10, 1:4]
flowers[c(1,5,12,30),c(1,3,6,8)]
flowers[1:8,]
big_flowers <- flowers[flowers$height > 12,]
big_flowers
tim_med_2 <- subset(flowers, treat =="tip"& nitrogen == "medium" & block ==2)
tim_med_2
height_ord <- flowers[order(flowers$height),]
height_ord
#exercise 3 intro2r
whale <- read.table(file="C:/Users/mloui/Downloads/all_data/data/whaledata.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
whale
head(whale)
names(whale)
str(whale)
# 'data.frame':	100 obs. of  8 variables:
#   $ month          : Factor w/ 2 levels "May","October": 1 1 1 1 1 1 1 1 1 1 ...
# $ time.at.station: int  1344 1633 743 1050 1764 580 459 561 709 690 ...
# $ water.noise    : Factor w/ 3 levels "high","low","medium": 2 3 3 3 3 1 3 3 3 3 ...
# $ number.whales  : int  7 13 12 10 12 10 5 8 11 12 ...
# $ latitude       : num  60.4 60.4 60.5 60.3 60.4 ...
# $ longitude      : num  -4.18 -4.19 -4.62 -4.35 -5.2 -5.22 -5.08 -5 -4.64 -4.84 ...
# $ depth          : int  520 559 1006 540 1000 1000 993 988 954 984 ...
# $ gradient       : int  415 405 88 409 97 173 162 162 245 161 ...

#the dataset has 100 observations and 8 variables, month and water.noise are factors
summary(whale)
#number.whales has one NA
whale.sub <- whale[1:10,1:4]
whale.num <- whale[,c("month", "water.noise", "number.whales")]
whale.may <- whale[1:50,]
whaleminrow <- whale[-(1:10),-(8)]
whaleminrow2 <- whale[-(1:10),"gradient"]
whale_deep <- whale[whale$depth>1200,]
whale_steep <- whale[whale$gradient>200,]
str(whale_steep)
whale_quiet <- whale[whale$water.noise == "low",]
whale_loud_may <- whale[whale$water.noise == "high" & whale$month == "May",]
whale_132 <- whale[whale$month == "October"& whale$water.noise == "low"& whale$gradient>132,]
whale_lat <- whale[whale$latitude>60.0 & whale$latitude<61.0 & whale$longitude<(-4.0) & whale$longitude>(-6.0),]
whale_noise <- whale[whale$water.noise != "medium",]
whale_median <- whale[whale$month == "October"& whale$water.noise == "low"& whale$gradient>median(whale$gradient),]
whale_median == whale_132
whale_great <- whale[whale$depth>1500 & whale$number.whales>mean(whale$number.whales),]
str(whale_great)
mean_whales <- mean(whale$number.whales, na.rm = TRUE)
#remove na (na.rm) for calculating mean to fix problem with output
mean_whales
whale_great <- whale[whale$depth>1500 & whale$number.whales>mean(whale$number.whales, na.rm = TRUE),]
whale_great
str(whale_great)
#12
whale_sub <- subset(whale, month == "May" & time.at.station<1000 & depth>1000)
whale_sub2 <- subset(whale, month == "October" & latitude>61, select = c("month", "latitude", "longitude", "number.whales"))
str(whale_sub2)
#13
whale.depth.sort <- whale[order(whale$depth),]
whale_sort2 <- whale[order(-xtfrm(whale$water.noise), whale$depth),]
str(whale_sort2)
whale_sort_rev <- whale[order(-xtfrm(whale$water.noise), -whale$depth),]
str(whale_sort_rev)
whale_sort_rev
mean(whale$time.at.station)
whale_mean_noise <- tapply(whale$number.whales, whale$water.noise, mean, na.rm = TRUE)
whale_mean_noise
whale_mean_noise_month <- tapply(whale$number.whales, list(whale$water.noise, whale$month), mean, na.rm = TRUE)
whale_mean_noise_month
?aggregate
whale_aggregate <- aggregate(whale[c("time.at.station","number.whales", "depth", "gradient")], by = list(water.noise = whale$water.noise), FUN = mean, na.rm = TRUE)
whale_aggregate
whale_aggregate2 <- aggregate(whale[c("time.at.station","number.whales", "depth", "gradient")], by = list(water.noise = whale$water.noise, month = whale$month), FUN = mean, na.rm = TRUE)
whale_aggregate2
#try to order levels from low to high
whale$water.noise <- factor(whale$water.noise, levels = c("low", "medium", "high"))
whale_aggregate3 <- aggregate(whale[c("time.at.station","number.whales", "depth", "gradient")], by = list(water.noise = whale$water.noise, month = whale$month), FUN = mean, na.rm = TRUE)
whale_aggregate3
str(whale_aggregate3)
table(whale$water.noise)
table(whale$water.noise, whale$month)
write.table(whale.num, file = "C:/Users/mloui/Downloads/all_data/data/whale_num.txt", col.names = TRUE, row.names = FALSE, sep = "\t")

#exercise 4 intro2r
squid <- read.table(file="C:/Users/mloui/Downloads/all_data/data/squid1.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
str(squid)
summary(squid)
#squid dataframe has 519 observations and 13 variables
squid$year_factor <- factor(squid$year)
squid$month_factor <- factor(squid$month)
squid$maturity.stage_factor <- factor(squid$maturity.stage)
str(squid)
table(squid$month_factor, squid$year_factor)
#we do not have data for each month each year, 1990 and 1991 have the most observations
dotchart(squid$DML)
dotchart(squid$weight)
dotchart(squid$nid.length)
dotchart(squid$ovary.weight)
#there is an unusually large observation for the nid length
par(mfrow = c(2, 2))
dotchart(squid$DML)
dotchart(squid$weight)
dotchart(squid$nid.length)
dotchart(squid$ovary.weight)
par(mfrow = c(1, 1))

pdf(file="C:/Users/mloui/Downloads/all_data/data/squid_plots.pdf")
par(mfrow = c(2, 2))
dotchart(squid$DML)
dotchart(squid$weight)
dotchart(squid$nid.length)
dotchart(squid$ovary.weight)
dev.off()

which(squid$nid.length > 400)
squid$nid.length[11]
squid$nid.length[11] <- 43.2
dotchart(squid$nid.length)

par(mfrow = c(2, 2))
hist(squid$DML)
hist(squid$weight)
hist(squid$nid.length)
hist(squid$ovary.weight)

pdf(file="C:/Users/mloui/Downloads/all_data/data/squid_histograms.pdf")
par(mfrow = c(2, 2))
hist(squid$DML)
hist(squid$weight)
hist(squid$nid.length)
hist(squid$ovary.weight)
dev.off()

par(mfrow = c(1, 1))

hist(squid$DML)
max(squid$DML)
min(squid$DML)

brk <- seq(from = 88, to = 323, by = 1)
hist(squid$DML, breaks = brk, main = "squid DML")

plot(x=squid$DML, y=squid$weight)
plot(x=squid$DML, y=log(squid$weight))
plot(x=squid$DML, y=sqrt(squid$weight))

squid$weight_log <- log(squid$weight)
squid$weight_sqrt <- sqrt(squid$weight)

pdf(file="C:/Users/mloui/Downloads/all_data/data/squid_scatterplots.pdf")
par(mfrow = c(2, 2))
plot(x=squid$DML, y=squid$weight)
plot(x=squid$DML, y=squid$weight_log)
plot(x=squid$DML, y=squid$weight_sqrt)
dev.off()
#sqrt transformation linearizes the data best

str(squid)

# 'data.frame':	519 obs. of  18 variables:
#   $ sample.no            : int  105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 ...
# $ specimen             : int  1002 1003 1005 1007 1008 1009 1011 1013 1014 1017 ...
# $ year                 : int  1989 1989 1989 1989 1989 1989 1989 1989 1989 1989 ...
# $ month                : int  12 12 12 12 12 12 12 12 12 12 ...
# $ weight               : num  152 106 138 141 126 ...
# $ sex                  : int  2 2 2 2 2 2 2 2 2 2 ...
# $ maturity.stage       : int  3 1 2 2 3 1 2 3 3 4 ...
# $ DML                  : int  174 153 169 175 169 116 135 192 170 205 ...
# $ eviscerate.weight    : num  87.5 62.6 79.4 83.1 72.2 ...
# $ dig.weight           : num  4.648 3.138 0.307 4.123 3.605 ...
# $ nid.length           : num  39.4 24.1 39 41.4 39.8 20 14 55 44 53 ...
# $ nid.weight           : num  2.46 0.319 1.169 1.631 2.03 ...
# $ ovary.weight         : num  1.68 0.103 0.289 0.252 0.86 ...
# $ year_factor          : Factor w/ 3 levels "1989","1990",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ month_factor         : Factor w/ 12 levels "1","2","3","4",..: 12 12 12 12 12 12 12 12 12 12 ...
# $ maturity.stage_factor: Factor w/ 5 levels "1","2","3","4",..: 3 1 2 2 3 1 2 3 3 4 ...
# $ weight_log           : num  5.02 4.66 4.93 4.95 4.84 ...
# $ weight_sqrt          : num  12.3 10.3 11.8 11.9 11.2 ...

boxplot(DML~maturity.stage_factor,data = squid, ylab = "dorsla mantle length (DML)", xlab = "maturity stage")

#exercise making a boxplot wit ggplot2

install.packages("tidyverse")
library(tidyverse)
mpg_df <- mpg
head(mpg_df)
str(mpg_df)

ggplot(mpg_df, aes(x= class, y= displ)) + geom_boxplot() + geom_jitter() + theme_classic()

#pipes example
mpg_df %>% head()

#chapter 5 tidyverse

table1
#compute rate per 10000
table1 %>% mutate(rate = cases / population * 10000)
#compute total cases per year
table1 %>% group_by(year) %>% summarize(total_cases = sum(cases))

ggplot(table1, aes(x = year, y = cases)) + 
  geom_line(aes(group = country), color = "grey50") + 
  geom_point(aes(color = country, shape = country)) +
  scale_x_continuous(breaks = c(1999, 2000))

#chapter 5 exercises

#in table 1 each column is a variable and each row is an observation split into countries and years
#in table 2 there is a row for each observation of cases and population, and into countries and years
# in table 3 there are two observations per cell, rows split into country and year but rate is combined into one

#calculate rate based on counts for table 2, use rate and get result for table 3

billboard

#tidy billboard data by pivoting
billboard %>% 
  pivot_longer(
    cols = starts_with("wk"),
    names_to = "week",
    values_to = "rank"
  )

#get rid of na values
billboard %>% 
  pivot_longer(
    cols = starts_with("wk"),
    names_to = "week",
    values_to = "rank",
    values_drop_na = TRUE
  )

#convert week values from character strings into numbers
billboard_longer<- billboard %>% 
  pivot_longer(
    cols = starts_with("wk"),
    names_to = "week",
    values_to = "rank",
    values_drop_na = TRUE
  ) %>% 
  mutate(
    week = parse_number(week)
  )
billboard_longer

#visualize data using ggplot
billboard_longer %>% 
  ggplot(aes(x = week, y = rank, group = track)) +
  geom_line(alpha = 0.25) +
  scale_y_reverse()

#example of pivoting
df <- tribble(
  ~id,  ~bp1, ~bp2,
  "A",  100,  120,
  "B",  140,  115,
  "C",  120,  125
)

#pivot df longer
df |> 
  pivot_longer(
    cols = bp1:bp2,
    names_to = "measurement",
    values_to = "value"
  )

who2
who2 |> 
  pivot_longer(
    cols = !(country:year),
    names_to = c("diagnosis", "gender", "age"), 
    names_sep = "_",
    values_to = "count"
  )

household
household |> 
  pivot_longer(
    cols = !family, 
    names_to = c(".value", "child"), 
    names_sep = "_", 
    values_drop_na = TRUE
  )

cms_patient_experience
cms_patient_experience |> 
  distinct(measure_cd, measure_title)

cms_patient_experience |> 
  pivot_wider(
    names_from = measure_cd,
    values_from = prf_rate
  )

cms_patient_experience |> 
  pivot_wider(
    id_cols = starts_with("org"),
    names_from = measure_cd,
    values_from = prf_rate
  )

df <- tribble(
  ~id, ~measurement, ~value,
  "A",        "bp1",    100,
  "B",        "bp1",    140,
  "B",        "bp2",    115, 
  "A",        "bp2",    120,
  "A",        "bp3",    105
)

df |> 
  pivot_wider(
    names_from = measurement,
    values_from = value
  )

df |> 
  distinct(measurement) |> 
  pull()

df |> 
  select(-measurement, -value) |> 
  distinct()

df |> 
  select(-measurement, -value) |> 
  distinct() |> 
  mutate(x = NA, y = NA, z = NA)

df <- tribble(
  ~id, ~measurement, ~value,
  "A",        "bp1",    100,
  "A",        "bp1",    102,
  "A",        "bp2",    120,
  "B",        "bp1",    140, 
  "B",        "bp2",    115
)

df |>
  pivot_wider(
    names_from = measurement,
    values_from = value
  )

df |> 
  group_by(id, measurement) |> 
  summarize(n = n(), .groups = "drop") |> 
  filter(n > 1)

#exercise 6

#make function to calculate area of circle
area_circle <- function(r){pi*(r^2)}
area_fun <- area_circle(3.4)
#test if function works
area_manually <- pi*(3.4^2)
area_manually
identical(area_fun,area_manually)
vec <- c(3.5,4.2,1.2,1.2,4.8)
area_circle(vec)

#make a function to convert fahrenheit to centigrade
f_to_c <- function(fahrenheit){
  centigrade <- (fahrenheit-32)*5/9
  cat("Fahrenheit:", fahrenheit , "is equivalent to ",centigrade, "centigrade.")
}
f_to_c (270)

#make function to calculate mean, median, range of vector
vec <- rnorm(100, mean=35, sd=15)
vec_fun <- function(i){
  cat("mean:", mean(i), "median:", median(i), "range:", range (i))
  hist(i, freq = FALSE)
  dens <- density(i)
  lines(dens)
}
vec_fun(vec)

#make function to calculate median of vectors
median_fun <- function(vec){
  if (length(vec) %% 2 == 0) {
    (vec[length(vec)/2] + vec[length(vec)/2+1])/2
  } else {
    vec[length(vec)/2+1]
  }
}
#make vectors to test function with and use them
vec_even <- c(1,2,3,4,5,6,7,8)
vec_odd <- c(1,2,3,4,5,6,7)
median_fun(vec_even)
median_fun(vec_odd)

#test if function works correctly
identical(median(vec_odd), median_fun(vec_odd))
identical(median(vec_even), median_fun(vec_even))

#make function to investigate ricker model
ricker_model <- function(nzero, r, time, K = 100){
  N <- numeric(time + 1)
  N[1] <- nzero
  for (t in 1:time){
    N[t+1] <- N[t]*exp(r*(1-N[t]/K))
  }
  return(N)
}
#test function
ricker_model(20,0.5,30,100)

#example making a heatmap
### make a data frame with gene expression values. The gene IDs should be in row names.
df <- data.frame(row.names = c("geneA", "geneB", "geneC"),
                 sample1 = runif(3, min = -3, max = 3),
                 sample2 = runif(3, min = -3, max = 3),
                 sample3 = runif(3, min = -3, max = 3),
                 sample4 = runif(3, min = -3, max = 3))

### make a simple heatmap
# load library
install.packages("pheatmap")
library(pheatmap) 

# make the heatmap
pheatmap(df)

# re-draw the heatmap without clustering genes
pheatmap(df,
         cluster_cols = FALSE)

#exercise making a heatmap

df2 <- data.frame(row.names = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF", "geneG", "geneH", "geneI", "geneJ"),
                 condition1 = runif(10, min = -3, max = 3),
                 condition2 = runif(10, min = -3, max = 3),
                 condition3 = runif(10, min = -3, max = 3),
                 condition4 = runif(10, min = -3, max = 3),
                 condition5 = runif(10, min = -3, max = 3))
pheatmap(df2, border_color = NA)
