# Exercise 3
# 14/02/25

#Question 5: Import file into R
## Importing the ‘whaledata.txt’ file into R
whale <- read.table(file = 'data/whaledata.txt', header = TRUE, sep = "\t", stringsAsFactors = TRUE)
# because it is a relative pathway. Refers to the working directory 
# sep = "t\" means separator for tab-delimited filed is a tab (\t) 

## Question 6: to see first few rows of data 
head(whale) ##function to display the first 5 rows of your dataframe.
names(whale) ##which will return a vector of variable names from your dataset.
str(whale) ##function to display the structure of the dataset and a neat summary of your variables.

## How many observations does this dataset have? How many variables are in this dataset? What type of variables are month and water.noise?
## number of observations = 100
## number of variables = 8
##month = factor
## water noise = factor

# Number of observations in this dataset = # rows
nrow(whale)

# Number of variables in this dataset = # columns
ncol(whale)

#To see the type of variables
class(whale$month)
class(whale$water.noise)

##Question 7:  get another useful summary of your dataframe
summary(whale)

##number of missing values (NA in the data)
## 1

## Question 8: Manipulating the dataframe
## Part a -Extract all the elements of the first 10 rows and the first 4 column
whale.sub <- whale[1:10, 1:4]

# Part b
whale.num <- whale[,c("month", "water.noise", "number.whales")]

#Part c
whale.may <- whale [1:50, ]

#Part d
whale.subset <- whale [-(1:10), -ncol(whale)]
# -(1:10) excludes the 1st 10 rows
# -ncol(whale) excludes the last column

## Question 9: Conditional statements to select particular rows based on some logical criteria

#part a - at depths greater than 1200 m
depth_greaterthan_1200 <- whale[whale$depth > 1200, ]

#part b - gradient steeper than 200 degrees
gradient_steeper_than_200 <- whale[whale$gradient > 200, ]

#part c - water noise level of ‘low’
low_waternoise_level <- whale[whale$water.noise == "low", ]

#part d - water.noise level of ‘high’ in the month of ‘May’
may_high_waternoise <- whale[whale$water.noise == "high" & whale$month == "May", ]

#part e - month of ‘October’, water noise level of ‘low’ and gradient greater than the median value of gradient (132)
Oct_levellow_132 <- whale[whale$month == "October" & whale$water.noise == "low" & whale$gradient > 132, ]

#part f - all observations from between latitudes 60.0 and 61.0 and longitudes -6.0 and -4.0
Lat60_61_long6_4 <- whale[whale$latitude >= 60.0 &  whale$latitude <= 61.0 & whale$longitude >= -6.0 & whale$longitude <= -4.0, ]

# part g - all rows that do not have a water noise level of medium
Water_noise_not_medium <- whale[whale$waternoise != "medium", ]

##Question 10: Include R functions
Oct_levellow_132 <- whale[whale$month == "October" & whale$water.noise == "low" & whale$gradient > median(whale$gradient), ]

### Question 11: Problem with output
depth_morethan_1500_spotted_more_average <- whale[whale$depth > 1500 & whale$number.whales > mean(whale$number.whales), ]

### Question 12: subset functions
#part a: All rows in May
whale_may <- subset(whale, month == "May" & time.at.station < 1000 & depth >1000)

#part b
whale_october <- subset(whale, month == "October" & latitude > 61, select = c(month,latitude, longitude, number.whales)) 

### Question 13: Order dataframes
whale.depth.sort <- whale[order(whale$depth), ]

### Question 14: Order by ascending & descending ###
#part a= ascending
whale.depth.noise <- whale[order(whale$water.noise, whale$depth), ]

#part b= descending
whale.depth.noise <- whale[order(whale$water.noise, -whale$depth), ]

### Question 15: Summarize variables
mean(whale$time.at.station)     # mean time at station
median(whale$depth)             # median depth
length(whale$number.whales)     # number of observations

# mean of whale at each water noise level
mean_whales_noiselevel <- tapply(whale$number.whales, whale$water.noise, mean)

# median of whale at each water noise level
median_whales_noiselevel <- tapply(whale$number.whales, whale$water.noise, median)

# median of whale for each month
median_whales_month <- tapply(whale$number.whales, whale$month, median)

### Question 16: Summarize variables
#part a:FUN = mean. It applies the mean function to each group #water.noise. Groups the data by the water.noise
mean_water_noise <- aggregate(cbind(time.at.station, number.whales, depth, gradient) ~ water.noise, data = whale, FUN = mean, na.rm = TRUE)

#part b:
mean_water_noise_month <- aggregate(cbind(time.at.station, number.whales, depth, gradient) ~ water.noise + month, data = whale, FUN = mean, na.rm = TRUE)

### Question 17: Table and Tabs functions
#part a =  Determine number of observations for each combo
num_obs_water_noise <- table(whale$water.noise)

#part b =  Display number of observations
num_obs_water_noise_month <- table(whale$water.noise, whale$month)

#part c =  Replace with xtabs
num_obs_wn_xtab <- xtabs(~ water.noise, data = whale)
num_obs_wnmonth_xtab <- xtabs(~ water.noise + month, data = whale)

### Question 18: Export the dataframe
exported_file_path <- "data/whale_num.txt"
write.table(whale.num, file = exported_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)

