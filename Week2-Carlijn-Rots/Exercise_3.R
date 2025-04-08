setwd('/Users/carlijnrots/Documents/Universiteit/Jaar 4/MBLS-307    Omics for the Life Sciences/R')

## 5 
whales <- read.table(file='whaledata.txt', header = TRUE, sep = "\t", stringsAsFactors = TRUE)

## 6 
whales
head(whales)
str(whales)
# Later reference
# 'data.frame':	100 obs. of  8 variables:
#   $ month          : Factor w/ 2 levels "May","October": 1 1 1 1 1 1 1 1 1 1 ...
# $ time.at.station: int  1344 1633 743 1050 1764 580 459 561 709 690 ...
# $ water.noise    : Factor w/ 3 levels "high","low","medium": 2 3 3 3 3 1 3 3 3 3 ...
# $ number.whales  : int  7 13 12 10 12 10 5 8 11 12 ...
# $ latitude       : num  60.4 60.4 60.5 60.3 60.4 ...
# $ longitude      : num  -4.18 -4.19 -4.62 -4.35 -5.2 -5.22 -5.08 -5 -4.64 -4.84 ...
# $ depth          : int  520 559 1006 540 1000 1000 993 988 954 984 ...
# $ gradient       : int  415 405 88 409 97 173 162 162 245 161 ...
## this dataset has 100 observations, 8 variables.

## 7 
summary(whales)

## 8 
whalessub <- whales[1:10, 1:4]
whalesnum <- whales[, c("month","water.noise","number.whales")]
whales.may2 <- whales[1:50,]
whales.may <- whales[whales$month == "May",]
whales[11:100, 1:7]
whales[11:100, c("month", "time.at.station", "water.noise", "number.whales", "latitude", "longitude", "depth")]

## 9 
whales.deep <- whales[whales$depth > 1200, ]
whales.steep <- whales[whales$gradient > 200,]
whales.low <- whales[whales$water.noise == "low",]
whales.high.may <- whales[whales$water.noise == "high" & whales$month == "May", ]
whales.low.oct.highmed <- whales[whales$water.noise == "low" & whales$month == "October" & whales$gradient > median(whales$gradient),]
whales.specific <- whales[whales$latitude >= 60.0 & whales$latitude <= 61.0 & whales$longitude >= -6.0 & whales$longitude <= -4.0 , ]
whales.specific
whales.no_medium <- whales[whales$water.noise != "medium", ]

## 10
# - see above

## 11
whales.problem <- whales[whales$depth > 1500 & whales$number.whales > mean(na.rm = TRUE, whales$number), ]
whales.problem #this problem is now solved. 

## 12 
subset_1 <- subset(whales, month == "May" & time.at.station < 1000 & depth > 1000)
subset_2 <- subset(whales, month == "October" & latitude > 61, select = c("month", "latitude", "longitude", "number.whales"))
subset_2

## 13
whale.depth.sort <- order(whales$depth)

## 14
whales$water.noise <- factor(whales$water.noise, levels = c("low", "medium", "high"))
subset_14 <- whales[order(whales$water.noise,whales$depth),]
subset_14
subset_14.2 <- whales[order(whales$water.noise, -whales$depth),]
subset_14.2

## 15 
mean.whales <- aggregate(number.whales ~ water.noise, FUN = mean, data = whales)
mean.whales
mean.whales.low <- mean(whales$number.whales[whales$water.noise == "low"], na.rm=TRUE)
mean.whales.low
mean.whales.2 <- aggregate(number.whales ~ water.noise + month, FUN = mean, data = whales)
mean.whales.2
mean.whales.low.may <- mean(whales$number.whales[whales$water.noise == "low" & whales$month == "May"], na.rm = T)
mean.whales.low.may

## 16
whales.table <- aggregate(whales[,c("time.at.station", "number.whales", "depth", "gradient")], by = list(noise = whales$water.noise), FUN = mean, na.rm = T)
whales.table
whale.table.month <- aggregate(whales[,c("time.at.station", "number.whales", "depth", "gradient")], by = list(noise = whales$water.noise, month = whales$month), FUN = mean, na.rm = T)
whale.table.month[,-c(1,2)] <- round(whale.table.month[,-c(1,2)], 2)
whale.table.month

## 17 
# nr of observations for water noise
table(whales$water.noise)
# nr of observations for each combi of noise and month
table(whales$water.noise, whales$month)
xtabs(~water.noise + month, data = whales)

## 18
whalesnum
write.table(whalesnum, file = 'whale.num.txt', col.names = T, row.names = F, sep = "\t")

