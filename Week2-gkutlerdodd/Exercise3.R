#Exercise 3
library(tidyverse)
library(readxl)

#3.5
whale <- read_excel("Data/whaledata.xls", na="NA")

#3.6
glimpse(whale)

#3.7
summary(whale)
#number.whales has one missing value

#3.8
whale.sub <- whale[1:10, 1:4]
whale.num <- whale[,c("month", "water.noise", "number.whales")]
whale.may <- whale[1:50,]
whale.subset <- whale[-10:-1,-8]

#3.9
whale.deep <- whale[whale$depth > 1200,]
whale.steep <- whale[whale$gradient > 200,]
whale.quiet <- whale[whale$water.noise == "low",]
whale.loud_may <- subset(whale, water.noise == "high" & month == "May")
whale.many_conditions <- subset(whale, month == "October" & water.noise == "low" & gradient > median(gradient))
whale.location <- subset(whale, latitude >= 60 & latitude <= 61 & longitude >= -6 & longitude <= -4)
whale.not_medium <- whale[whale$water.noise != "medium",]

#3.11
subset(whale, depth > 1500 & number.whales > mean(number.whales)) #Returns empty tibble
subset(whale, depth > 1500 & number.whales > mean(number.whales, na.rm=TRUE))

#3.12
subset(whale, month == "May" & time.at.station < 1000 & depth > 1000)
subset(whale, month == "October" & latitude > 61, c("month", "latitude", "longitude", "number.whales"))

#3.13
whale.depth.sort <- whale[order(whale$depth),]

#3.14
whale <- whale %>% mutate(water.noise = factor(water.noise, c("low", "medium", "high")))
whale.noise_depth_asc.sort <- whale[order(whale$water.noise, whale$depth),]
whale.noise_depth_des.sort <- whale[order(whale$water.noise, -whale$depth),]

#3.15
tapply(whale$number.whales, whale$water.noise, mean, na.rm=TRUE)
tapply(whale$number.whales, whale$water.noise, median, na.rm=TRUE)
tapply(whale$number.whales, whale$month, median, na.rm=TRUE)

#3.16
aggregate(whale[,c("time.at.station", "number.whales", "depth", "gradient")], by=list(noise=whale$water.noise), mean, na.rm=TRUE)
aggregate(whale[,c("time.at.station", "number.whales", "depth", "gradient")], by=list(noise=whale$water.noise, month=whale$month), mean, na.rm=TRUE)

#3.17
table(whale$water.noise)
table(whale$water.noise, whale$month)
xtabs(~water.noise + month, data=whale)

#3.18
write.table(whale.num, file="whale_num.txt", row.names=FALSE)