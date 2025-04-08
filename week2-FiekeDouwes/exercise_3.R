#Exercise 3
#https://alexd106.github.io/intro2R/exercise_3.html
# Wed Feb 12 15:04:28 2025 ------------------------------

# > str(whale)
# 'data.frame':	100 obs. of  8 variables:
#   $ month          : Factor w/ 2 levels "May","October": 1 1 1 1 1 1 1 1 1 1 ...
# $ time.at.station: int  1344 1633 743 1050 1764 580 459 561 709 690 ...
# $ water.noise    : Factor w/ 3 levels "high","low","medium": 2 3 3 3 3 1 3 3 3 3 ...
# $ number.whales  : int  7 13 12 10 12 10 5 8 11 12 ...
# $ latitude       : num  60.4 60.4 60.5 60.3 60.4 ...
# $ longitude      : num  -4.18 -4.19 -4.62 -4.35 -5.2 -5.22 -5.08 -5 -4.64 -4.84 ...
# $ depth          : int  520 559 1006 540 1000 1000 993 988 954 984 ...
# $ gradient       : int  415 405 88 409 97 173 162 162 245 161 ...

#5
whale <- read.table('data/whaledata.txt', header = TRUE, sep="\t", stringsAsFactors = TRUE)

#6
head(whale)
names(whale)
str(whale)

#7
summary(whale)

#8
whale.sub <- whale[1:10, 1:4]
whale.num <- whale[ , c("month", "water.noise", "number.whales")]
whale.may <- whale[1:50,]
whale.idk <- whale[-(1:10),-8]

#9
whale.1200m <- whale[whale$depth > 1200,]
whale.200deg <- whale[whale$gradient > 200,]
whale.lownoise <- whale[whale$water.noise == "low",]
whale.highnoise_may <- whale[whale$water.noise == "high" & whale$month == 'May',]
whale.oct_lownoise_steep <- whale[whale$water.noise == "low" 
                                  & whale$month == 'October' 
                                  & whale$gradient > 132,]
whale.lat_long <- whale[whale$latitude >= 60 & whale$latitude <= 61 
                        & whale$longitude >= -6 & whale$longitude <= -4, ]
whale.lat_long <- whale[whale$water.noise != "medium",]

#10
whale.oct_lownoise_steep <- whale[whale$water.noise == "low" 
                                  & whale$month == 'October' 
                                  & whale$gradient > median(whale$gradient),]

#11
whale.deep_many <- whale[whale$depth > 1500 
                         & whale$number.whales > mean(whale$number.whales,na.rm = TRUE),]

#12
whale.may_short_deep <- subset(whale, month == "May" & time.at.station < 1000 & depth > 1000)
whale.oct_lat_61 <- subset(whale, month == "October" & latitude > 61, 
                           select = c("month", "latitude", "longitude", "number.whales"))

#13
whale.depth.sort <- whale[order(whale$depth),]

#14
whale.noise.depth.sort <- whale[order(whale$water.noise,whale$depth),]
whale.noise.depth.sort.descend <- whale[order(whale$water.noise,-whale$depth),]

#15
tapply(whale$number.whales, whale$water.noise, mean, na.rm=TRUE)
with(whale, tapply(number.whales, list(water.noise, month),median, na.rm=TRUE))

#16
aggregate(whale[ , c(2, 4, 7, 8)], 
          by = list(noise = whale$water.noise), 
          FUN = mean, na.rm = TRUE)
aggregate(whale[ , c(2, 4, 7, 8)], 
           by = list(noise = whale$water.noise, month = whale$month), 
           FUN = mean, na.rm = TRUE)

#17
table(whale$water.noise)
table(whale$water.noise,whale$month)

#18
write.table(whale.num, "output/whale_num.txt", sep="\t", row.names = FALSE, col.names = TRUE)
