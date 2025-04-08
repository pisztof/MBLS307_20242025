#5
whale <- read.csv("data/whaledata.csv")

#6
str(whale)
summary(whale)
#100 observations, 8 variables character, and character

#7
#number.whales, 1 missing value

#8
whale_sub <- whale[1:10, 1:4]
whale_num <- whale[,c("month", "water.noise","number.whales")]
whale_may <- whale[whale$month == "May",]
whale_less <- whale[-1:-10, -8]

#9
whale.1200 <- whale[whale$depth > 1200, ]
whale.200 <- whale[whale$gradient > 200, ]
whale.low <- whale[whale$water.noise == "low", ]
whale.h.may <- whale[whale$water.noise == "high" & whale$month == "May", ]
whale.subset <- whale[whale$month == "October" & whale$water.noise == "low" & 
                        whale$gradient > 132, ]
whale.lat.long <- whale[whale$latitude > 60 & whale$latitude < 61 & 
                          whale$longitude > -6 & whale$longitude <-4, ]
whale.nomed <- whale[whale$water.noise != "medium", ]

#10
whale.subset <- whale[whale$month == "October" & whale$water.noise == "low" 
                      & whale$gradient > median(whale$gradient), ]

#11
whale_test <- whale[whale$depth >= 1500 & whale$number.whales > 
                      mean(whale$number.whales, na.rm = TRUE),]

#12
subset(whale, month == "May" & time.at.station < 1000 & depth > 1000)
subset(whale, month == "October" & latitude > 61, 
       select = c("month", "latitude", "longitude", "number.whales"))

#13
whale_sort_depth <- whale[order(whale$depth),]

#14
whale$water.noise <- factor(whale$water.noise, levels=c("low", "medium","high"))
whale_sort_depth <- whale[order(whale$water.noise, -whale$depth),]

#15
mean(whale$depth)
tapply(whale$number.whales, whale$water.noise, median,na.rm=TRUE)
whale$month <- factor(whale$month)
tapply(whale$number.whales, list(whale$month, whale$water.noise), median, na.rm =TRUE)

#16
aggregate(cbind(time.at.station, number.whales, depth) ~ water.noise+month, 
          data = whale, 
          FUN = function(x) {round(mean(x, na.rm = TRUE),2)})

#17
xtabs(~water.noise + month, data=whale,subset=depth <1000)

#18
write.table(whale_num, file = "data/whale_num.txt", col.names=TRUE, row.names=FALSE,
            sep = "\t")
