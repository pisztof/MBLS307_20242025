whale <- read.table(file='whaledata.txt', header=TRUE, sep="\t", stringsAsFactors = TRUE)
head(whale)
names(whale)

#7-count
str(whale)
#there are 100 observations of 8 variables

summary(whale)
sum(is.na(whale))
colSums(is.na(whale))
#there is 1 missing value in the column number.whales

#8
#a)
whale.sub <- whale[1:10, 1:4]
whale.sub

#b)data from columns
whale.num <- whale[, c("month","water.noise","number.whales")]
whale.num

#c)
whale.may <- whale[1:50,]
whale.may

#d)
whale.exc <- whale[10:100, 1:7] 

#9 and 10
#at depths > 1200m
whale.dep <- whale[whale$depth > 1200, ]
whale.dep

#at gradient > 200m
whale.gra <- whale[whale$gradient > 200,]
whale.gra

#water noise = low
whale.no <- whale[whale$water.noise=="low",]
whale.no

#water noise = high in MAY
whale.noma <- whale[whale$water.noise == "high" & whale$month == "May",]
whale.noma

#in OCTOBER, WATER NOISE=low and gradient>median
whale.oclo <- whale[whale$month=="October" & whale$water.noise=="low" & whale$gradient>median(whale$gradient,),]
whale.oclo

#latitude between 60 and 61 AND longitude between -6.0 and 4.0
whale.obs <- whale[whale$latitude<=61 & whale$latitude>=60 & whale$longitude<=-4.0 & whale$longitude>=-6.0,]
whale.obs

#all rows withOUT a water noise = medium
whale.nomed <- whale[whale$water.noise != "medium",]
whale.nomed

#11
whale.one <- whale[whale$depth>1500 & whale$number.whales>mean(whale$number.whales),]
whale.one
help(mean)
#the output has NA everywhere, this is probably because the value that is missing is in the number.whales column, meaning that the mean cannot be performed unless we make a little correction
whale.two <- whale[whale$depth>1500 & whale$number.whales>mean(whale$number.whales,na.rm = TRUE),]
whale.two

#12
?subset
whale.el <- subset(whale, month=="May" & time.at.station<1000 & depth>1000)
whale.el

whale.ti <- subset(whale, month=="October" & latitude>61, select = c("month", "latitude", "longitude", "number.whales"))
whale.ti

#13
whale.order <- order(whale$depth) 
whale.names <- whale[whale.order,]
whale.names

#14
whale.order2 <- order(whale$water.noise, whale$depth)
whale.names2 <- whale[whale.order2,]
whale.names2

whale.order3 <- order(whale$water.noise, whale$depth, decreasing = TRUE)
whale.names3 <- whale[whale.order3,]
whale.names3

#15
whale.low <- mean(subset(whale, water.noise == "low")$number.whales, na.rm = TRUE)
whale.low
whale.med <- mean(subset(whale, water.noise == "medium")$number.whales, na.rm = TRUE)
whale.med
whale.high <- mean(subset(whale, water.noise == "high")$number.whales, na.rm = TRUE)
whale.high

whale.low.median <- median(subset(whale, water.noise == "low")$number.whales, na.rm = TRUE)
whale.med.median <- median(subset(whale, water.noise == "medium")$number.whales, na.rm = TRUE)
whale.high.median <- median(subset(whale, water.noise == "high")$number.whales, na.rm = TRUE)

unique(whale$month)
whale.may.median <- median(subset(whale, month == "May")$number.whales, na.rm = TRUE)
whale.october.median <- median(subset(whale, month == "October")$number.whales, na.rm = TRUE)

#16
?aggregate
whale.means <- aggregate(cbind(time.at.station, number.whales, depth, gradient) ~ water.noise, data = whale, FUN = mean, na.rm = TRUE)
whale.means
whale.means.month <- aggregate(cbind(time.at.station, number.whales, depth, gradient) ~ water.noise + month, data = whale, FUN = mean, na.rm = TRUE)
whale.means.month

#17
whale.table <- table(whale$water.noise)
whale.table.month <- table(whale$water.noise, whale$month)
whale.table.month

whale.tab <- xtabs(~ water.noise, data=whale)
whale.tab.month <- xtabs(~ water.noise + month, data=whale)
whale.table.month

#18
write.table(whale.num, file="whale_num.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
