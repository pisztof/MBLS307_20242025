#5
rawdata <- read.table(file = "whaledata.txt", header=TRUE, sep="\t", stringsAsFactors=TRUE)
whaledf <- rawdata[complete.cases(rawdata),]

#6
whaledf
head(whaledf)
names(whaledf)
str(whaledf)
#there are 8 variables, of which month and water.noise are discrete data strings

#7
summary(whaledf)
#Only number.whales has one missing value

#8
whale.sub = whaledf[1:10,1:4]
summary(whale.sub)

whale.num= whaledf[, c("month", "water.noise", "number.whales")]
summary(whale.num)

whale.may = whaledf[1:50,]
summary(whale.may)

whale.rand = whaledf[-(1:10),-(8)]
summary(whale.rand)

#9
whale.deep <- whaledf[whaledf$depth>1200,]
summary(whale.deep)

whale.steep <- whaledf[whaledf$gradient>200,]
summary(whale.steep)

whale.low <- whaledf[whaledf$water.noise=="low",]
summary(whale.low)

whale.mayhigh <- whaledf[whaledf$water.noise=="high"&whaledf$month=="May",]
summary(whale.mayhigh)

whale.octlowsteep <- whaledf[whaledf$month=="October"&whaledf$water.noise=="low"&whaledf$gradient>median(whaledf[,"gradient"]),]
summary(whale.octlowsteep)

whale.geospot <- whaledf[whaledf$latitude >= 60 & whaledf$latitude <= 61 & whaledf$longitude <= -4 & whaledf$longitude >= -6,]
summary(whale.geospot)

whale.highlow <- whaledf[whaledf$water.noise != "medium",]
summary(whale.highlow)

#10
#done already

#11
whale.bigpopdeep <- whaledf[whaledf$depth > 1500&whaledf$number.whales > mean(whaledf[,"number.whales"]),]
summary(whale.bigpopdeep)
whale.verydeep <- whaledf[whaledf$depth > 1500,]
summary(whale.verydeep)
mean(whaledf[,"number.whales"])

#12
whale.mayshortdeep <- subset(whaledf,month=="May"&time.at.station<1000&depth>1000,)
summary(whale.mayshortdeep)

whale.octlocnum <- subset(whaledf,month=='October'&latitude>61,c('month','latitude','longitude','number.whales'))
summary(whale.octlocnum)

#13
whale.depth.sort <- whaledf[order(whaledf$depth),]
head(whale.depth.sort)

#14
whale.noise.depth.sort <- whaledf[order(whaledf$water.noise,-whaledf$depth),]
head(whale.noise.depth.sort)

#15
medians <- list()
for (i in unique(whaledf$water.noise)){
  water.noise.median <- median(whaledf[whaledf$water.noise == i,"number.whales"])
  medians[[i]] <- water.noise.median
}
medians

#16
rounded.mean <- function(lis){
  return(round(mean(lis),2))
}
aggregate(whaledf[c('time.at.station','number.whales','depth','gradient')],list(whaledf[,'water.noise']),rounded.mean)

#17
table(whaledf$water.noise)
table(whaledf[,c('water.noise','month')])

#18
write.table(whale.num, file = 'whale_num.txt', sep = '\t', row.names = FALSE)
