#Sergey Yudintsev, 1527320
print('Hello world')

#This R script covers exercises 2,3,4 and 6.
#The last task in exercise 4 was omitted due to the remark that the whole
#exercise 4 needs to be done only until the task 10. 

#exercise2
x <- 12.43
vactor1 <- c(log(x), log10(x), sqrt(x), exp(x))
a <- array(vactor1)
print(a)
area_circle <- 20*20*pi
print(area_circle)
d <- (14*0.51)^(1/3)
print(d)
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)
w_a <- c('mean', mean(weight), 'var', var(weight), 'sd', 
         +sd(weight), 'len', length(weight))
print(w_a)
first_five <- weight[c(1, 2, 3, 4, 5)]
print(first_five)
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
dd <- summary(height)
print(dd)
some_child <- height[c(2, 3, 9, 10)]
print(some_child)
shorter_child <- height[height > 98]
print(shorter_child)
BMI <- vector( "numeric" , 10 )
for (k in 1:length(height)) {
  BMI[k] = height[k]/(weight[k]*weight[k])
}
print(BMI)
seq1 <- seq(from = 0, to = 1, by = 0.1)
print(seq1)
seq2 <- rev(seq(from = 1, to = 10, by = 0.5))
print(seq2)
seq11 <- rep(1:3, times = 3)
seqvec <- c("a", "c", "e", "g")
seq22 <- rep(seqvec, each = 3)
seq33 <- rep(seqvec, times = 3)
seq44 <- rep(1:3, each=3, times=3)
seqvec2 <- 1:5
seq55 = vector( "numeric" , 15 )
j <- 1
for (i in 1:length(seqvec2)) {
  print(i)
  print(seqvec2[i])
  for (k in 1:(6-seqvec2[i])) {
    seq55[j] <- seqvec2[i]
    j <- j + 1
  }
}
seqvec3 <- c(7,2,8,1)
seq66 <- rep(seqvec3, times=c(4,3,1,4))
print(seq11)
print(seq22)
print(seq33)
print(seq44)
print(seq55)
print(seq66)
height_sorted <- sort(height)
print(height_sorted)
height_rev_sorted <- rev(sort(height))
print(height_rev_sorted)
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", 
                "Louise", "Mary", "Ronald", "William")
height_ord <- order(height)
weight_ord <- order(weight)
names_sort <- child_name[height_ord]
print(names_sort)
weight_rev <- rev(child_name[weight_ord])
print(weight_rev)



#exercise3
getwd()
#setwd("/Users/sergey/MBLS307_20242025/week2-clayunder")
whale <- read.table("data/data/whaledata.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE)
head(whale)
names(whale)
str(whale)

#'data.frame':	100 obs. of  8 variables:
#$ month          : Factor w/ 2 levels "May","October": 1 1 1 1 1 1 1 1 1 1 ...
#$ time.at.station: int  1344 1633 743 1050 1764 580 459 561 709 690 ...
#$ water.noise    : Factor w/ 3 levels "high","low","medium": 2 3 3 3 3 1 3 3 3 3 ...
#$ number.whales  : int  7 13 12 10 12 10 5 8 11 12 ...
#$ latitude       : num  60.4 60.4 60.5 60.3 60.4 ...
#$ longitude      : num  -4.18 -4.19 -4.62 -4.35 -5.2 -5.22 -5.08 -5 -4.64 -4.84 ...
#$ depth          : int  520 559 1006 540 1000 1000 993 988 954 984 ...
#$ gradient       : int  415 405 88 409 97 173 162 162 245 161 ...

summary(whale)
whale.sub <- whale[1:10, 1:4]
whale.num <- whale[,c("month","water.noise","number.whales")]
whale.may <- whale[1:50,]
whale.exceptfrom <- whale[-(1:10),-8]
#print(whale.sub)
#print(whale.num)
#print(whale.may)
#print(whale.exceptfrom)
whale.deeperthan1200 <- whale[whale$depth > 1200,]
whale.steeperthan200 <- whale[whale$gradient > 200,]
whale.waternoiselow <- whale[whale$water.noise=="low",]
whale.maywaternoisehigh <- whale[whale$water.noise=="high"&whale$month=='May',]
whale.octwaternoiselowsteeperthan132 <- whale[whale$month=='October'&whale$water.noise=="low" 
                                              & whale$gradient > median(whale$gradient),]
whale.specificlonglat<- whale[whale$latitude>=60 & whale$latitude<=61 
                        & whale$longitude>=-6 & whale$longitude<=-4,]
whale.waternoisenotmedium <- whale[whale$water.noise!="medium",]
#print(whale.deeperthan1200)
#print(whale.steeperthan200)
#print(whale.waternoiselow)
#print(whale.maywaternoisehigh)
#print(whale.octwaternoiselowsteeperthan132)
#print(whale.specificlonglat)
#print(whale.waternoisenotmedium)
whale.deepandmorethanmean <- whale[whale$depth>1500 
                         & whale$number.whales>mean(whale$number.whales, na.rm=TRUE),]
print(whale.deepandmorethanmean)
whale.maytimelessthan1000depthmorethan1000 <- subset(whale, month=="May" & time.at.station<1000 & depth>1000)
whale.octlatitudemorethan61selected <- subset(whale, month=="October" & latitude>61, 
                                       select = c("month","latitude","longitude","number.whales"))
#print(whale.maytimelessthan1000depthmorethan1000)
#print(whale.octlatitudemorethan61selected)
whale.depth.sort <- whale[order(whale$depth),]
#print(whale.depth.sort)
whale.depth.noise.sort <- whale[order(whale$water.noise, whale$depth),]
whale.depth.noise.sort_descend <- whale[order(whale$water.noise, -whale$depth),]
#print(whale.depth.noise.sort)
#print(whale.depth.noise.sort_descend)
tapply(whale$number.whales, whale$water.noise, mean, na.rm=TRUE)
with(whale, tapply(number.whales, list(water.noise,month), median, na.rm=TRUE))
aggregate(whale[,c(2, 4, 7, 8)], by=list(noise=whale$water.noise), FUN=mean, na.rm=TRUE)
aggregate(whale[,c(2, 4, 7, 8)], by=list(noise=whale$water.noise, month=whale$month), FUN=mean, na.rm=TRUE)
table(whale$water.noise)
table(whale$water.noise, whale$month)
write.table(whale.num, file='output/whaledata_output.txt', col.names = TRUE, row.names = FALSE, sep = "\t")
whale_output_check <- read.table("output/whaledata_output.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE)
head(whale_output_check)
summary(whale_output_check)



#exercise4
squid <- read.table("data/data/squid1.txt", header=TRUE, sep="\t", stringsAsFactors=TRUE)
str(squid)
summary(squid)
squid$year_f <- factor(squid$year)
squid$month_f <- factor(squid$month)
squid$maturity.stage_f <- factor(squid$maturity.stage)
str(squid)
with(squid,xtabs(~year_f+month_f))
ftable(xtabs(~maturity.stage_f+year_f+month_f, data=squid))

#pdf("output/dotchatrtssquid.pdf")
#par(mfrow=c(2, 2))
#dotchart(squid$DML, main="DML")
#dotchart(squid$weight, main="Weight")
#dotchart(squid$nid.length, main="Nid length")
#dotchart(squid$ovary.weight, main="Ovary weight")
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format as dotchatrtssquid.pdf') 
which(squid$nid.length>400)
print(squid$nid.length[11])
squid$nid.length[11] <- 43.2

#pdf("output/histogramssquid.pdf")
#par(mfrow=c(2,2))
#hist(squid$DML, breaks=seq(from=50, to=350, by=20), main="DML")
#hist(squid$weight, main="Weight")
#hist(squid$eviscerate.weight, main="Eviscerate weight")
#hist(squid$ovary.weight, main="Ovary weight")
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format, as histogramssquid.pdf')

#pdf("output/DML_vs_weight.pdf")
#par(mfrow=c(2,2))
#plot(squid$DML, squid$weight)
#squid$weight.log <- log(squid$weight)
#squid$weight.sqrt <- sqrt(squid$weight)
#plot(squid$DML, squid$weight.log)
#plot(squid$DML, squid$weight.sqrt)
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format, as DML_vs_weight.pdf')
print('Sqrt linearises the relationship the best')

#pdf("output/boxplotDMLmaturitystages.pdf")
#par(mfrow=c(1,1))
#boxplot(squid$DML~squid$maturity.stage_f, ylab="DML", xlab="Maturity Stage")
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format, as boxplotDMLmaturitystages.pdf')
library(vioplot)

#pdf("output/vioplotDMLmaturitystages.pdf")
#par(mfrow=c(1,1))
#vioplot(squid$DML~squid$maturity.stage_f, ylab="DML", xlab="Maturity Stage")
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format, as vioplotDMLmaturitystages.pdf')

#pdf("output/coplotDMLmaturitystages.pdf")
#coplot(weight.sqrt ~ DML | maturity.stage_f, data=squid)
#dev.off()

print('The results for this section can be seen in the output folder in the pdf format, as coplotDMLmaturitystages.pdf')

panel.hist <- function(x, ...)
{
  #from help of pairs
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y)
}
pairs(squid[, c("DML", "weight", "eviscerate.weight", "ovary.weight", "nid.length", "nid.weight")], diag.panel = panel.hist)




