#4
squid <- read.csv("data/squid1.csv")
str(squid)
summary(squid)
squid$Fmaturity <- factor(squid$maturity.stage)
squid$Fmonth <- factor(squid$month) 
squid$Fyear <- factor(squid$year)

#5
xtabs(~Fmonth+Fyear+Fmaturity, data=squid)
ftable(xtabs(~Fmonth+Fyear+Fmaturity, data=squid))

#6
pdf("data/badFigureExercise4.6.pdf")
par(mfrow = c(2, 2))
boxplot(squid$DML, main="DML")
boxplot(squid$weight, mean="weight")
boxplot(squid$nid.length, main="nid.length")
boxplot(squid$ovary.weight, main="overay.weight")
par(mfrow = c(1, 1))
dev.off()
#weight, nid.length and ovary.weight

#7
which(squid$nid.length > 400)
squid$nid.length[11] <- 43.2

#8
pdf("data/badFigureExercise4.8.pdf")
par(mfrow = c(2, 2))
hist(squid$DML, main="DML",breaks=100)
hist(squid$weight, main="weight")
hist(squid$nid.length, main="nid.length")
hist(squid$ovary.weight, main="overay.weight")
par(mfrow = c(1, 1))
dev.off()

#9
squid$weightlog <- log(squid$weight)
squid$weightsqrt <- sqrt(squid$weight)
pdf("data/badFigureExercise4.9.pdf")
par(mfrow = c(3, 1))
plot(squid$DML, squid$weight)
plot(squid$DML, squid$weightlog)
plot(squid$DML, squid$weightsqrt)
par(mfrow = c(1, 1))
dev.off()

#10
boxplot(DML ~ Fmaturity, data=squid, ylab="DML", xlab="maturity")
install.packages("vioplot")
library(vioplot)
vioplot(DML ~ Fmaturity, data = squid, xlab = "maturity stage", ylab = "DML",
        col = "lightblue")
