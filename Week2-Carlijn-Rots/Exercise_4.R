# from intro2R chapter 4 

## 4
squid <- read.table("squid1.txt", header = TRUE)
str(squid)
summary(squid)
squid$year <- factor(squid$year)
squid$month <- factor(squid$month)
squid$maturity.stage <- factor(squid$maturity.stage)
str(squid)

## 5 
xtabs(~ month + year, data = squid)
table(squid$month, squid$year)
ftable(xtabs(~ month + year + maturity.stage, data = squid))

## 6
dotchart(squid$DML)
dotchart(squid$weight)
dotchart(squid$nid.length)
dotchart(squid$ovary.weight)

## 7 
which(squid$nid.length > 400)
squid$nid.length[11]
squid$nid.length[11] <- 43.2
dotchart(squid$nid.length)

## 8 
hist(squid$DML)
hist(squid$weight)
hist(squid$eviscerate.weight)
hist(squid$ovary.weight)
brk <- seq(from = 75, to = 350, by = 1)
hist(squid$DML, breaks = brk)

## 9 
plot(x = squid$DML, y = squid$weight)
squid$weight_sqrt <- sqrt(squid$weight)
plot(x = squid$DML, y = squid$weight_sqrt)
squid$weight_log <- log(squid$weight)
plot(x = squid$DML, y = squid$weight_log)

## 10
boxplot(DML ~ maturity.stage, data = squid, xlab = "state of maturity", ylab = "Length of outer body (DML)", main = "Size of squid in maturity stages")
install.packages("vioplot")
vioplot(DML ~ maturity.stage, data = squid,  xlab = "state of maturity", ylab = "Length of outer body (DML)", main = "Size of squid in maturity stages")

