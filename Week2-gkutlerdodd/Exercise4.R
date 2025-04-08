#Exercise 4
library(tidyverse)
library(readxl)
library(vioplot)

#4.4
squid <- read_excel("Data/squid1.xlsx")
glimpse(squid)
squid <- squid %>% mutate(
  year_factor = factor(year),
  month_factor = factor(month),
  maturity.stage_factor = factor(maturity.stage)
)
glimpse(squid)

#4.5
xtabs(~year_factor + month_factor, data=squid)

#4.6
pdf("Plots/Exercise4_6.pdf")
par(mfrow=c(2,2))
dotchart(squid$DML, xlab="DML")
dotchart(squid$weight, xlab="weight")
dotchart(squid$nid.length, xlab="nid.length")
dotchart(squid$ovary.weight, xlab="ovary.weight")
dev.off()

#4.7
which(squid$nid.length > 400)
squid$nid.length[11] <- 43.2
par(mfrow=c(1,1))
dotchart(squid$nid.length, xlab="nid.length")

#4.8
pdf("Plots/Exercise4_8.pdf")
par(mfrow=c(2,2))
hist(squid$DML, xlab="DML")
hist(squid$weight, xlab="weight")
hist(squid$eviscerate.weight, xlab="eviscerate.weight")
hist(squid$ovary.weight, xlab="ovary.weight")
dev.off()

#4.9
png("Plots/Exercise4_9.png", width=2000, height=666, res=100)
par(mfrow=c(1,3))
squid <- squid %>% mutate(
  weight_log = log(weight),
  weight_sqrt = sqrt(weight)
)
plot(squid$DML, squid$weight, main="No transform")
plot(squid$DML, squid$weight_log, main="Log transform")
plot(squid$DML, squid$weight_sqrt, main="Sqrt transform")
dev.off()

#4.10
par(mfrow=c(1,1))
boxplot(DML ~ maturity.stage_factor, data=squid)
vioplot(DML ~ maturity.stage_factor, data=squid)