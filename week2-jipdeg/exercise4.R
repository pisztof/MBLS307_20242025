#4
squid <- read.table('squid1.txt', header=TRUE, sep = '\t')
str(squid)
summary(squid)
#519 observations with 13 variables

for (i in c('year','maturity.stage')){
  squid[[i]] <- factor(squid[[i]])
}
squid[['month']] <- factor(squid[['month']], levels = 1:12,labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
str(squid)
summary(squid)

#5
table(squid[['month']], squid[['year']])
#There is data from dec 1989 till Juli 1991, also there is way more data in september and november. 
ftable(xtabs(~year+maturity.stage+month, data=squid))

#6
par(mfrow = c(2,2))
dotchart(squid$DML)
dotchart(squid$weight)
dotchart(squid$nid.length)
dotchart(squid$ovary.weight)
#small variables are sometimes hard to see, but it doesn't strike as too bad. There is one big value in nid.length
pdf('dotcharts_squid')

#7
#value was written down incorrectly with a 0 inbetween
which(squid$nid.length>400)
squid$nid.length[11] <- 43.2
squid$nid.length[11]
par(mfrow = c(1,1))
dotchart(squid$nid.length)

#8
par(mfrow = c(2,2))
hist(squid$DML, breaks = 20)
hist(squid$weight)
hist(squid$nid.length)
hist(squid$ovary.weight)

#9
par(mfrow=c(1,1))
scatter.smooth(squid$DML,squid$weight)
#this is a quadratic or logarithmic relationship, but which is hard to see. 
squid$weight.log <- log(squid$weight)
squid$weight.sqrt <- sqrt(squid$weight)
scatter.smooth(squid$DML,squid$weight.log)
scatter.smooth(squid$DML,squid$weight.sqrt)
#sqrt gives the best linear approximation


#10
boxplot(DML~maturity.stage, data=squid)

#11
#STOP, don't need to continue