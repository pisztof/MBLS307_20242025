##Exercise 4: Visualising data using base R and lattice graphics

## Importing the ‘squid1.txt’ file into R
squid <- read.table(file = 'data/squid1.txt', header = TRUE, sep = "\t", stringsAsFactors = TRUE, dec = ",")

###Question 4:
summary(squid)
str(squid)
squid$year_factor <- factor(squid$year)
squid$month_factor <- factor(squid$month)
squid$maturity_stage_factor <- factor(squid$maturity.stage)
str (squid)


### Question 5:
## part a: How many observations are there per month and year combination?
xtabs(~ year_factor + month_factor, data = squid)

##part b: Do you have data for each month in each year? 
# No, because there is no data for Jan - Nov in 1989, no data for Feb & June 1990 & no data between August - Dec 1991

##part c: Which years have the most observations?
yearly_counts <- rowSums(table(squid$year_factor, squid$month_factor))

###Question 6:
##saving dotplots in pdf version
pdf("data/squid_dotplots.pdf")
## Creating dotcharts 
dotchart(squid$DML, main = "Dotplot of DML", xlab = "DML")
dotchart(squid$weight, main = "Dotplot of Weight", xlab = "Weight")
dotchart(squid$nid.length, main = "Dotplot of Nid.length", xlab = "Nid.length")
dotchart(squid$ovary.weight, main = "Dotplot of Ovary.weight", xlab = "Ovary.Weight")
dev.off()

###Question 7:
## Function to identify which observation is greater than 400
which(squid$nid.length > 400)
## 11th observation is greater than 400

##to confirm if 11th observation is actually greater than 400
squid$nid.length [11]

##to change the value to 43.2
squid$nid.length [11] <-43.2

##re-create dot plot with new values of observation 11
dotchart(squid$nid.length, main = "Corrected dotchart of nid length", xlab = "Nid.length")


### Question 8: 
## To create histograms with different values and export to pdf with first creating a pdf and then putting information inside
pdf("data/squid_histogram.pdf")
hist(squid$DML, main = "Histogram of DML", xlab = "DML")
hist(squid$weight, main = "Histogram of Weight", xlab = "Weight")
hist(squid$eviscerate.weight, main = "Histogram of eviscerate.weight", xlab = "eviscerate.weight")
hist(squid$ovary.weight, main = "Histogram of Ovary.weight", xlab = "Ovary.Weight")
dev.off()

###Question 9 
##Plot the relationship between DML on the x axis and weight on the y axis.
plot(x = squid$DML, y = squid$weight)
# the relationship appears somewhat linear but also a bit exponential

## transforming the weight variable with either a natural log (log()) or square root (sqrt()) transformation
squid$log_weight <- log(squid$weight)
squid$sqrt_weight <- sqrt (squid$weight)
# Doing sqrt weight linearizes there relationship.

#Scatterplot of DML vs sqrt(weight)
plot(x = squid$DML, y = squid$sqrt_weight)

##Scatterplot of DML vs log(weight)
plot(x = squid$DML, y = squid$log_weight)

##saving plot as pdf 
pdf("data/squid_scatterplot log weight.pdf")
plot(x = squid$DML, y = squid$log_weight)
dev.off()

pdf("data/squid_scatterplot sqrt weight.pdf")
plot(x = squid$DML, y = squid$sqrt_weight)
dev.off()

###Question 10: Making box plot
boxplot(squid$DML ~ squid$maturity.stage, data = squid, 
        ylab = "DML", xlab = "Maturity stage")

### installed the package vioplot using instal package (vioplot)
vioplot(squid$DML ~ squid$maturity.stage, data = squid, 
        ylab = "DML", xlab = "Maturity stage",
        col = "lightblue")

