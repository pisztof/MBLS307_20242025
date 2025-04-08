#4
squid <- read.table(file='squid1.txt', header=TRUE, sep="\t", stringsAsFactors = TRUE)
str(squid)
#there are 519 observations if 13 variables

squid$year.factor <- factor(squid$year)
squid$month.factor <- factor(squid$month)
squid$maturity.stage.factor <- factor(squid$maturity.stage)
str(squid)

#5
squid.month.year <- xtabs(~ year.factor + month.factor, data=squid)
squid.month.year
#there are some months without data
#the years with most observations are 1990 amd 1991

squid.month.year.maturity <- xtabs(~ year.factor + month.factor + maturity.stage.factor, data = squid)
squid.month.year.maturity.ftable <- ftable(squid.month.year.maturity)

#6
pair.plot <-pairs(squid[, c("DML", "weight", "nid.length", "ovary.weight")])

pdf("dml_dot.pdf")
dotchart(squid$DML)
dev.off()

pdf("weight_dot.pdf")
dotchart(squid$weight)
dev.off()

pdf("nid_length_dot.pdf")
dotchart(squid$nid.length)
dev.off()

pdf("ovary_weight_dot.pdf")
dotchart(squid$ovary.weight)
dev.off()

#7
which(squid$nid.length > 400)
squid$nid.length[11] <- 43.2 
squid$nid.length[11]

#8
brk <- seq(from=0, to=900, by=50)

pdf("weight_hist.pdf")
hist(squid$weight, breaks = brk)
dev.off()

brk.dml <- seq(from=70, to=350, by=25)

pdf("dml_hist.pdf")
hist(squid$DML, breaks = brk.dml)
dev.off()

pdf("eviscerate_weight_hist.pdf")
hist(squid$eviscerate.weight)
dev.off()

pdf("ovary_weight_hist.pdf")
hist(squid$ovary.weight)
dev.off()

#9
pdf("dml_weight_scatter.pdf")
plot(x=squid$DML, y=squid$weight)
dev.off()
#it seems to be exponential

weight_log <- log(squid$weight)
pdf("dml_log_weight_scatter.pdf")
plot(x=squid$DML, y=weight_log)
dev.off()


weight_sqrt <- sqrt(squid$weight)
pdf("dml_sqrt_weight_scatter.pdf")
plot(x=squid$DML, y=weight_sqrt)
dev.off()
#with the square root of weight we manage to linearise the scatter plot

#10 boxplot to visualise dml differences
#at each maturity stage (factors-Q4)
#x and y labels
boxplot(DML ~ maturity.stage.factor, data = squid, ylab = "DML", xlab = "Maturity stage")

install.packages("vioplot")
install.packages("sm")

library(vioplot)
vioplot(DML ~ maturity.stage.factor, data = squid, ylab = "DML", xlab = "Maturity stage", col="pink")
