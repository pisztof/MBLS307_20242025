#Exercise 4
#https://alexd106.github.io/intro2R/exercise_4.html
# Wed Feb 12 20:28:51 2025 ------------------------------

#4
squid <- read.table("data/squid1.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE)
str(squid)
summary(squid)
squid$year_factor <- factor(squid$year)
squid$month_factor <- factor(squid$month)
squid$maturity.stage_factor <- factor(squid$maturity.stage)
str(squid)

#5
with(squid,table(year_factor,month_factor))
with(squid,xtabs(~year_factor + month_factor))
ftable(xtabs(~ maturity.stage_factor + year_factor + month_factor, data = squid))

#6
pdf("output/squid_dotcharts.pdf")
par(mfrow = c(2,2))
dotchart(squid$DML, main = "DML")
dotchart(squid$weight, main = "weight")
dotchart(squid$nid.length, main = "nid length")
dotchart(squid$ovary.weight, main = "ovary weight")
dev.off()

#7
which(squid$nid.length > 400)
squid$nid.length[11]
squid$nid.length[11] <- 43.2

#8
pdf("output/squid_histograms.pdf")
par(mfrow = c(2,2))
hist(squid$DML, breaks = seq(from = 50, to = 350, by = 15), main = "DML")
hist(squid$weight, main = "weight")
hist(squid$nid.length, main = "nid length")
hist(squid$ovary.weight, main = "ovary weight")
dev.off()

#9
png("output/scatter_DML_weight.png")
par(mfrow = c(2,2))
plot(squid$DML,squid$weight)
squid$log_weight <- log(squid$weight)
squid$sqrt_weight <- sqrt(squid$weight)
plot(squid$log_weight,squid$weight)
plot(squid$sqrt_weight,squid$weight)
dev.off()

#10
png("output/box_DML_maturity.png")
par(mfrow = c(1,1))
boxplot(squid$DML~squid$maturity.stage_factor, xlab = "Maturity Stage", ylab = "DML")
dev.off()

#11
png("output/coplot_squid.png")
coplot(sqrt_weight ~ DML | maturity.stage_factor, data = squid)
dev.off()

#12
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
png("output/squid_pairplot.png")
pairs(squid[,c('DML', 'weight', 'eviscerate.weight', 
               'ovary.weight', 'nid.length', 'nid.weight')],
      upper.panel = panel.cor, lower.panel = panel.smooth, diag.panel = panel.hist)
dev.off()

#13
png("output/squid_DML_ovary.png")
plot(squid$DML,sqrt(squid$ovary.weight), type = "n", xlab = "DML", ylab = "sqrt(ovary weight)")
for (i in 1:5){
points(squid$DML[squid$maturity.stage_factor == i],
       sqrt(squid$ovary.weight[squid$maturity.stage_factor == i]), col = i)
}

leg_lab <- c("1", "2", "3","4","5")
legend(x = 90, y = 7, col = 1:5,pch = rep(1,5), 
       legend = leg_lab, bty = "n", 
       title = "Maturity stage")
dev.off()
