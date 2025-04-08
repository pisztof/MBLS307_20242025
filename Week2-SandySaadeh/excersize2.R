##Find the natural log, the natural log is log() without specifying the base
log(12.43)
##log to the base 10, 
log10(12.43)
##log to the base 2, 
log2(12.43)
##square root
sqrt(12.43)
##the natural antilog of 12.43
exp(12.43)

area_circle <- pi*sqrt(20)
area_circle
nthroot(14*0.51,3)
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)
weight
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
height
height1 <- height/100
height1
#weight (in kg) divided by the square of the height (in meters)
bmi <- weight/(height1)^2
bmi
seq1 <- seq(from = 0, to = 1, by = 0.1)
seq1
seq2 <- rep(1:3, times = 3)
seq2
seq3 <- rep(c("a","b","c"), each=3)
seq3
seq4 <- rep(c("a", "c","e","g" ), times=3)
seq4
seq5 <- rep(1:3, times=2, each=3)
seq5
#111112222333445
seq6<- rep(1:5, times=c(5,4,3,2,1))
seq6
#7 7 7 7 2 2 2 8 1 1 1 1 1
seq7 <- rep(c(7,2,8,1), times=c(4,3,1,5))
seq7
height_sorted <- sort(height, decreasing = TRUE)
height_sorted
height_descending <- sort(height)
height_descending
vec1 <- c("Alfred","Barbara","James","Jane","John","Judy","Louise","Mary","Ronald","William")
vec1
order_height <- order(height_sorted)
order_height
ordering_height_descending <- order(height_descending)
ordering_height_descending
names_sort <- vec1[order_height]
names_sort
weight_rev <- vec1[rev(order_height)]
weight_rev
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mydata
mean(mydata, na.rm= TRUE)
ls()
rm(seq1)
ls()
rm(seq1.rm)
ls()
