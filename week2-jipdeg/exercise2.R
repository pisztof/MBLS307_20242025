#2
n <- 12.43
natural <- log(n)
print(natural)

base10 <- log10(n)
print(base10)

base2 <- log2(n)
print(base2)

sqroot <- sqrt(n)
print(sqroot)

natantilog <- exp(n)
print(natantilog)

#3
d <- 20
r <- d/2
area <- pi*r^2
print(area)

#4
base <- 14*0.51
cuberoot <- base^(1/4)
print(cuberoot)

#5
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)
print(weight)

#6
meankg <- mean(weight)
print(meankg)

varkg <- var(weight)
print(varkg)

stdkg <- stats::sd(weight)
print(stdkg)

rangekg <- range(weight)
print(rangekg)

lnkg <- length(weight)
print(lnkg)

first_five <- weight[0:5]
print(first_five)

#7
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
summary(height)

some_child <- height[c(2,3,9,10)]
print(some_child)

shorter_child <- height[height <= 99]
print(shorter_child)

#8
bmi <- weight / height^2
print(bmi)

#9
seq1 <- seq(0,1,0.1)
print(seq1)

#10
seq2 <- seq(10,1,-0.5)
print(seq2)

#11
print(rep(seq(1,3),3))
print(rep(c("a", "c", "e", "g"), each = 3))
print(rep(c('a','c','e','g'),3))
print(rep(rep(seq(1,3),each=3),2))
print((rep(seq(1,5),seq(5,1,-1))))
print(rep(c(7,2,8,1),c(4,3,1,5)))

#12
height_sorted <- sort(height)
print(height_sorted)
height_revsort <- rev(height_sorted)
print(height_revsort)

#13
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")
print(child_name)

#14
names_sort <- child_name[order(height)]
print(names_sort)
#Louise is shortest, Ronald the longest

#15
weight_rev <- rev(child_name[order(weight)])
print(weight_rev)

#16
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
print(mean(mydata))
print(mean(mydata,0,TRUE))

#17
ls()
rm(seq1)
ls()

