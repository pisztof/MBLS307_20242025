## 2
print(log10(12.34))
print(log2(12.34))
print(log(12.34))
print(sqrt(12.34))

## 3
area_circle <- (1/2*20)^2 * pi
print(area_circle)

## 4 
cube_root <- (14*0.51)^(1/3)
print(cube_root)

## 5
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)

## 6 
mean(weight)
var(weight)
sd(weight)
length(weight)
range(weight)
first_five <- weight[1:5]
first_five

## 7
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
summary(height)
some_child <- height[c(2,3,9,10)]
shorter_child <- height[height <= 99]
shorter_child

## 8 
BMI <- weight/(height/100)^2
BMI

## 9 
seq1 <- seq(from=0, to=1, by=0.1)

## 10
seq2 <- rev(seq(from=1, to=10, by=0.1))

## 11
seq3 <- rep(1:3, times=3)
seq3.2 <- rep(c(1, 2, 3), times=3)
seq3
seq3.2
seq4 <- rep(c("a", "c", "e", "g"), each=3)
seq4
seq5 <- rep(c("a", "c", "e", "g"), times=3)
seq5
seq6 <- rep(c(rep(c(1, 2, 3), each=3), times=2))
seq6        
seq7 <- rep(1:5, times=c(5:1))
seq7
seq8 <- rep(c(7,2,8,1), times=c(4, 3,1,4))
seq8

## 12
height_sorted <- sort(height)
height_LtoS <- sort(height_sorted, decreasing = TRUE)

## 13 
child_name <- c('Alfred', 'Barbara', 'James', 'Jane', 'John', 'Judy', 'Louise', 'Mary', 'Ronald', 'William')

## 14
height_order <- order(height)
names_sort <- child_name[height_order]
print(paste("The shortest child is",names_sort[1]))
print(paste("The tallest child is",names_sort[length(height)]))

## 15
weight_order <- order(weight)
names_sort_w <- rev(child_name[weight_order])
print(paste("The heaviest child is",names_sort_w[1]))
print(paste("The lightest child is",names_sort_w[length(weight)]))

## 16
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mean(mydata)
mean(mydata, na.rm=TRUE)

## 17 
#?


