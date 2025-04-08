#2-find natural log, log to 10 and to 2 and the antilog of 12.43
log(12.43)
log10(12.43)
log2(12.43)
exp(12.43)

#3-area of circle with d=20cm
area_circle <- pi*(10**2)
area_circle

#4-cube root of 14x0.51
cube_root <- (14*0.51)**(1/3)
cube_root

#5-weight vector
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)

#6-mean, variance, standard dev and range
mean_weight <- mean(weight)
variance_weight <- var(weight)
std_weight <- sd(weight)
range_weight <- range(weight)
number_weight <- length(weight)

mean_weight
variance_weight
std_weight
range_weight
number_weight

first_five <- weight[1:5]
first_five

#7-height vector
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
summary(height)

some_child <- height[c(2,3,9,10)]
some_child

shorter_child <- height[height<=99]
shorter_child

#8-BMI
bmi <- weight/((height/100)**2)
bmi

#9-seq1
seq1 <- seq(from = 10, to = 1, by=-0.1)
seq1

#10-seq2
seq2 <- seq(from=10, to=1, by=-0.5)
seq2

#11-rep()
nums <- c(1,2,3)
nums_rep <- rep(nums, times=3)
nums_rep

aceg <- c(rep("a", times=3), rep("c", times=3), rep("e", times=3), rep("g", times=3))
aceg
aceg2 <- c("a", "c", "e", "g")

aceg2_rep <- rep(aceg2, times=3)
aceg2_rep

nums2 <- c(rep(1, times=3), rep(2, times=3), rep(3, times=3))
nums2_rep <- rep(nums2, times=2)
nums2_rep

nums3 <- 1:5
nums3_2 <- 5:1
nums3_rep <- rep(nums3, times=nums3_2)
nums3_rep

nums4_rep <- rep(c(7,2,8,1), times=c(4,3,1,5))
nums4_rep

#12-order
height_sorter <- sort(height)
height_sorter
height_sorter_inverse <- sort(height, decreasing = TRUE)
height_sorter_inverse

#13-names
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")

#14-order
height_ord <- order(height)
height_ord

names_sort <- child_name[height_ord]
names_sort_ascending <- sort(names_sort, decreasing = FALSE)
names_sort_ascending

#the shortest child is Alfred and the tallest is William

#15-descending
weight_ord <- order(weight)
names_sort_weight <- child_name[weight_ord]
names_sort_weight_descending <- sort(names_sort_weight, decreasing = TRUE)
names_sort_weight_descending

#the heaviest is William and the lightest is Alfred

#16-NA
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mean(mydata)
help(mean)

mean(mydata, na.rm = TRUE)

#17-list all variables
ls()
rm(seq1)
ls()
