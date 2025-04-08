#Exercises from: https://alexd106.github.io/intro2R/exercise_2.html
install.packages("styler")

#2
a <- c(log(12.43), log10(12.43), log2(12.43), sqrt(12.43), exp(12.43))
print(a)

#3
diameter <- 20
area_circle <- pi*(0.5*diameter)^2

#4
(14 * 0.51)^(1/3)

#5
weights <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)

#6
mean(weights)
sd(weights)
range(weights)
length(weights)
var(weights)
first_five <- weights[1:5]

#7
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
summary(height)
some_child <- height[c(2,3,9,10)]
shorter_child <- height[height<=99]

#8
mbi <- weights/((height/100)^2)

#9
seq1 <- seq(0,1,0.1)

#10
seq2 <- seq(10,1,-0.5)

#11
rep(c(1,2,3),times = 3)
rep(c("a","c","e","g"), each = 3)
rep(c("a","c","e","g"), times = 3)
rep(rep(c(1,2,3),each = 3),times=2)
rep(1:5, times=5:1)
rep(c(7, 2, 8, 1), times = c(4, 3, 1, 5))
    
#12
height_sorted <- order(height)
height_sorted_decreasing <- order(height, decreasing = TRUE)

#13
child_name <- c("Alfred", "Barbara", "James", 
                "Jane", "John", "Judy", "Louise", "Mary", "Ronald", 
                "William")

#14
names_sort <- child_name[height_sorted]

#15
weight_rev <- child_name[order(weights, decreasing = TRUE)]

#16
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mean(mydata, na.rm = TRUE)

#17
ls()
rm(seq1)
ls() 
