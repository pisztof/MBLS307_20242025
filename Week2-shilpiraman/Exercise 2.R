#natural log
log (12.43)

#Log to the base 10
log10(12.43)

#Log to the base 2
log2(12.43)

#Square root of the value 
sqrt(12.43)

#natural antilog of 12.43.
exp(12.43)

#To calculate area of a circle
# Area = pi*r2 (square)
d <- 20
r <- d/2
area_circle <- pi*r^2
print (area_circle)

#to calculate the cube root of 14 x 0.51
value <- 14*0.51 
print (value)
return(value^(1/3))
cube_root <- function (value)
print (cube_root)



#to create a vector called weight containing the weight (in kg) of 10 children: 
c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66)

#In the code below we’ve created an object called my_vec and assigned it a value using the function c(). 
my_vec <-c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66) 

#to calculate mean of the vector values
mean(my_vec)

#to calculate variance of the vector values
var(my_vec)

#to calculate standard deviation of the vector values
sd (my_vec)

# returns the number of elements in my_vec
length(my_vec)  


my_seq2 <- seq(from = 0, to = 10, by = 2)
print (my_seq2)


# extract the weights for the first five children using Positional indexes
my_vec 
my_vec [1:5]

# if you want to store this value in another object
first_five <- my_vec[1:5]
first_five

#to create a vector called heigh containing the height (in cms) of 10 children: 
c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112)
my_vec2 <-c (112, 102, 83, 84, 99, 90, 77, 112, 133, 112)

summary()

#extract the height of the 2nd, 3rd, 9th and 10th child and assign these heights to a variable called some_child
my_vec [c(2,3,9,10)]
some_child <-my_vec [c(2,3,9,10)]

#extract all the heights of children less than or equal to 99 cm and assign to a variable called shorter_child
my_vec2 [my_vec2 < 99]
my_vec2 [my_vec2 <= 99]
shorter_child <- my_vec2 [my_vec2 <= 99]

#calculate the body mass index (BMI) for each child
# formula for BMI = weight (in kg) divided by the square of the height (in meters)
weight <- my_vec
height_cms <- my_vec2
height_ms <- my_vec2/100
bmi <- my_vec/height_ms^2
print (bmi)


#create a sequence of numbers ranging from 0 to 1 in steps of 0.1
my_seq <- seq(from=1, to=10, by=0.1)
print (my_seq)
seq1 <- my_seq <- seq(from=1, to=10, by=0.1)

#create a sequence of numbers ranging from 0 to 1 in steps of 0.1
my_seq <- seq (from=1, to=10, by=0.5)
seq2 <-my_seq <- seq (from=1, to=10, by=0.5)
print (seq2)

#create repeating seq
my_seq4 <- rep("123", times = 3)    # repeats ‘123’ 3 times 
my_seq5 <- rep(1:3, times = 3)  # repeats the series 1 to 3, 3 times

#to sort the height and assigned to height_sorted
vec_sort <- sort(my_vec2)
height_sorted <-vec_sort <- sort(my_vec2)

#naming children 
c ("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")

#Create a new variable called names_sort to store the names of the children ordered by child height (from shortest to tallest)
height_ord <- order (my_vec2)
height_ord
names_ord <- child_name[height_ord]

#to sort in descending order
vec_sort2 <- sort(my_vec, decreasing = TRUE)
weight_rev <-vec_sort2 <- sort(my_vec, decreasing = TRUE)

#treating missing values, na.rmm ignores the na from the list 
c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7)
mean (mydata)
mean_data <- mean(mydata, na.rm = TRUE)
mean_data

ls()
rm(seq1)
ls()
