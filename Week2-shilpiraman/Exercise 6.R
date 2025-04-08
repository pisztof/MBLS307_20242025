##Exercise 6: Basic programming in R


## Question 1: 1. Create a function to calculate the area of a circle. Test the function by finding the area of a circle with a diameter of 3.4 cm. 

diameter <- 3.4
radius <-diameter/2
area <-pi * radius^2
print (area)

circle.area <- function(d){ 
  pi * (d/2)^2
}
# to use new function
circle.area(3.4)
# 9.079203

# to test on a vector of diameter
# create a vector with diameters ranging from 0 to 10 in steps of 2

cir.diam <- seq(from = 0, to = 10, by = 2)

# testing function
circle.area(cir.diam)


## Question 2: Function to convert farenheit into centegrade (oC = (oF - 32) x 5/9)

far.cent <- function(a) {
  val <- (a - 32) * 5/9
  print(paste("Fahrenheit: ", round(a, digits = 3), "oF", sep = " "), quote = FALSE)  # round 3dp
  print(paste("Centigrade: ", round(val, digits = 3), "oC", sep = " "), quote = FALSE)  # round 3dp
}


## Question 3: create a vector of normally distributed data
data_vector <-rnorm(1:100, mean=35, sd = 15)
#function to calculate mean, median and range of vector
data_mean <-mean(data_vector)
print (data_mean)

#function to calculate median
data_median<-median (data_vector)
print(data_median)

#function to calculate range of vector
data_range <-range(data_vector)
print (data_range)

#plot histogram with density curve
hist(data_vector, probability = TRUE, col = "lightblue", border = "black", main = "Histogram with density curve", xlab = "values")
lines(density(data_vector), col = "red", lwd=2)

## Question 4: calculate the median value of a vector of numbers 
c(1:10)
# calculating median
ourmedian <- function(x){
  n <- length(x)
  if (n %% 2 == 1)      # odd numbers
    sort(x)[(n + 1)/2]  # find the middle number by adding 1 to length and div 2
  else {                # even numbers
    middletwo <- sort(x)[(n/2) + 0:1]   #find the two middle numbers
    mean(middletwo)
  }
}

# use the function
mydat <- sample(1:50, size = 10, replace = TRUE )

# our function
ourmedian(mydat)

# R median function
median(mydat)

###Question 5:
# function to simulate Ricker model

Ricker.model <- function(nzero, r, time, K = 1) {
  # sets initial parameters
  N <- numeric(time + 1)  # creates a real vector of length time+1 to store values of Nt+1
  N[1] <- nzero  # sets initial population size in first element of N
  for (i in 1:time) {
    # loops over time
    N[i + 1] <- N[i] * exp(r * (1 - N[i]/K))
  }
  Time <- 0:time  # creates vector for x axis
  plot(Time, N, type = "o", xlim = c(0, time), xlab = "Time", ylab = "Population size (N)", main = paste("r =",
                                                                                                         r, sep = " "))  # plots     output
}

# To run play around with the different parameters, especially r
Ricker.model(nzero = 0.1, r = 1, time = 100)
