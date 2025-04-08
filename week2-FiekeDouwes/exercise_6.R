#Exercise 6
#https://alexd106.github.io/intro2R/exercise_6.html

# Fri Feb 14 13:45:06 2025 ------------------------------

#1
circle_area_func <- function(r) pi*r**2
circle_diam <- 3.4
circle_area_func(circle_diam/2)
circle_diams <- 1:4
circle_area_func(circle_diams/2)

#2
F_to_C_converter <- function(Fahr) {
  C <- (Fahr-32)*5/9
  sprintf("Fahrenheit: %s is equivalent to %s centigrade.", Fahr, C)
}
F_to_C_converter(-40)

#3
norm_vec = rnorm(100, mean = 35, sd = 15)
vec_stats = function(my_vector) {
  vec_mean <- mean(my_vector, na.rm = TRUE)
  vec_median <- median(my_vector, na.rm = TRUE)
  vec_range <- c(min(my_vector, na.rm = TRUE), max(my_vector, na.rm = TRUE))
  sprintf("Mean: %s, Median: %s, range: (%s, %s)", 
          vec_mean, vec_median, vec_range[1], vec_range[2])
  hist(my_vector, freq = FALSE)
  lines(density(my_vector))
}
vec_stats(norm_vec)

#4
my_median <- function(my_vector){
  sorted_vector <- sort(my_vector)
  if (length(sorted_vector)%%2 == 1)
    return(sorted_vector[length(sorted_vector) %/% 2 + 1])
  else
    middle <- length(sorted_vector) / 2
    return(mean(sorted_vector[middle:(middle + 1)]))
}
my_median(c(1,7,6,3,5,4,2))
my_median(c(8,3,5,2,1,7,6,4))

#5
ricker <- function(nzero, r, time, K = 100) {
  n <- c(nzero)
  for (i in 1:time)
    n[i+1] <- round(n[i]*exp(r*(1-n[i]/K)))
  return(n)
}
ricker(200,2,5)
ricker(400,2,5)
ricker(50,2,5)
ricker(50,4,5)
