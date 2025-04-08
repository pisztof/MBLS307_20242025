#Exercise 6

#6.1
area_func <- function(d) pi*(d/2)^2
area_func(3.4)

#6.2
fah_to_cen <- function(temp) paste("Fahrenheit:",temp,"is equivalent to",round((temp-32)*(5/9),2),"centigrade.")
fah_to_cen(0)

#6.3
norm_data <- rnorm(100, 35, 15)
summary_func <- function(data) {
  cat(paste("Mean:",mean(data),"\nMedian:",median(data),"\nRange:",max(data)-min(data)))
  hist(data,freq=FALSE)
  lines(density(data))
}
summary_func(norm_data)

#6.4
median_func <- function(data) {
  data <- sort(data)
  if (length(data) %% 2 == 0) {
    return((data[length(data)/2]+data[(length(data)/2)+1])/2)
  }
  else {
    return(data[(length(data)+1)/2])
  }
}

data1 <- c(1,3,2,4)
data2 <- c(1,3,5,4,2)
median_func(data1)
median_func(data2)

#6.5
ricker <- function(nzero, r, time, K=100) {
  data <- c(nzero)
  for (i in seq(time)) {
    data <- c(data, data[i]*exp(r*(1-(data[i]/K))))
  }
  plot(seq(0, time), data, type="l", xlab="Time", ylab="Population size")
}
ricker(10,2,100)