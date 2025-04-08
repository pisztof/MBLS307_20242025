#1 
area_circle <- function(d) {
  area <- pi*(d/2)**2
  return(area)
}

area_circle(d=3.4)

#with a vector
vector_d <- c(2,3,4,5)
area_circle(d=vector_d)

#2 
far_cen <- function(f) {
  cen <- (f-32)*5/9 
  return(paste("Farenheit:", f, "is equivalent to", round(cen,1), "centigrade."))
}
far_cen(f=87)
#also made it so the result is rounded up to 1 decimal place

?hist
#3
random <- rnorm(100, mean=35, sd=15)

statis <- function(v) {
  mean_v <- mean(v)
  median_v <- median(v)
  range_v <- range(v)
  hist(v)
  lines(v, col="red")
  return(paste("mean:", mean_v, "median:", median_v, "range:", range_v))
}

statis(random)

#4
#the %% means the reminder!!
median_fun <- function(nums) {
  if (length(nums) %% 2 == 1){
    pos1 <- nums[(length(nums)/2)-0.5]
    pos2 <- nums[(length(nums)/2)+0.5]
    result <- (pos1+pos2)/2
  }
  else {
    result <- nums[length(nums)/2]
  }
  return(result)
}

vector_odd <- c(1,2,3,4,5,6)
median_fun(vector_odd)

vector_even <- c(1,2,3,4,5)
median_fun(vector_even)

#5
rick_fun <- function(r,n0,k=100,t){
  n_vector<-t+1
  n_vector[1]<-n0
  for (i in 1:t) {
    n_vector[i+1] <- n_vector[i]*exp(r*(1-n_vector[i]/k))
  }
  t_vector <- 0:t
  plot(t_vector,n_vector, xlab="time", ylab="population size")
}

rick_fun(r=1, n0=0.1, t=100)
