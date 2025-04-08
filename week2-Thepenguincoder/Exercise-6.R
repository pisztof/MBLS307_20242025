#1
area_calculator <- function(diameter) return(pi*(0.5*diameter)*(0.5*diameter))
test <- area_calculator(c(3.4, 5.0 ,6.0, 7,8))
print(test)

#2
degrees_converter <- function(farenheit){
  celcius = (farenheit-32)*(5/9)
  print(paste(toString(farenheit), "Farenheit =", toString(round(celcius,1)), "Celsius"))
}
degrees_converter(3)

#3
norm <- rnorm(100, 35, 15)
information_maker <- function(info){
  mean = sum(info)/length(info)
  if (length(info)%%2 == 0){
    median = (info[length(info)/2]+info[length(info)/2+1])/2
  } else{
    median = info[(length(info)+1)/2]
  }
  range = c(min(info), max(info))
  
  hist(info,freq=FALSE)
  lines(density(info))
  print(c(mean,median,range))
  return(c(mean,median,range))
}
test1<-information_maker(norm)

#4
#lol I already did this in the last assignment
values <- rnorm(100,50,25)
median2 <- function(info){
  if (length(info)%%2 == 0){
    median = (info[length(info)/2]+info[length(info)/2+1])/2
  } else{
    median = info[(length(info)+1)/2]
  }
  print(median)
  return(median)
}

wow <- median2(values)

#5
ricker_model <- function(N0, r, K=100, tend){
  t=1
  Ns = list(N0)
  while(t< tend){
    Ns[t+1] = Ns[t]
    Ns[t+1] = Ns[[t]]*exp(r*(1-(Ns[[t]]/K)))
    t=t+1
  }
  plot(seq(1,tend,1),Ns)
  lines(seq(1,tend,1),Ns)
}

ricker_model(10,1.2,100,100)



