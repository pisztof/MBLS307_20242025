#1
circle.area <- function(r){
  return(pi*r*r)
}
circle.area(3.4)
circle.area(c(1,2,3,4))

#2
FahrToCels <- function(F){
  cels <- toString((F-32)*5/9)
  print(paste('Farenheid : value of oF is equivalent to ',cels,' oC centrigrade'))
}
FahrToCels(5)

#3
info.distr <- function(distr){
  print(paste('Mean: ', mean(distr)))
  print(paste('Median: ', median(distr)))
  print(paste('Range: ', min(distr), '-', max(distr)))
  hist(distr, probability=TRUE)
  d <- density(distr)
  lines(d$x,d$y,col='red',lwd=2)
}
distr.norm <- rnorm(100,mean=35,sd=15)
info.distr(distr.norm)

#4 (median without if statement)
median.own <- function(vect){
  vect.sorted <- sort(vect)
  med1 = ((length(vect)+1)+(length(vect)-1)%%2)/2
  med2 = ((length(vect)+1)-(length(vect)-1)%%2)/2
  medi = (vect.sorted[med1]+vect.sorted[med2]) / 2
  return(medi)
}

#5
ricker <- function(nzero, r, t, K=100){
  tlist = c(0:t)
  Nlist = c(nzero)
  for (i in tlist[1:t]){
    Nlist <- append(Nlist, tail(Nlist,1)* (r * (1-tail(Nlist,1)/K)))
  }
  return(data.frame(tlist, Nlist))
}
plot(ricker(5,20/19,100))
