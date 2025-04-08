#importing the data
data <- read.table("C:\\Users\\sandy\\Downloads\\whaledata.txt", header= TRUE, sep= "\t")
data
head(data)
names(data)
str(data)
#month is a character and water.noise is also the same
summary(data)
whale.sub <- data[1:10, 1:4]
whale.sub
whale.num <- data[,c(1,3,4)]
whale.num
rows_all_coloumns <- data[c(1:50),]
rows_all_coloumns
rows_coloumns <- data[c(11:101),c(1:7)]
rows_coloumns
#9
data_depth <- data[depth> 1200]
data_depth
