print("hellooo")

#object, assignment/gets operator, value
my_obj <- 4555
my_obj2 <- 1000
my_obj3 <- my_obj + my_obj2

#vector(data structure), concatenate, series of values
my_vec <- c(3,4,5)
my_vec2 <- c(3,5,9)
summary(my_vec)
my_seq <- 10:1 #vector 10 to 1
my_seq2 <- seq(from = 20, to = 5, by = -0.7)
my_seq3 <- rep(4:8, each = 4)
my_vec3 <- my_vec * my_vec2
my_vec3

#extract data
my_value <- my_seq3[c(3,8)]
my_val2 <- my_seq3[my_seq3 > 6| my_seq3 < 5]
my_val2[my_val2 <= 4] <- 503
my_val2

#order data
sorted_vec <- rev(sort(my_val2,decreasing = TRUE))
group1 <- c("henk", "gert", "bob", "dirk")
grade1 <- c(10,4,6,8)
grade1_ord <- order(grade1)
group1_ord <- group1[grade1_ord]
group1_ord

#missing data
my_data <- c(50,NA,80,3,202,30)
mean(my_data,na.rm = TRUE)

#exercise 2
log(12.43)      #1`
log10(12.43) 
sqrt(12.43)
exp(12.43)
area_circle <- pi*10^2 #3
(14*0.51)^(1/3) #4
weight <- c(69, 62, 57, 59, 59, 64, 56, 66, 67, 66) #5
mean(weight) #6
var(weight)
sd(weight)
length(weight)
range(weight)
first_five <- weight[1:5]
height <- c(112, 102, 83, 84, 99, 90, 77, 112, 133, 112) #7
summary(height)
some_child <- height[c(2,3,9,10)]
shorter_child <- height[height <= 99]
bmi <- weight * (height/100)^2 #8
seq1 <- seq(from = 0, to = 1, by = 0.1) #9
seq2 <- seq(from = 10, to = 1, by = -0.1) #10
seq3 <- rep(1:3,times = 3)
seq4 <- rep(c("a","c","e","g"), each = 3)
seq5 <- rep(c("a","c","e","g"), times = 3)
seq6 <- rep(1:3,times = 2,each = 3)
vec10 <- c(5:1)
seq7 <- rep(1:5, times = vec10)
vec11 <- c(7,2,8,1)
mult <- c(4,3,1,5)
seq8 <- rep(vec11,times = mult)
height_sorted <- sort(height) #11
height_sorted2 <- rev(height_sorted)
child_name <- c("Alfred", "Barbara", "James", "Jane", "John", "Judy", "Louise", "Mary", "Ronald", "William")
height_order <- order(height) #14
names_sort = child_name[height_order]
weight_rev = child_name[rev(order(weight))]
weight_rev
print(paste("the heaviest is", weight_rev[1], ", and the lightest is ", weight_rev[10]))
mydata <- c(2, 4, 1, 6, 8, 5, NA, 4, 7) #16
mean(mydata, na.rm = TRUE)      
ls() #17
rm(seq1)

#data types
#numeric (decimal), integer (whole nr), logical (true or false or na), character ("string") 
new <- c(3.5,6.7,4,"hi")
class(new)
is.numeric(new)
new_char <- as.character(new)
new_num <- as.numeric(new_char)
new_num
#scalar is a vector with only one value

#matrix is 2D, array is multipleD
my_mat <- matrix(1:16, nrow = 4, byrow = TRUE)
my_ar <- array(1:16, dim = c(2,4,2))
colnames(my_ar) <- c("hi","bye","yes","no")
my_mat_t <- t(my_mat)
my_mat_diag <- diag(my_mat)
mat.1 <- matrix(1:4, nrow = 2)
mat.2 <- matrix(4:1, nrow = 2)
mat.1 * mat.2
mat.3 <- mat.1 %*% mat.2

#lists
my_list <- list(hellow = c(2,3,4),
              bye = c("hi, bye", "yes"),
              ciao = matrix(1:8,nrow =2))
names(my_list) <- c("djfd","jdkjfa", "fjdkaf")
my_list

#dataframe
#2D similar to matrix but can hold different types of data
p.height <- c(180, 155, 160, 167, 181)
p.weight <- c(65, 50, 52, 58, 70)
p.names <- c("Joanna", "Charlotte", "Helen", "Karen", "Amy")
dataf <- data.frame(height = p.height, weight = p.weight, names = p.names,stringsAsFactors = TRUE)
dim(dataf)
str(dataf)
rownames(dataf) <- p.names
dataf

#importing data
flowers <- read.table(file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/flower.txt', header = TRUE, sep = "\t",
                      stringsAsFactors = TRUE)
install.packages("data.table")
library(data.table)
all_data <- fread(file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/flower.txt')
install.packages('tinyverse')

#wrangling data
summary(flowers$height)
flowers[1:10,c(1,3,4)]
flowers$height[1]
flowers[1:10,]  #[row,collumn] -> empty means all of them
flowers[1:10,-1:-3] #can exclude values using - (also with -c(1,2,3))
flowers[1:5,c('weight',"shootarea")]
big_flowers <- flowers[flowers$height > 12,]
nit_high <- flowers[flowers$nitrogen != "high",1:4]
low_notip_heigh6 <- flowers[flowers$height >= 6 & flowers$nitrogen == "medium" &
                              flowers$treat == "notip", ]
tip_med_2 <- subset(flowers, treat == "tip" & nitrogen == "medium" & block == 2,
                    select = c("treat", "nitrogen", "leafarea"))

#order data
ordered_height <- flowers[order(flowers$block,-flowers$height),] #order based on two variables
flowers$nitrogen <- factor(flowers$nitrogen, 
                           levels = c("low", "medium", "high"))  #set order of factors  
nit_ord <- flowers[order(flowers$nitrogen),]

#add to data
df1 <- data.frame(id = 1:4, height = c(120, 150, 132, 122),
                  weight = c(44, 56, 49, 45))
df2 <- data.frame(id = 5:6, height = c(119, 110),
                  weight = c(39, 35))
df3 <- data.frame(id = 1:4, height = c(120, 150, 132, 122),
                  weight = c(44, 56, 49, 45))
df4 <- data.frame(location = c("UK", "CZ", "CZ", "UK"))
df_rcomb <- rbind(df1, df2)
df_ccomb <- cbind(df3,df4)
df_rcomb$height_log10 <- log10(df_rcomb$height)
df_rcomb$Fid <- factor(df_rcomb$id)
str(df_rcomb)

#merge data
taxa <- data.frame(GENUS = c("Patella", "Littorina", "Halichondria", "Semibalanus"),
                   species = c("vulgata", "littoria", "panacea", "balanoides"),
                   family = c("patellidae", "Littorinidae", "Halichondriidae", "Archaeobalanidae"))
zone <- data.frame(genus = c("Laminaria", "Halichondria", "Xanthoria", "Littorina", 
                             "Semibalanus", "Fucus"),
                   species = c("digitata", "panacea", "parietina", "littoria", 
                               "balanoides", "serratus"),
                   zone = c( "v_low", "low", "v_high", "low_mid", "high", "low_mid"))
taxa_zone <- merge(x = taxa, y = zone, by.x = c("species", "GENUS"), 
                   by.y = c("species", "genus"), all = TRUE)
long_data <- data.frame(
  subject = rep(c("A", "B", "C", "D"), each = 3),
  sex = rep(c("M", "F", "F", "M"), each =3),
  condition = rep(c("control", "cond1", "cond2"), times = 4),
  measurement = c(12.9, 14.2, 8.7, 5.2, 12.6, 10.1, 8.9,
                  12.1, 14.2, 10.5, 12.9, 11.9))
wide_data <- data.frame(subject = c("A", "B", "C", "D"),
                        sex = c("M", "F", "F", "M"),
                        control = c(12.9, 5.2, 8.9, 10.5),
                        cond1 = c(14.2, 12.6, 12.1, 12.9),
                        cond2 = c(8.7, 10.1, 14.2, 11.9))

#reshape data
install.packages("reshape2")
library(reshape2)
my_long_df <- melt(data = wide_data, id.vars = c("subject", "sex"),
                   measured.vars = c("control", "cond1", "cond2"),
                   variable.name = "condition", value.name = "measurement")
my_long_df
my_wide_df <- dcast(data = long_data, subject + sex ~ condition, 
                    value.var = "measurement")

#exercise 3
whale <- read.table(file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/data/whaledata.txt', header = TRUE)
str(whale)
# data.frame':	100 obs. of  8 variables:
#  $ month          : chr  "May" "May" "May" "May" ...
#  $ time.at.station: int  1344 1633 743 1050 1764 580 459 561 709 690 ...
#  $ water.noise    : chr  "low" "medium" "medium" "medium" ...
#  $ number.whales  : int  7 13 12 10 12 10 5 8 11 12 ...
#  $ latitude       : num  60.4 60.4 60.5 60.3 60.4 ...
#  $ longitude      : num  -4.18 -4.19 -4.62 -4.35 -5.2 -5.22 -5.08 -5 -4.64 -4.84 ...
#  $ depth          : int  520 559 1006 540 1000 1000 993 988 954 984 ...
#  $ gradient       : int  415 405 88 409 97 173 162 162 245 161 ...
summary(whale) #missing one data point in number.whales
whale.sub <- whale[1:10,1:4]
whale.num <- whale[,c("month","water.noise","number.whales")]
whale.may <- whale[1:50,]
whale.some <- whale[-1:-10,-8]
whale.deep <- whale[whale$depth >1200,] #9
whale.steep <- whale[whale$gradient >200,]
whale.lownoise <- whale[whale$water.noise== "low",]
whale.highmay <- whale[whale$water.noise == "high" & whale$month == "May",]
whale.oct <- whale[whale$month == "October" & whale$water.noise == "low" & whale$gradient >median(whale$gradient),]
whale.latlon <- whale[whale$latitude >=60.0 & whale$latitude <= 61.0 & whale$longitude >=-6 & whale$longitude <=-4,]
whale.nomed <- whale[whale$water.noise != "medium",]
whale.deepgreat <- whale[whale$depth > 1500 & whale$number.whales > mean(whale$number.whales, na.rm =TRUE), ]
whale.shortmay <- subset(whale, month == "May" & time.at.station < 1000 & depth > 1000,)
whale.lat <- subset(whale, month == "October" & latitude > 61, 
                    select = c("month","latitude","longitude","number.whales"))
whale.depth.sort <- whale[order(whale$depth),]
whale$water.noise <- factor(whale$water.noise,
                            levels = c('low','medium','high'))
whale.depth.sortinnoise <- whale[order(whale$water.noise, whale$depth),]
tapply(whale$number.whales, whale$water.noise, mean,na.rm = TRUE)
tapply(whale$number.whales, whale$water.noise, median, na.rm = TRUE)
tapply(whale$number.whales, whale$month, median, na.rm = TRUE)
aggregate(cbind(whale$time.at.station, whale$number.whales, whale$depth, whale$gradient), by = list(whale$water.noise), data = whale ,FUN = mean, na.rm = TRUE)
table(whale$water.noise,whale$month)
write.table(whale.num, file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/data/whalenum.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#graphics with base R
library(lattice)
plot(flowers$weight,flowers$shootarea)
my_x <- 1:10
my_y <- seq(from = 1, to = 20, by = 2)
par(mfrow = c(1, 1))
plot(my_x, my_y, type = "l")
plot(my_x, my_y, type = "b")
plot(my_x, my_y, type = "o")
plot(my_x, my_y, type = "c")

#histogram
hist(flowers$height)
brk <- seq(from = 0, to = 18, by = 1)
hist(flowers$height, breaks = brk, main = "petunia height",freq = FALSE)
dens <- density(flowers$height)
lines(dens)

#boxplot
boxplot(flowers$weight, ylab = "weight in g")
flowers$nitrogen <- factor(flowers$nitrogen, 
                           levels = c("low", "medium", "high"))
boxplot(weight ~ nitrogen, data = flowers, 
        ylab = "weight (g)", xlab = "nitrogen level")
boxplot(weight ~ nitrogen * treat, data = flowers, 
        ylab = "weight (g)", xlab = "nitrogen level",cex.axis = 0.7)

#violinplot
install.packages('vioplot')
library(vioplot)
vioplot(weight ~ nitrogen, data = flowers, 
        ylab = "weight (g)", xlab = "nitrogen level",
        col = "lightblue")

#dot charts
dotchart(flowers$height, groups = flowers$nitrogen)

#panels and correlation
pairs(flowers[, c("height", "weight", "leafarea", 
                  "shootarea", "flowers")], panel = panel.smooth)
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
pairs(flowers[, c("height", "weight", "leafarea", 
                  "shootarea", "flowers")], 
      lower.panel = panel.cor)
pairs(flowers[, c("height", "weight", "leafarea", 
                  "shootarea", "flowers")], 
      lower.panel = panel.cor,
      diag.panel = panel.hist,
      upper.panel = panel.smooth)

#exercise 4
squid1 <- read.table(file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/data/squid1.txt', header = TRUE)
str(squid1)
#'data.frame':	519 obs. of  13 variables:
# $ sample.no        : int  105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 105128901 ...
# $ specimen         : int  1002 1003 1005 1007 1008 1009 1011 1013 1014 1017 ...
# $ year             : int  1989 1989 1989 1989 1989 1989 1989 1989 1989 1989 ...
# $ month            : int  12 12 12 12 12 12 12 12 12 12 ...
# $ weight           : num  152 106 138 141 126 ...
# $ sex              : int  2 2 2 2 2 2 2 2 2 2 ...
# $ maturity.stage   : int  3 1 2 2 3 1 2 3 3 4 ...
# $ DML              : int  174 153 169 175 169 116 135 192 170 205 ...
# $ eviscerate.weight: num  87.5 62.6 79.4 83.1 72.2 ...
# $ dig.weight       : num  4.648 3.138 0.307 4.123 3.605 ...
# $ nid.length       : num  39.4 24.1 39 41.4 39.8 20 14 55 44 53 ...
# $ nid.weight       : num  2.46 0.319 1.169 1.631 2.03 ...
# $ ovary.weight     : num  1.68 0.103 0.289 0.252 0.86 ..
summary(squid1) 
#591 obersvations, 13 variables, 
squid1$Fyear <- factor(squid1$year)
squid1$Fmonth <- factor(squid1$month)
squid1$Fmaturity.stage <- factor(squid1$maturity.stage)
str(squid1)
table(squid1$year,squid1$month)
dotchart(squid1$DML)
dotchart(squid1$weight)
dotchart(squid1$nid.length)
pdf(file = 'C:/Users/Tjits/OneDrive/Bureaublad/Studie/Omics/data/dmldot.pdf')
dotchart(squid1$DML)
dev.off()
which(squid1$nid.length > 400)
squid1[11,11] <- 43.2
dotplot(squid1$weight~squid1$DML)


#ggplot
install.packages('tidyverse')
library(tidyverse)
squid1 %>% ggplot(aes(x = weight, y = nid.length)) +
  geom_smooth(method='lm')+
  geom_point(aes(colour = maturity.stage)) +
#  facet_wrap(~month)+
  xlim(0,800)+ ylim(0,150)+
  xlab('weight (g)') +
  ylab('length')+
  ggtitle(label = 'oh wow', subtitle = 'so pretty')

#density plot
ggplot(flowers) +
  geom_density(aes(x = flowers, fill = nitrogen), alpha = 0.5) +
  labs(y = "Density", x = "Number of Flowers", fill = "Nitrogen\nConcentration") +
  scale_fill_manual(labels = c("High", "Medium", "Low"),
                    values = c("#DB24BC", "#BCDB24", "#24BCDB"))

#histogram
ggplot(flowers) +
  geom_histogram(aes(x = flowers, fill = factor(block)), colour = "black", bins = 20,
                 position = "dodge") +
  labs(y = "Count", x = "Number of Flowers", fill = "Greenhouse\nBlock") +
  scale_fill_manual(labels = c("1", "2"),
                    values = c("#8CD926", "#7326D9"))

#frequency polygon
ggplot(flowers) +
  geom_freqpoly(aes(x = flowers, colour = factor(block)), size = 1, bins = 20) +
  labs(y = "Count", x = "Number of Flowers", colour = "Greenhouse\nBlock") +
  scale_colour_manual(labels = c("1", "2"),
                      values = c("#8CD926", "#7326D9"))

#boxplot
ggplot(flowers) +
  geom_boxplot(aes(y = flowers, x = nitrogen)) +
  labs(y = "Number of Flowers", x = "Nitrogen Concentration")

#violing plot
ggplot(flowers) +
  geom_violin(aes(y = flowers, x = nitrogen, fill = nitrogen), 
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(aes(y = flowers, x = nitrogen), colour = "black", height = 0, 
              width = 0.1, alpha = 0.5) +
  scale_fill_manual(labels = c("High", "Medium", "Low"),
                    values = c("#5f7f5c", "#749770", "#9eb69b")) +
  labs(y = "Number of Flowers", x = "Nitrogen Concentration") +
  scale_y_log10()

#barchart
ggplot(flowers) +
  geom_bar(aes(x = factor(block), fill = nitrogen)) +
  scale_fill_manual(labels = c("High", "Medium", "Low"),
                    values = c("#2613EC", "#9313EC", "#EC13D9")) +
  labs(y = "Count", x = "Block")

#add quantile lines
install.packages('quantreg')
library(quantreg)
ggplot(aes(x = weight, y = shootarea), data = flowers) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_quantile(colour = "darkgrey", size = 1) +
  labs(y = "Shoot Area", x = "Weight")

#heatmap
ggplot(aes(x = leafarea, y = shootarea), data = flowers) +
  geom_bin2d() +
  labs(y = "Shoot Area", x = "Leaf Area") +
  coord_cartesian(xlim = c(5,25))

#hexmap
install.packages('hexbin')
library(hexbin)
ggplot(aes(x = leafarea, y = shootarea), data = flowers) +
  geom_hex() +
  labs(y = "Shoot Area", x = "Leaf Area") +
  coord_cartesian(xlim = c(5,25))

#contourmap
ggplot(aes(x = leafarea, y = shootarea), data = flowers) +
  geom_density2d() +
  labs(y = "Shoot Area", x = "Leaf Area") +
  coord_cartesian(xlim = c(5,25))

ggplot(aes(x = leafarea, y = shootarea), data = flowers) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
  labs(y = "Shoot Area", x = "Leaf Area") +
  coord_cartesian(xlim = c(5,25))

#cleveland dotplot
ggplot(flowers) + 
  geom_hline(aes(yintercept = as.numeric(rownames(flowers))), linetype = 2,
             colour = "lightgrey") +
  geom_point(aes(x = leafarea, y = as.numeric(rownames(flowers)))) +
  labs(y = "Order in Dataset", x = bquote("Area of leaf"~(cm^2)))

#pairs plot
install.packages('GGally')
library(GGally)
flowers$block <- factor(flowers$block)

ggpairs(flowers, aes(colour = nitrogen), cardinality_threshold = NULL, progress = FALSE)


#exercise
mpg_df <- mpg
str(mpg_df)

ggplot(mpg_df) +
  geom_boxplot(aes(y = displ, x = class)) +
  labs(y = "displ", x = "class") +
  geom_jitter(aes(y = displ, x = class), colour = "black", height = 0, 
              width = 0.1, alpha = 0.5) +
  theme_classic()

#pipe?
#ctrl shift M
mpg_df %>% head()

#tidyverse
table1
table1 %>%
  mutate(rate = cases / population * 10000) %>% 
  group_by(year) %>% 
    summarize(total_cases = sum(cases))
ggplot(table1, aes(x = year, y = cases)) +
  geom_line(aes(group = country), color = "grey50") +
  geom_point(aes(color = country, shape = country)) +
  scale_x_continuous(breaks = c(1999, 2000)) # x-axis breaks at 1999 and 2000


#exercises:
table1
#each column is a variable, each row an observation
table2
#each column is a variable, each row is an observation per type
table3
#each column is a variable but rate has two, and each row is an observation

billboard
billboard %>% pivot_longer(
  cols = starts_with("wk"), 
  names_to = "week", 
  values_to = "rank",
  values_drop_na = TRUE
  )

billboard_longer <- billboard |> 
  pivot_longer(
    cols = starts_with("wk"), 
    names_to = "week", 
    values_to = "rank",
    values_drop_na = TRUE
  ) |> 
  mutate(
    week = parse_number(week)
  )

billboard_longer |> 
  ggplot(aes(x = week, y = rank, group = track)) + 
  geom_line(alpha = 0.25) + 
  scale_y_reverse()          

#pivoting
df <- tribble(
  ~id,  ~bp1, ~bp2,
  "A",  100,  120,
  "B",  140,  115,
  "C",  120,  125
)
df |> 
  pivot_longer(
    cols = bp1:bp2,
    names_to = "measurement",
    values_to = "value"
  )

#separate based on names?
who2 |> 
  pivot_longer(
    cols = !(country:year),
    names_to = c("diagnosis", "gender", "age"), 
    names_sep = "_",
    values_to = "count"
  )

household |> 
  pivot_longer(
    cols = !family, 
    names_to = c(".value", "child"), 
    names_sep = "_", 
    values_drop_na = TRUE
  )

#pivot wider
cms_patient_experience |> 
  distinct(measure_cd, measure_title)

cms_patient_experience |> 
  pivot_wider(
    id_cols = starts_with("org"),
    names_from = measure_cd,
    values_from = prf_rate
  )
df <- tribble(
  ~id, ~measurement, ~value,
  "A",        "bp1",    100,
  "B",        "bp1",    140,
  "B",        "bp2",    115, 
  "A",        "bp2",    120,
  "A",        "bp3",    105
)

df |> 
  pivot_wider(
    names_from = measurement,
    values_from = value
  )
    