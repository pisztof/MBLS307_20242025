###Exercise 5:
###14/02/2025


###Question 1: importing the file and assigning a name
prawn <- read.table(file = 'data/prawnGR.CSV', header = TRUE, sep = ",", stringsAsFactors = TRUE, dec = ",")

##To look at the data
str(prawn)
summary(prawn)

prawn$GRate <- as.numeric(prawn$GRate)

#Boxplot to see growth rate Vs diet
boxplot(GRate ~ diet, data = prawn, xlab = "Diet", ylab = "Growth Rate")

#To identify observations in each diet treatment?
ncol (prawn)
nrow(prawn)

###Question 2: To assess the normality and equal variance assumptions

##  Shapiro-Wilk normality test
shapiro.test(prawn$GRate)
##W = 0.95523, p-value = 0.02761

shapiro.test(prawn$GRate[prawn$diet == "Artificial"])
##W = 0.97929, p-value = 0.8064

shapiro.test(prawn$GRate[prawn$diet == "Natural"])
##W = 0.91742, p-value = 0.023

##data seems to be normally distributed

##to check for variance
var.test(prawn$GRate ~ prawn$diet)
##F = 0.64466, num df = 29, denom df = 29, p-value = 0.2431
##alternative hypothesis: true ratio of variances is not equal to 1
##95 percent confidence interval:
## 0.3068334 1.3544189
##sample estimates:
##ratio of variances 
## 0.6446557 
##As the p value is greater than 0.05, there is no evidence to suggest that the variances are unequal. No evidence to reject null hypothesis

##Question 3: two sample t-test 
t.test(prawn$GRate ~ prawn$diet)
##t = 0, df = 55.413, p-value = 1
##alternative hypothesis: true difference in means between group Artificial and group Natural is not equal to 0
##95 percent confidence interval:
##-9.112781  9.112781
##sample estimates:
##mean in group Artificial    mean in group Natural 
###30.5                     30.5 

###to perform the t-test assuming equal variances
t.test(prawn$GRate ~ prawn$diet, var.equal = TRUE)
##t = 0, df = 58, p-value = 1
##alternative hypothesis: true difference in means between group Artificial and group Natural is not equal to 0
##95 percent confidence interval:
##-9.103724  9.103724
##sample estimates:
##mean in group Artificial    mean in group Natural 
##30.5                     30.5 
# No evidence to reject the Null hypothesis, therefore no difference in growth rate of prawns fed on either artificial or natural diet.

###Question 4: 
## To fit a linear model with GRate as the response variable and diet as an explanatory variable
GRate ~ diet
growth.lm < - Grate ~ diet
growth.lm <- lm(GRate ~ diet, data = prawn)

###Question 5:
###ANOVA table for the fitted model
anova(growth.lm)
##Analysis of Variance Table
##Response: GRate
##Df Sum Sq Mean Sq F value Pr(>F)
##diet       1      0    0.00       0      1
##Residuals 58  17995  310.26       
## p-value is same as in t-test 

###Question 6:
##Splitting the plotting device into 2 rows and 2 columns
par(mfrow=c(2,2))
(plot(growth.lm))
##plot the residuals to assess normality and equal variance 
## divide the plotting device into 2 rows and 2 columns to get all the graphs on one device


###Question 7:
##importing the file and assigning a name
gigartina <- read.table(file = 'data/Gigartina.CSV', header = TRUE, sep = ",", stringsAsFactors = TRUE, dec = ",")
str(
  gigartina
  )
table(
  gigartina$diatom.treat
  )
xtabs(
  ~diatom.treat, 
  data = gigartina
  )

# boxplot of the data
boxplot(
  gigartina$diameter ~ gigartina$diatom.treat,
  data = gigartina,
  ylab = "diameter",
  xlab = "diatom.treat"
  )

# boxplot of the data
boxplot(
  gigartina$diameter ~ gigartina$diatom.treat,
  data = gigartina,
  ylab = "diatom.treat",
  xlab = "diameter"
)

###Question 8: 
##What is your null hypothesis?
# The null hypothesis is that: there is no difference in mean diameter between the different treatment groups

###Question 9:
###conduct the ANOVA using the linear model function 
gigartina.lm <-lm (diameter ~ diatom.treat, data = gigartina)

###Question 10: Using anova
anova(gigartina.lm) 
##Analysis of Variance Table

##Response: diameter
##Df Sum Sq Mean Sq F value    Pr(>F)    
##diatom.treat  3 1880.3  626.76  22.775 1.929e-08 ***
##Residuals    36  990.7   27.52                      
# reject the null hypothesis, therefore there is a significant difference in the diameter 
# between the treatment groups (F 3 and 36 = 22.78, p < 0.001)

###Question 11:
# Splitting the plotting device into 2 rows and 2 columns
par(mfrow = c(2, 2))

# Plotting the residuals for the gigartina.lm model
plot(gigartina.lm)
# The plots look fine

###Question 13:
##installed package mosaic
install.packages("mosaic")
library(mosaic)
##use the function TukeyHSD()
TukeyHSD(gigartina.lm)  

##Tukey multiple comparisons of means
##95% family-wise confidence level

##Fit: aov(formula = x)

##$diatom.treat
            ##diff         lwr       upr     p adj
##Sdecl-ASGM  -14.3 -20.6184102  -7.98159 0.0000030
##Sexpo-ASGM   -8.9 -15.2184102  -2.58159 0.0029489
##Sstat-ASGM  -18.3 -24.6184102 -11.98159 0.0000000
##Sexpo-Sdecl   5.4  -0.9184102  11.71841 0.1165421
##Sstat-Sdecl  -4.0 -10.3184102   2.31841 0.3360087
##Sstat-Sexpo  -9.4 -15.7184102  -3.08159 0.0016145

####Sexpo-Sdecl and Sstat-Sdecl are significantly different

### Question 14:
plot(TukeyHSD(gigartina.lm))

### Question 15:
##importing the file and assigning a name
temora <- read.table(file = 'data/TemoraBR.CSV', header = TRUE, sep = ",", stringsAsFactors = TRUE, dec = ",")
str(temora)
## 45 observations and 3 variables
## type of variables: temperature =numeric, beat_rate = numeric , acclimitisation_temp = integer
## response (dependent) variable = heat_beat , explanatory (independent) variables = acclimitisation temp


## Questionn 16:
temora$Facclimitisation_temp <- factor(temora$acclimitisation_temp)
##boxplot to visualize
boxplot(
  temora$acclimitisation_temp ~ temora$beat_rate,
  data = temora,
  ylab = "acclimitisation_temp",
  xlab = "beat_rate"
)

# scatter plot using the with function
with(temora, plot(beat_rate ~ temp, xlab = "temperature", ylab = "beat rate"))

# using a coplot
coplot(beat_rate ~ temp | Facclimitisation_temp, data = temora)

###Question 17:
##analyse these data using an Analysis of Covariance (ANCOVA) 
## the slope between beat rate and temp looks different for each acclimatization temp

###Question 18:
##fit the model
temora$beat_rate <- as.numeric((temora$beat_rate))
temora$temp <- as.numeric((temora$temp))

temora.lm <- lm(beat_rate ~ temp + Facclimitisation_temp + temp:Facclimitisation_temp, data = temora)
temora.lm <- lm(beat_rate ~ temp * Facclimitisation_temp, data = temora)

###Question 19:
##summary using ANOVA
summary(temora.lm)
anova(temora.lm)
##Analysis of Variance Table

##Response: beat_rate
##                            Df Sum Sq Mean Sq F value Pr(>F)
##temp                        1  272.0  272.04  1.6048 0.2127
##Facclimitisation_temp       2  371.8  185.91  1.0967 0.3440
##temp:Facclimitisation_temp  2  334.9  167.46  0.9879 0.3815
##Residuals                  39 6611.2  169.52       

# there is a significant interaction between temp and Facclimitisation_temp therefore there is a significant 
# relationship between beat_rate and temp, and this relationship 
# is different depending on the level of Facclimitisation_temp.
# Therefore we should not interpret the main effect of temp or Facclimitisation_temp

###Question 20:Assess the assumptions of normality and equal variance
par(mfrow = c(2,2))
plot(temora.lm)
# the variance increases with the fitted values.

