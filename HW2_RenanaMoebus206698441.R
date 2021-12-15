# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### Libraries ###

library(ggplot2)
library(dplyr)
library(ggseqlogo)
library(RColorBrewer)
library(ComplexHeatmap)
library(stats)
library(plyr)
library(e1071)
library(gridExtra)

### Functions ###

## Function - Calculation of p-values (exact, mid and approximation) ##
pValCalc <- function(x, n, p){
  
  if (x > p*n){
    
    # Exact p-value test
    exactPval <- 2*(pbinom(x-1, n, p, lower.tail=FALSE))
    
    # Mid p-value test
    midPval <- 2*(dbinom(x, n, p)/2 + pbinom(x, n, p, lower.tail=FALSE))
    
    # Binomial approximation test
    biPval <- pnorm(x+0.5, mean = n*p, sd = sqrt(p*(1-p)*n)) - pnorm(x-0.5, mean = n*p, sd = sqrt(p*(1-p)*n))
    
  } else {
    
    # Exact p-value test
    exactPval <- 2*(pbinom(x-1, n, p, lower.tail=TRUE))
    
    # Mid p-value test
    midPval <- 2*(dbinom(x, n, p)/2 + pbinom(x, n, p, lower.tail=TRUE))
    
    # Binomial approximation test
    biPval <- pnorm(x+0.5, mean = n*p, sd = sqrt(p*(1-p)*n)) - pnorm(x-0.5, mean = n*p, sd = sqrt(p*(1-p)*n))
    
  }
  
  return(c(exactPvalue = exactPval, midPvalue = midPval, biPvalue = biPval)) 
  
}


### Question 1 ###

# Create vectors for each category
S <- c(50, 176*4/16) # Stripes 
D <- c(41, 176*3/16) # Dots
SAD <- c(85, 176*9/16) # Stripes and dots

# Add vectors to dataframe
data1 <- data.frame(S, D, SAD)

# Calculate Chi-Square test
chisq.test(data1)



### Question 2 ###
Data <-read.csv('HW2_8h.csv')

# Check normality of the data with shapiro test
shapiro.test(Data$tlv5)
shapiro.test(Data$tlv6)
shapiro.test(Data$rg5)
shapiro.test(Data$rg6)

## Section 1 ## - Perform t test to check difference between cities in 2020
wilcox.test(Data$tlv5, Data$rg5, conf.level = 0.99)

## Section 2 ## - Perform t test to check difference between cities in 2021
wilcox.test(Data$tlv6, Data$rg6, conf.level = 0.99)

## Section 3 ## - Perform t test to check difference between RG and the population
Data$diffRG <- Data$rg6 - Data$rg5 # Calculate difference between the years in RG
wilcox.test(Data$diffRG, mu = 0.03, conf.level = 0.99, exact=FALSE) # Perform t-test

## Section 4 ## - Perform t test to check difference between TLV and the population
Data$diffTLV <- Data$tlv6 - Data$tlv5 # Calculate difference between the years in TLV
wilcox.test(Data$diffTLV, mu = 0.03, conf.level = 0.99, exact=FALSE) # Perform t-test

## Section 5 ## - Perform t test to check difference between the cities (height difference)
wilcox.test(Data$diffTLV, Data$diffRG, conf.level = 0.99, exact=FALSE)



### Question 3 ###


## Section 1 ### - Fair coin toss --> 3-10 tosses 

# Generate number of coin tosses
tossNum <- floor(runif(1, min = 3, max = 11))

# Define fair p = 0.5
pF <- 0.5

# Generate coin tosses
suc <- rbinom(1, tossNum, pF)

# Calculate p-values with 'pValCalc' function
resultsF <- pValCalc(suc, tossNum, pF)
print(resultsF)


## Section 2 ### - Unfair coin toss --> 3-10 tosses 

# Define unfair p = 0.8
pNF <- 0.8

# Generate coin toss
suc <- rbinom(1, tossNum, pNF)

# Calculate p-values with 'pValCalc' function
resultsNF <- pValCalc(suc, tossNum, pF)
print(resultsNF)


## Section 3 ## - Fair and unfair coin toss --> 3-10 tosses and 10,000 repetitions

# Define number of repetitions
rep <- 10000

# Define matrices for the fair and unfair tosses to save p-Values
fairMat <- matrix(0, nrow = rep, ncol = 3)
unfairMat <- matrix(0, nrow = rep, ncol = 3)

# Get p-values 10000 times
for (i in 1:rep){
  
  # Generate coin tosses fair and unfair
  sucF <- rbinom(1, tossNum, pF)
  sucNF <- rbinom(1, tossNum, pNF)
  
  # Calculate p-values and add to matrices
  fairMat[i,] <- pValCalc(sucF, tossNum, pF)
  unfairMat[i,] <- pValCalc(sucNF, tossNum, pF)
}


## Section 4 ## - Plot ROC curve for the data from the matrices

FP <- c(mapply(function(x){sum(fairMat[,1] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(fairMat[,2] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(fairMat[,3] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)))
TP <- c(mapply(function(x){sum(unfairMat[,1] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(unfairMat[,2] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(unfairMat[,3] < x)/10000}, x = seq(from = 0, to = 1, by = 0.001)))
Test <- rep(c('Exact p-value', 'Mid p-value', 'Normally approximated p-value'), each = 1001)
xAlpha <- 0.05
yAlpha <- 0.05
df <- data.frame(Test, FP, TP)

png(filename = 'ROC - pValues.png')
ggplot(data = df, aes(x = FP, y = TP, group = Test)) +
  geom_line(aes(color = Test))
dev.off()


## Section 5 ## - Add '+' where alpha = 0.05


## Section 6 ## - Two people toss fair and unfair coin toss --> 5000 tosses

# Generate coin toss - Person 1
suc_fair_P1 <- rbinom(5000, 1, pF)
suc_unfair_P1 <- rbinom(5000, 1, pNF)

# Generate coin toss - Person 1
suc_fair_P2 <- rbinom(5000, 1, pF)
suc_unfair_P2 <- rbinom(5000, 1, pNF)

## Function for normal approximation ##

normalApp <- function(x, n, p){
  return(pnorm(x+0.5, mean = n*p, sd = sqrt(p*(1-p)*n)) - pnorm(x-0.5, mean = n*p, sd = sqrt(p*(1-p)*n)))
}

# Create empty matrices for person 1 and person 2
matFair <- matrix(0, nrow = 5000, ncol = 2)
matUnfair <- matrix(0, nrow = 5000, ncol = 2)

# Normal approximation test - Fair coin - Person 1
matFair[,1] <- sapply(1:5000, function(x){normalApp(sum(suc_fair_P1[1:x]), x, pF)})

# Normal approximation test - Fair coin - Person 1
matUnfair[,1] <- sapply(1:5000, function(x){normalApp(sum(suc_unfair_P1[1:x]), x, pF)})

# Normal approximation test - Fair coin - Person 1
matFair[,2] <- sapply(1:5000, function(x){normalApp(sum(suc_fair_P2[1:x]), x, pF)})

# Normal approximation test - Fair coin - Person 1
matUnfair[,2] <- sapply(1:5000, function(x){normalApp(sum(suc_unfair_P2[1:x]), x, pF)})

# Get total p-Values
pVal_fair_P1 <-sum(matFair[,1])/5000
pVal_fair_P2 <-sum(matFair[,2])/5000
pVal_unfair_P1 <-sum(matUnfair[,1])/5000
pVal_unfair_P2 <-sum(matUnfair[,2])/5000

# ROC curve

FP <- c(mapply(function(x){sum(matFair[,1] < x)/5000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(matFair[,2] < x)/5000}, x = seq(from = 0, to = 1, by = 0.001)))
TP <- c(mapply(function(x){sum(matUnfair[,1] < x)/5000}, x = seq(from = 0, to = 1, by = 0.001)), mapply(function(x){sum(matUnfair[,2] < x)/5000}, x = seq(from = 0, to = 1, by = 0.001)))
Person <- rep(c('Person 1', 'Person 2'), each = 1001)
df <- data.frame(Person, FP, TP)

png(filename = 'ROC - persons.png')
ggplot(data = df, aes(x = FP, y = TP, group = Person)) +
  geom_line(aes(color = Person)) 
  # geom_point(data = hightlight_df, aes(x = FP, y = TP, group = Person))
dev.off()


