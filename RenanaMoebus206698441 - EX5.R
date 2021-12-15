# Set working directory
#setwd("~/R/RenanaMoebus206698441 - EX5")

# Load libraries
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(Metrics)

### Question 2 ###

# Load the data
DataQue2 <- data.frame(Haemoglobin=unlist(strsplit('11.10 10.70 12.40 14.00 10.5 9.6 13.10 12.50 13.50 13.90 15.10 13.90 16.20 16.30 16.80 17.10 17.00 16.90 15.70 16.50',' ')), Age=unlist(strsplit('20 22 25 28 28 31 32 35 38 40 45 49 54 55 57 60 62 63 65 67',' ')))

# Remove women at age 31 and 32
indices <- c(which(DataQue2$Age == '31'), which(DataQue2$Age == '32'))
testSet <- DataQue2[indices,] # Data including only the women at age 31 and 32
trainSet <- DataQue2[-indices,] # Data excluding the women at age 31 and 32

## Section 2 ##

# Create a linear model where Haemoglobin is the response and Age is the explanatory
Model <- lm(formula = as.numeric(Haemoglobin) ~ as.numeric(Age), data = trainSet)

## Section 3 and 4 ##

# Print summary of the linear model for t-test and ANOVA test
summary(Model)

## Section 5 ##

# Predict hemoglobin level fo age 31 using 'prediction interval'
predict(Model, testSet[1,], interval = "prediction")

# Predict hemoglobin level fo age 31 using 'prediction interval'
predict(Model, testSet[2,], interval = "confidence")




### Question 3 ###

## Section 1 ##

# Set seed to 42
set.seed(42)

## Section 2 ##

# Load data from the file
Data <- read.csv(file = 'HW5.csv')

## Section 3 ##

# Get number of observations
rowNum <- nrow(Data)

# Get training set and test set
trainSet <- Data[1:(rowNum/2),]
testSet <- Data[((rowNum/2)+1):rowNum,]

## Section 4 ##

# Fit the two linear regression models
Model_1 <- lm(y ~ a + b + c + d, data = trainSet) 
Model_2 <- lm(y ~ a + b + c + d + I(a^2) + I(b^2) + I(c^2), data = trainSet)

## Section 5 ##

# Find number of parameters for each model
parsNum_1 <- Model_1$rank - 1
parsNum_2 <- Model_2$rank - 1

# Calculate RMSE for training set for both models
trainRMSE_1 <- rmse(trainSet$y, predict(Model_1, trainSet))
trainRMSE_2 <- rmse(trainSet$y, predict(Model_2, trainSet))

# Calculate RMSE for test set for both models
testRMSE_1 <- rmse(testSet$y, predict(Model_1, testSet))
testRMSE_2 <- rmse(testSet$y, predict(Model_2, testSet))

## Section 7 ##

# Define vector for RMSE results
ttestRMSE_1 <- matrix(0, nrow = 1, ncol = 10)
ttestRMSE_2 <- matrix(0, nrow = 1, ncol = 10)

# Loop to get 30 RMSE results from 30 trained models
for (i in 1:10){
  
  # Split data to training and test set
  indices <- sample(1:1000, 500, replace = F)
  trainSet <- Data[indices,]
  testSet <- Data[-indices,]
  
  # Train models on training set
  Model_1 <- lm(y ~ a + b + c + d, data = trainSet) 
  Model_2 <- lm(y ~ a + b + c + d + I(a^2) + I(b^2) + I(c^2), data = trainSet)
  
  # Calculate RMSE for test set for both models
  ttestRMSE_1[i] <- rmse(testSet$y, predict(Model_1, testSet))
  ttestRMSE_2[i] <- rmse(testSet$y, predict(Model_2, testSet))
}

# Run t-test on RMSE results for both models
pValModels <- t.test(ttestRMSE_1, ttestRMSE_2, mu = 0, alt = 'two.sided', paired = T, conf.level = 0.99)

# Print the p-value got from the t-test
print(pValModels$p.value)