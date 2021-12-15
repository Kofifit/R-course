### Question 1 - cherry blossom data set ###

# Load data from csv file
Data <- read.csv('C:\\Users\\moebu\\OneDrive\\Documents\\R\\Cherry-Blossom-2008.csv')

# Find observations number
RowNum <- nrow(Data)

# Find number of females
FNum = sum(Data$sex =='F')

# Find number of males
MNum = sum(Data$sex =='M')

# Find age range
MaxAge = max(Data$age)
MinAge = min(Data$age)

# Extract age and run time vectors for females and males separately 
AgeF <- Data[Data$sex == 'F', 5]
AgeM <- Data[Data$sex == 'M', 5]
GunF <- Data[Data$sex == 'F', 8]
GunM <- Data[Data$sex == 'M', 8]

# Calculate the mean and SD for age and run time for females and makes separately 
ageSummaryF <- c(mean(AgeF), sd(AgeF))
gunSummaryF <- c(mean(GunF), sd(GunF))
ageSummaryM <- c(mean(AgeM), sd(AgeM))
gunSummaryM <- c(mean(GunM), sd(GunM))

# Create a new data set consisting only females
FemalesOnlyData <- Data[Data$sex =='F',]

# Sort new data set by country and then position. The remove gender column
FemalesOnlyData <- FemalesOnlyData[with (FemalesOnlyData, order(place, position)), ]
FemalesOnlyData <- subset(FemalesOnlyData, select = -sex)

# Calculate distance in meters
dis <- 1609.34*10

# Add speed column and relocate it as the first column
FemalesOnlyData$speed = dis/(FemalesOnlyData$gun)
FemalesOnlyData <- subset(FemalesOnlyData, select = c(9 ,1:8))

# Calculate the average speed for Males from Kenya
AveSpeedMKenya <- dis/mean(Data[Data$sex == 'M'  & Data$place == 'Kenya',8])

# Sort data set by age
DatabyAge <- Data[with (Data, order(age)), ]

# Calculate average speed of 10 youngest participants
Ave10Youngest <- dis/mean(head(DatabyAge[,8], 10))

# Calculate average speed of 10 oldest participants
Ave10Oldest <- dis/mean(tail(DatabyAge[,8], 10))




### Question 2 - mtcars data set ###

# Attach mtcars data set
attach(mtcars)

# Find number of observation and number of variables
observationNum <- nrow(mtcars)
variableNum <- ncol(mtcars)

# Create new data frame with gear, mpg, am, and disp
newData <- data.frame(rownames(mtcars), mtcars$gear, mtcars$mpg, mtcars$am, mtcars$disp)

# Sort data by gear and then mpg
newData <- newData[with (newData, order(mtcars$gear, mtcars$mpg)), ]

# Find name of first car in data set
firstCarName <- newData[1,1]

# Save new data frame to csv file
pathName <- sprintf('C:\\Users\\moebu\\OneDrive\\Documents\\R\\%s_206698441.csv', firstCarName)
write.csv(newData, pathName, row.names = FALSE)



### Question 3 - msleep data set %%%

# load and attach msleep data set
library(ggplot2)
attach(msleep)

# open pdf file
pdf("Number of observations based on diet type plot.pdf") 

# plot the bar graph for diet type
ggplot(msleep, aes(x = as.factor(vore), fill = as.factor(vore) )) + 
  geom_bar(xlab = 'Diet type', ylab = 'Number of observations' ) +
  scale_fill_brewer(palette = "Purples") +
  labs(title = 'Number of observations based on diet type', x = 'Diet type', y = 'Number of observations', fill = 'Diet type') +
  theme(legend.position = "left")

# Close the pdf file
dev.off() 

# open pdf file
pdf('Total time of sleep based on weight of the brain plot.pdf')

# plot the scatter graph for brain weight and total sleep time
ggplot(msleep, aes(x = brainwt, y = sleep_total)) + 
  geom_point(shape = 2, fill = 'green', color = 'green', alpha = 0.7) +
  labs(title = 'Total time of sleep based on weight of the brain', x = 'Brain weight', y = 'Total time of sleep (hours)')

# close the pdf file
dev.off()