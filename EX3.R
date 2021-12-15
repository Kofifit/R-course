### Libraries ###

library(magrittr)
library("plyr")
library(stringr)
library(ggplot2)

# Set working directory
setwd("~/R/RenanaMoebus206698441EX3")

### Functions ###

# Function to translate DNA into amino acids
TranslateDNA <- function(DNAvec, dictionary){
  
  proVec <- vector(mode = "list", length = 0) # Define empty list for protein
  index <- 0 # Define index for protein vector
  l <- nchar(DNAvec)
  for(i in seq(from=1, to=nchar(DNAvec)-2, by=3)){ # Loop to translate DNA to AA
    too <- nchar(DNAvec)-2
    index <- index + 1
    proVec[index] <- dictionary[[substr(DNAvec,i,i+2)]]# Find codon in dictionary and append to protein vector
    
    if (proVec[index] == 'STOP'){ # if codon is STOP, return protein vector
      output <- paste(unlist(proVec[1:index-1]), collapse = "")
      return(output)}
  }
  output <- paste(unlist(proVec), collapse = "")
  return(output) # Return protein vector
}

# Function that finds the most frequent AA from a string of all AA in the data
mostFrequentAA <- function(x) {
  tab <- table(strsplit(x, '')[[1]])
  names(tab)[tab == max(tab)]
}

# Function that finds the most frequent 3-mer AA from a list of AA strings
MerFunc <- function(AAlist, wantedOutput){
  
  # Create list of all AA 3-mer 
  merAA <- vector(mode="list", length = 0)
  counter <- 1
  
  # run for each cell (row) in the data
  for (i in 1:length(AAlist)){if (is.character(AAlist[[i]])){
    
    currStr <- AAlist[[i]] # Takes current AA string
    
    # Add 3-mer strings to the list "merAA"
    for (j in 1:(nchar(currStr)-2)){
      merAA[counter] <- substring(currStr, j, j+2)
      counter <- counter + 1
    }}}
  
  # Find the most frequent 3-mer AA 
  mostFrequentMerAAMat <- count(unlist(merAA)) # Create frequency matrix for 3-mer
  mostFrequentMerAAMat <- mostFrequentMerAAMat[order(mostFrequentMerAAMat$freq, decreasing = T),] # Order the frequency with descending order 
  mostFrequentMerAA <- mostFrequentMerAAMat$x[1] # Take the most frequent 3-mer 
  
  # Return the wanted output
  if (wantedOutput == 'Sin') {return(mostFrequentMerAA)} # Output - the most frequent 3-mer
  else if (wantedOutput == 'Mat') {return(mostFrequentMerAAMat)} # Output - the frequency matrix
  
}


### Question 1 ###

## Section 1 ##

# Load data from the file
Data <- read.table(file = 'HW2_ds1.tsv', sep = '\t', header = TRUE, fill = TRUE)

## Section 2 ## 

# Add a new column to the data with the length of the junction
lenVec <- lapply(Data$junction, function(x) {nchar(x)}) # Calculate the length of each junction
Data$Junction_length <- lenVec # Add the column to the data

## Section 3 ## 

# Filter out rows that the junction length is not modulo 3
wantedIndex <- sapply(Data$Junction_length, function(x){(x%%3) == 0})
Data <- Data[wantedIndex,]

## Section 4 ##

# Number of observations left after applying the filter
obsNumLeft <- nrow(Data)
sprintf('The Data has %d observations left', obsNumLeft)

## Section 5 ## 

isotypes <- unique(Data$c_call) # find all isotypes
aveLen <- vector(mode = "list", length = length(isotypes)) # Define vector for averages

# Loop to calculate average length for each isotype
for (i in 1:length(isotypes)){
  currIso <- grep(isotypes[i],Data$c_call)# find the indices of current isotype
  aveLen[i] <- mean(unlist(Data$Junction_length[currIso]))# Calculate mean length for current isotype
}

aveLen <- unlist(aveLen) # unlist the variable for the dataframe
averageIsoLen <- data.frame(isotypes, aveLen) # Create dataframe for average length of isotypes

## Section 6 ##

# Create dictionary for translation of DNA into amino acids
AA <- c('F','F','L','L','L','L','L','L','I','I','I','M','V','V','V','V','S','S','S','S','P','P','P','P','T','T','T','T','A','A','A','A','Y','Y','STOP','STOP','H','H','Q','Q','N','N','K','K','D','D','E','E','C','C','STOP','W','R','R','R','R','S','S','R','R','G','G','G','G')
names(AA) <- c('TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG')

# Translate DNA to AA
junction_aa <- lapply(Data$junction, function(x){if(str_count(x, "[ATGC]") == nchar(x)){TranslateDNA(x, AA)}})

# Add AA strings to new column in Data frame
Data$junction_aa <- junction_aa

## Section 7 ##

# Create a long string that contains all AA sequences
allAA <- paste(Data$junction_aa, collapse = '')

# Print most frequent AA that was found
print(mostFrequentAA(allAA))

## Section 8 ##

# Find the most frequent 3-mer for the current data
mostFreqMerAA <- MerFunc(Data$junction_aa, 'Sin')



### Question 2 ###

# Load data 
DataQue2 <- read.table(file = 'HW2_ds.tsv', sep = '\t', header = TRUE, fill = TRUE)

## Section 1 ## 

# Calculate length of string and add as a new column
lenVec2 <- lapply(DataQue2$junction, function(x) {nchar(x)})
DataQue2$Junction_length <- lenVec2

# Filter out rows that the junction length is not modulo 3
wantedIndex <- sapply(DataQue2$Junction_length, function(x){(x%%3) == 0})
DataQue2 <- DataQue2[wantedIndex,]

# Translate DNA to AA
junction_aa2 <- lapply(DataQue2$junction, function(x){if(str_count(x, "[ATGC]") == nchar(x)){TranslateDNA(x, AA)}})

# Add AA strings to new column in Data frame
DataQue2$junction_aa <- junction_aa2

# Find all the 3-mer combinations for the col names
matColNames <- MerFunc(DataQue2$junction_aa, 'Mat')
matColNames <- matColNames$x

# Find names of samples
samplesName <- unique(DataQue2$subject)

# Initialize frequency matrix
freqMatrix <- data.frame(matrix(0,ncol = length(matColNames), nrow = 9))
colnames(freqMatrix) <- matColNames

# Loop to find 3-mer frequency for every sample
for (i in 1:length(samplesName)){
  
  currSampleIndex <- grep(samplesName[i],DataQue2$subject) # Find indices of current sample
  currFreqMat <- MerFunc(DataQue2$junction_aa[currSampleIndex], 'Mat') # Run function and get frequency matrix for current sample
  currRow <- i # Define current row for general frequency matrix ('freqMatrix')
  
  # Run loop for all the 3-mer found for current sample
  for (j in 1:nrow(currFreqMat)){
    currCol <- grep(currFreqMat$x[j], colnames(freqMatrix)) # Find the matching column for current 3-mer
    freqMatrix[currRow, currCol] <- freqMatrix[currRow, currCol] + currFreqMat$freq[j] # Add the current 3-mer's frequency to matrix
    
  }}

## Section 2 ## 

# Create matrix of euclidean distance
distEuc <- dist(freqMatrix, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distEuc <- as.matrix(distEuc)

## Section 3 ## 

# Open png file for euclidean heatmap
png('HMdistEuc.png')

# Represent the euclidean distance matrix as a heatmap
HMdistEuc <- heatmap(distEuc)
dev.off()

## Section 4 ##  

# Create matrix of manhattan distance
distMan <- dist(freqMatrix, method = "manhattan", diag = FALSE, upper = FALSE, p = 2)
distMan <- as.matrix(distMan)

## Section 5 ## 

# Open png file for manhattan heatmap
png('HMdistMan.png') 

# Represent the manhattan distance matrix as a heatmap
HMdistMan <- heatmap(distMan)
dev.off()

## Section 6 ##

# Create frquency matrix for VJ combination
VJvec <- paste(Data$v_call, Data$j_call, sep=",") # Paste the two wanted columns
VJcounts <- count(VJvec) # Create frequency matrix
VJcounts <- VJcounts[order(VJcounts$freq, decreasing = T),] # Order the frequency with descending order

# Create matrices of euclidean and manhattan distance for the VJ combinations
distEucVJ <- dist(VJcounts$freq, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distManVJ <- dist(VJcounts$freq, method = "manhattan", diag = FALSE, upper = FALSE, p = 2)
distEucVJ <- as.matrix(distEucVJ)
distManVJ <- as.matrix(distManVJ)

# Represent the matrices of VJ combinations as a heatmaps
png('HMdistEucVJ.png') # Open png file for euclidean heatmap
HMdistEucVJ <- heatmap(distEucVJ) # Create heatmap
dev.off()

png('HMdistManVJ.png') # Open png file for manhattan heatmap
HMdistManVJ <- heatmap(distManVJ) # Create heatmap
dev.off()

## Section 7 & 8 ##

# Get indices for each sample
indices <- sapply(samplesName, function(x){grep(x, DataQue2$subject)})

# Boxplot the length variation for each sample with and without the outliers
for (x in indices){
  currData <- unlist(DataQue2$Junction_length[unlist(x)]) # Get data of current sample
  png(paste(DataQue2$subject[x[[1]][1]],'with.png'))
  # Create boxplot of data with outliers and find indices of outliers 
  indicesToRemove <- which(currData %in% boxplot(currData, main=DataQue2$subject[x[[1]][1]], ylab="Junction Length (bp)")$out)
  dev.off()
  currData <- currData[-indicesToRemove] # Remove outliers from Data
  png(paste(DataQue2$subject[x[[1]][1]],'without.png'))
  # Create boxplot of data without outliers
  boxplot(currData, main=DataQue2$subject[x[[1]][1]], ylab="Junction Length (bp)")
  dev.off()
  }




### Question 3 ###

## Section 1 ##

