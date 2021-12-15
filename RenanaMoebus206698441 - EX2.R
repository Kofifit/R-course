library(magrittr)
library(plyr)
library(stringr)
library(ggplot2)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Question 1 ##

## Section 1 ##

# Function to create DNA strand
DNAseq <- function(DNAlen) {
  
  nuc <- list('A','T','G','C') # Define 4 nucleotides
  seqVec <- vector(mode = "list", length = DNAlen) # Define empty list for strand
  
  for(i in 1:DNAlen){ # Loop to randomize nucleotides for strand
    seqVec[i] <- nuc[sample(1:4,1,replace = FALSE)]}
  
  return(seqVec) # return DNA strand
}

## Section 2 ## 

# Create 5 DNA strands into one vector 
sequences <- c(DNAseq(21),DNAseq(21),DNAseq(21),DNAseq(21),DNAseq(21))
sequences <- paste(sequences, collapse='')

## Section 3 ##

# Create dictionary for translation of DNA into amino acids
AA <- c('F','F','L','L','L','L','L','L','I','I','I','M','V','V','V','V','S','S','S','S','P','P','P','P','T','T','T','T','A','A','A','A','Y','Y','STOP','STOP','H','H','Q','Q','N','N','K','K','D','D','E','E','C','C','STOP','W','R','R','R','R','S','S','R','R','G','G','G','G')
names(AA) <- c('TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGA','TGG','CGT','CGC','CGA','CGG','AGT','AGC','AGA','AGG','GGT','GGC','GGA','GGG')

## Section 4 ##

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

## Section 5 ##

# Create protein vector from DNA sequence 
proteins <- TranslateDNA(sequences,AA)

## Section 6 ##

# Print both vector - DNA and protein
print(sequences)
print(proteins)



## Question 2 ##

## Section 1 ##

# Load data from the file
Data <- read.table(file = 'HW2_ds1.tsv', sep = '\t', header = TRUE, fill = TRUE)

## Section 2 ##

# Find number of observations in the data
obsNum <- nrow(Data)
sprintf('The Data has %d observations left', obsNumLeft)

## Section 3 ##

# Filter out the rows that have the different genes assignments
vCell <- Data[,1] # Get the first column from the data
vCell <- strsplit(vCell, "\\*|,") # Split the genes based on ','
vCell <- sapply(vCell,function(x){unique(x)}) # Remove duplicate strings
## Section 6 ## - only next line
Data$v_gene <- sapply(vCell, function(x){x[[1]][1]}) # Add v_gene column
## Section 3 ## - continue
vCell <- sapply(vCell,function(x){length(grep("IGH",x)) == 1}) # Find indices of rows with only one gene
Data <- Data[vCell,] # Clean the data

## Section 4 ##

# Most frequent VJ combination
VJvec <- paste(Data$v_call, Data$j_call, sep=",") # Paste the two wanted columns
VJcounts <- count(VJvec) # Create frequency matrix
VJcounts <- VJcounts[order(VJcounts$freq, decreasing = T),] # Order the frequency with descending order
VJ_most_frequent <- VJcounts$x[1:5] # Take only the first five most frequent combinations

## Section 5 ##

# Most frequent Genes and alleles (v_Call)
V_FAMILY <- vector(mode = "list", length = 0) # Define new vector for 'V_FAMILY'
V_GENE <- vector(mode = "list", length = 0) # Define new vector for 'V_GENE'
V_ALLELES <- vector(mode = "list", length = 0) # Define new vector for 'V_ALLELES'
counter = 1 # Define counter for loop

# Loop to find most frequent genes and alleles for each gene family
for (i in 1:7){
  
  # Find most frequent Genes
  currFam <- sprintf('IGHV%i',i) # Define current gene family
  currDataGene <- sapply(Data[,1],function(x){length(grep(currFam,x))==1}) # Find all genes from current family
  currDataGene <- Data[currDataGene,1] # Copy data of the current genes to a variable
  currDataGene <- sapply(currDataGene, function(x){strsplit(x, ",")}) # Split the string to genes based on ','
  currDataGene <- sapply(currDataGene,function(x){gsub("\\*.*","",x)}) # Remove the allele from the strings
  currDataGene <- sapply(currDataGene,function(x){unique(x)}) # Remove duplicate strings (Leaving only one gene in cell)
  currDataGene <- count(currDataGene) # Create frequency matrix for genes
  currDataGene <- currDataGene[order(currDataGene$freq, decreasing = T),] # Order the frequency with descending order
  
  # If loop to make sure there are at least three genes
  if (nrow(currDataGene) >= 3){run = 3
  } else {run = nrow(currDataGene)}
  
  
  V_GENE[counter:(counter+run-1)] <- currDataGene[(1:run),1] # Add most frequent genes to the vector
  V_FAMILY[counter:(counter+run-1)] <- currFam # Add the gene family to the vector
  
  # Loop to find the most frequent alleles for each gene
  for (j in 1:run){
    currDataAllele <- sapply(Data[,1],function(x){length(grep(V_GENE[[(counter+j-1)]],x))==1}) # Find indices of current gene
    currDataAllele <- Data[currDataAllele,1] # Copy data of current gene to new variable
    currDataAllele <- sapply(currDataAllele,function(x){sub(".*\\*","",x)}) # Remove the gene from the string (Leaving only the allele)
    currDataAllele <- count(currDataAllele) # Create frequency matrix for alleles
    currDataAllele <- currDataAllele[order(currDataAllele$freq, decreasing = T),] # Order the frequency with descending order 
    V_ALLELES[(counter+j-1)] <- paste(V_GENE[[(counter+j-1)]],'*',currDataAllele[[1,1]])}# Take the most frequent allele
  
  
  counter = counter + run # Update counter 
}

V_FAMILY <- unlist(V_FAMILY) # Unlist the variable
V_GENE <- unlist(V_GENE) # Unlist the variable
V_ALLELES <- unlist(V_ALLELES) # Unlist the variable
v_data <- data.frame(V_FAMILY, V_GENE, V_ALLELES) # Add variables to data.frame

## Section 7 ##

# Add a new column to the data with the length of the junction
lenVec <- lapply(Data$junction, function(x) {nchar(x)}) # Calculate the length of each junction
Data$Junction_length <- lenVec # Add the column to the data

## Section 8 ##

# Find all genes in data
ALLgenes <- unique(Data$v_gene)

# Get length variation for all genes
indices <- sapply(ALLgenes, function(x){grep(x, Data$v_gene)})

# Boxplot the length variation for each gene
for (x in indices){
  png(paste(Data$v_gene[x[[1]][1]],'.png'))
  boxplot(unlist(Data$Junction_length[unlist(x)]), main=Data$v_gene[x[[1]][1]], ylab="Junction Length (bp)")
  dev.off()}

# Check results of boxplot with quantile function
sapply(indices, function(x){print(quantile(unlist(Data$Junction_length[unlist(x)])))})



