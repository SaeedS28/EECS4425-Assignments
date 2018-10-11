#install.packages("seqinr")
# Question 2

library(seqinr)

# Part a
dnaRaw1 <- read.fasta(file="sequence1.fasta") # Import fasta file
seq1 <- getSequence(dnaRaw1[[1]], as.string = FALSE) # Extracts the nucleotides from fasta file
count(seq1,1, freq = TRUE) # outputs the percentage values of the four bases that make up the sequence
count(seq1,1, freq = FALSE) # outputs the freqency of bases


# Part b
count(seq1,2) # outputs the dimers for the sequence


# Part c
count(seq1,3) # outputs the dimers for the sequence



# Question 3

extractAndIndicate <- function(seqInFunc, startIndex, endIndex){
  # checking for illegal arguments
  if(is.character(seqInFunc)){
    print("Data passed in is a char vector")
  }
  else{
    stop('The data type passed is not a character')
  }
  if((startIndex < 1) || (endIndex > length(seqInFunc)) || (startIndex > length(seqInFunc)) || (endIndex < startIndex)){
    stop('Index out of bounds')
  }
  
  # fetches the subsequence from the bigger sequence
  subsequence <- seqInFunc[startIndex:endIndex]
  print(subsequence)
  
  # generate indicator sequence
  counter <- 1
  
  indicatorSequence <- numeric()
  
  while (counter <= length(subsequence)) {
    if(subsequence[counter] == "g"){
      indicatorSequence <- append(indicatorSequence,1)
    }
    else{
      indicatorSequence <- append(indicatorSequence,0)
    }
    counter <- counter + 1
  }
  print(indicatorSequence)
}

extractAndIndicate(seqInFunc = seq1, startIndex = 1, endIndex = 10)