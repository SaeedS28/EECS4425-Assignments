
#install.packages("seqinr")
library(seqinr)

# I have left the bases as lower case characters because they seem to work with the seqinr package

# Question 1


# Part 1
i <- 0

dnaSeq1 <- character()

while (i < 100) { # The loop nuns 400 times because 4 bases x 100 = 400 base strand
  dnaSeq1 <- append(dnaSeq1,"a")
  dnaSeq1 <- append(dnaSeq1,"c")
  dnaSeq1 <- append(dnaSeq1,"t")
  dnaSeq1 <- append(dnaSeq1,"g")
  i = i+1
}
print(dnaSeq1)
#length(dnaSeq1)


# Part 2
dnaSeq2 <- character()
randomValGen <- sample(1:4, 400, replace=T) #Generates 400 random integers between 1-4 (representing each base)
j <- 1
#print(randomValGen)

while(j<401){
  # Maps the generated random value to a base  
  if(randomValGen[j]==1){
    dnaSeq2 <- append(dnaSeq2,"a")
  }
  else if(randomValGen[j]==2){
    dnaSeq2 <- append(dnaSeq2,"c")
  }
  else if(randomValGen[j]==3){
    dnaSeq2 <- append(dnaSeq2,"t")
  }
  else{
    dnaSeq2 <- append(dnaSeq2,"g")
  }
  j = j+1
}
#count(dnaSeq2,1)
print(dnaSeq2)

#print(nchar(dnaSeq2))


# Part 3
dnaSeq3 <- character()
k <- 1
randVal=0;

while (k<601) {
  
  # this generates a random number based on the probabilities given in the question 
  # sheet only when the number is divisible by three.
  # goes to default cause of equal probabilities otherwise.
  if(k %% 3 ==0){
    #print(sprintf("%i is divisible by 3",k))
    randVal = sample(c(1,2,3,4), 1, replace=TRUE, prob=c(0.5,0.25,0.15,0.1))
    #print(randVal)
    #print(sprintf("%i is divisible by 3",k))
  }
  else{
    randVal = sample(c(1,2,3,4), 1, replace=TRUE)
  }
  
  if(randVal==1){
    dnaSeq3 <- append(dnaSeq3,"a")
  }
  else if(randVal==2){
    dnaSeq3 <- append(dnaSeq3,"c")
  }
  else if(randVal==3){
    dnaSeq3 <- append(dnaSeq3,"t")
  }
  else{
    dnaSeq3 <- append(dnaSeq3,"g")
  }
  k = k+1
}
#count(dnaSeq3,1)
print(dnaSeq3)

#print(nchar(dnaSeq3))

# -----------------------------------------------------------------------------------------------------------

# Question 2

# Part a
dnaRaw1 <- read.fasta(file="sequence1.fasta") # Import fasta file
seq1 <- getSequence(dnaRaw1[[1]], as.string = FALSE) # Extracts the nucleotides from fasta file
count(seq1,1, freq = TRUE) # outputs the percentage values of the four bases that make up the sequence
count(seq1,1, freq = FALSE) # outputs the freqency of bases


# Part b
count(seq1,2) # outputs the dimers for the sequence


# Part c
count(seq1,3) # outputs the dimers for the sequence

# -----------------------------------------------------------------------------------------------------------

# Question 3

# Part 1
extractAndIndicate <- function(seqInFunc, startIndex, endIndex){
  # "static" variable for graph naming convention  
  staticVarCounter <- attr(extractAndIndicate, "graphCount")
  if (is.null(staticVarCounter)) {
    staticVarCounter <- 1
  }
  
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
  counter <- 1
  print(indicatorSequence)
  
  # fast Fourier transform
  fftCoefficients <- fft(indicatorSequence)
  print(fftCoefficients)
  
  fname = sprintf("fourierAnalysisGraph-%d.jpg", staticVarCounter)
  jpeg(fname)
  plot(fftCoefficients)
  dev.off
  
  staticVarCounter <- staticVarCounter + 1
  attr(extractAndIndicate, "graphCount") <<- staticVarCounter
  
}

extractAndIndicate(seqInFunc = seq1, startIndex = 1, endIndex = 12)


# Part 2
dnaRaw2 <- read.fasta(file="sequence2.fasta") # Import fasta file for 2nd sequence
#print(dnaRaw2)
codingSequence1 <- getSequence(dnaRaw2[[1]], as.string = FALSE)

# According to the annotation file, the following range of values identify the area of the genome
# that codes for a protien with protein_id WP_010877517.1
extractAndIndicate(seqInFunc = codingSequence1, startIndex = 4290, endIndex = 4784)

# The following range of values identify the area of the genome that doesn't code for any protein
extractAndIndicate(seqInFunc = codingSequence1, startIndex = 1, endIndex = 500)

# The function is now used on the sequence generated for Question 1, part c
extractAndIndicate(seqInFunc = dnaSeq3, startIndex = 1, endIndex = length(dnaSeq3))
