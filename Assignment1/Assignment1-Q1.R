# Code for Part 1 of Question 1
i <- 0
#subSeq <- "ACTG"
dnaSeq1 <- character()

while (i < 100) { # The loop nuns 400 times because 4 bases x 100 = 400 base strand
  dnaSeq1 <- append(dnaSeq1,"A")
  dnaSeq1 <- append(dnaSeq1,"C")
  dnaSeq1 <- append(dnaSeq1,"T")
  dnaSeq1 <- append(dnaSeq1,"G")
  i = i+1
}
print(dnaSeq1)
length(dnaSeq1)


# Code for Part 2 of Question 1
dnaSeq2 <- ""
randomValGen <- sample(1:4, 400, replace=T) #Generates 400 random integers between 1-4 (representing each base)
j <- 1
#print(randomValGen)

while(j<401){
  # Maps the generated random value to a base  
  if(randomValGen[j]==1){
    dnaSeq2 = paste(dnaSeq2,'A',sep="")
  }
  else if(randomValGen[j]==2){
    dnaSeq2 = paste(dnaSeq2,"C",sep="")
  }
  else if(randomValGen[j]==3){
    dnaSeq2 = paste(dnaSeq2,"T",sep="")
  }
  else{
    dnaSeq2 = paste(dnaSeq2,"G",sep="")
  }
  j = j+1
}
print(dnaSeq2)

#print(nchar(dnaSeq2))


# Code for Part 3 of Question 1
dnaSeq3 <- ""
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
    dnaSeq3 = paste(dnaSeq3,'A',sep="")
  }
  else if(randVal==2){
    dnaSeq3 = paste(dnaSeq3,"C",sep="")
  }
  else if(randVal==3){
    dnaSeq3 = paste(dnaSeq3,"T",sep="")
  }
  else{
    dnaSeq3 = paste(dnaSeq3,"G",sep="")
  }
  k = k+1
}
print(dnaSeq3)

#print(nchar(dnaSeq3))