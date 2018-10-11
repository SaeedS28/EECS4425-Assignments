# I have left the bases as lower case characters because they seem to work with the seqinr package

# Code for Part 1 of Question 1
i <- 0
#subSeq <- "ACTG"
dnaSeq1 <- character()

while (i < 100) { # The loop nuns 400 times because 4 bases x 100 = 400 base strand
  dnaSeq1 <- append(dnaSeq1,"a")
  dnaSeq1 <- append(dnaSeq1,"c")
  dnaSeq1 <- append(dnaSeq1,"t")
  dnaSeq1 <- append(dnaSeq1,"g")
  i = i+1
}
#print(dnaSeq1)
#length(dnaSeq1)


# Code for Part 2 of Question 1
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
#print(dnaSeq2)

#print(nchar(dnaSeq2))


# Code for Part 3 of Question 1
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
#print(dnaSeq3)

#print(nchar(dnaSeq3))