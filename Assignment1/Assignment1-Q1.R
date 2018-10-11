# Code for Part 1 of Question 1
dnaSeq1 <- ""
i <- 0
subSeq <- "ACTG"

while (i < 100) { # The loop nuns 100 times because 4 bases x 100 = 400 base strand
  dnaSeq1 = paste(dnaSeq1,subSeq,sep="") #string concatenation
  i = i+1
}

print(dnaSeq1)
print(nchar(dnaSeq1))


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
  if(k %% 3 ==0){
    print(sprintf("%i is divisible by 3",k))
  }
  k = k+1
}