library('seqinr')

# Question 1
lengthExon <- integer()
i <- 1

# Puts the first value in the vector to get things going
lengthExon[i] <- sample(500:10000, replace=T)

# Now checks to see if the sum is still less than 100 000
while(sum(lengthExon) <= 100000){
  i <- i+1
  lengthExon[i] <- sample(500:10000, replace=T)
}
i <- 1
# If the sum is greater than 100 000, it subtracts the difference from the element with the highest value
# to get things back to 100 000.
maxIndex <- which.max(lengthExon)

difference <- sum(lengthExon) - 100000 # Sum will always be greater than or equal to 100000
lengthExon[maxIndex] <- lengthExon[maxIndex]-difference
print(sprintf("Sum of the vector: %s", sum(lengthExon)))

# Now randomize the process of finding the starting point of the exons

startPoints <-integer()
startPoints[i] <- sample(500:10000, replace=T) # Seed

while(i <= length(lengthExon)){
  i <- i+1
  startPoints[i] <- startPoints[i-1]+lengthExon[i-1]+sample(12000:18000, replace=T) # Randomizes the start points based on previous results
}

print(sprintf("Length of startVector: %s  Length of exons: %s", length(startPoints), length(startPoints)))

print(lengthExon)
print(startPoints)

dnaSequence <- rep('D',1e6)
#print(length(dnaSequence))

ind <- 1

while(ind <= 1e6){
  if(ind %in% startPoints){
    grab <- match(ind,startPoints)
    looper <- startPoints[grab]+lengthExon[grab]
    for (counter in grab:looper ) {
      if(counter %% 3 == 0){
        dnaSequence[ind] <- append(dnaSequence, sample(c('A','C','T','G'), 1, replace=TRUE, prob=c(0.5,0.25,0.15,0.1)))
      }
      else{
        dnaSequence[ind] <- sample(c('A','C','T','G'), 1, replace=TRUE, prob=c(0.25,0.25,0.25,0.25))
      }
      print(ind)
      ind <- ind + 1
    }
  }
  else{
    dnaSequence[ind] <- sample(c('A','C','T','G'), 1, replace=TRUE, prob=c(0.25,0.25,0.25,0.25))
    print(ind)
    ind <- ind + 1
  }
  
}
