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
startPoints[i] <- sample(500:10000, replace=T)

while(i <= length(lengthExon)){
  i <- i+1
  startPoints[i] <- startPoints[i-1]+lengthExon[i-1]+sample(12000:18000, replace=T)
}

print(sprintf("Length of startVector: %s  Length of exons: %s", length(startPoints), length(startPoints)))

print(lengthExon)
print(startPoints)