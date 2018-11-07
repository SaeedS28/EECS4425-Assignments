placeholder <- integer()
i <- 1

# Puts the first value in the vector to get things going
placeholder[i] <- sample(500:10000, replace=T)

# Now checks to see if the sum is still less than 100 000
while(sum(placeholder) <= 100000){
  i <- i+1
  placeholder[i] <- sample(500:10000, replace=T)
}

print(sum(placeholder))

# If the sum is greater than 100 000, it subtracts the difference from the element with the highest value
# to get things back to 100 000.
maxIndex <- which.max(placeholder)
print(placeholder[maxIndex])

difference <- sum(placeholder) - 100000 # Sum will always be greater than or equal to 100000
#print(difference)
placeholder[maxIndex] <- placeholder[maxIndex]-difference
print(sum(placeholder))