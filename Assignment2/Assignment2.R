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