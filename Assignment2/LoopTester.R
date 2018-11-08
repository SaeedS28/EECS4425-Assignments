test <- rep(1,1e6)

variable <- 1

while (variable <= length(test)) {
  test[variable] <- variable
  print(test[variable])
  variable <- variable+1
}