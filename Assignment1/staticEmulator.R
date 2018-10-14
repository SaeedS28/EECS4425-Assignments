f <- function(x) {
  y <- attr(f, "sum")
  if (is.null(y)) {
    y <- 0
  }
  print(sprintf("value of y = %i",y))
  y <- y+1
  attr(f, "sum") <<- y
  return(y)
}

k=0

while (k<20) {
  f(k)
  k <- k+1
}