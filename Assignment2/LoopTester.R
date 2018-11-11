# compute sums of windows of size 3 over a small sequence of integers, with windows jumping by 3

x = 1:30
window <- 3
myfunc<- function(y) print(x[y:(y+window-1)] )
z = seq(1,28,by=3)
zz <- sapply(z,myfunc)
zz

