
# compute sums of windows of size 5 over a small sequence of integers
 window <- 5
  x<- 1:100         
 myfunc<- function(y) sum(x[y:(y+window-1)] )
z<- sapply(x,myfunc)
# print the result
z[1:26]

#do the same with a loop
y<- rep(0,100)
for(i in 1:26) y[i]<- sum(x[i:(i+window-1)] )
#check that the two methods produce the same results by
#printing the elementwise differences

z[1:26]- y[1:26]
