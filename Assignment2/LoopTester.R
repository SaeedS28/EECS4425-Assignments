x<- sample(c('A','C','T','G'), size=600,replace=TRUE) 
rseq<- seq(3,600,3) 
x[rseq]<- sample(c('A','C','T','G'), size=200,prob=c(0.5,0.25,0.15,0.1),replace=TRUE)

print(length(x))

y <- seq(400,600,2)
xseq <- y[y%%4==0]
print(xseq)

