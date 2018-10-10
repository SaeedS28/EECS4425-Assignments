# Code for Part 1 of Question 1
dnaSeq <- ""
i <- 1
subSeq <- "ACTG"

while (i < 101) { # The loop nuns 100 times because 4 bases x 100 = 400 base strand
  dnaSeq = paste(dnaSeq,subSeq,sep="")
  i = i+1
}

print(dnaSeq)
print(nchar(dnaSeq))