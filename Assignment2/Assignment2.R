library('seqinr')
library(stringr)
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
startPoints[i] <- sample(5000:10000, replace=T) # Seed

while(i < length(lengthExon)){
  i <- i+1
  startPoints[i] <- startPoints[i-1]+lengthExon[i-1]+sample(12000:18000, replace=T) # Randomizes the start points based on previous results
}

print(sprintf("Length of startVector: %s  Length of exons: %s", length(startPoints), length(startPoints)))

#print(lengthExon)
#print(startPoints)

dnaSequence <- sample(c('A','C','T','G'), size= 1e6, replace = TRUE, prob = c(0.25,0.25,0.25,0.25)) #Fills the background with data
#print(length(dnaSequence))

ind <- 1


print("Building Sequence. Please wait...")

# Loops until the size of the startPoint vector since that's how many exon partitions there are
# Way faster than looping through 1000 000 times!
while(ind<=length(startPoints)){
  starts <- startPoints[ind] #sets the start point to the ones in the startPoints vector
  len <- lengthExon[ind]
  replacer <- seq(starts,starts+len)
  #print(replacer)
  subseq <-replacer[replacer %% 3 == 0] # extracts the values divisible by three
  dnaSequence[subseq]<- sample(c('A','C','T','G'), size=length(subseq),prob=c(0.5,0.25,0.15,0.1),replace=TRUE)
  #print(dnaSequence[subseq])
  ind <- ind+1
}

print("Sequence Built!")

# Question 2
# Part 1
indicatorGSeq <- character()
dnaCopy <- character()

dnaCopy <- dnaSequence[1:length(dnaSequence)]
dnaCopy <- tolower(dnaCopy)
indicatorGSeq <- dnaCopy

# Replaces g's with ones and the rest with zeroes
indicatorGSeq[dnaCopy!='g'] <- 0
indicatorGSeq[dnaCopy=='g'] <- 1
indicatorGSeq <- as.numeric(indicatorGSeq)
#print(indicatorGSeq)

# Compute FFTs for the sequence with a window size of 351
window <- 351 #characters displayed in every read
myfunc<- function(y) {
  (Mod(fft((indicatorGSeq[y:(y+window-1)]))))[117]
}
shiftVals <- seq(1,length(dnaSequence),by=3) # shift 3 after every read for the length of the sequence
# shiftVals

readWindows <- sapply(shiftVals,myfunc)
readWindows <- readWindows[!is.na(readWindows)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

# Define the threshold to be the midpoint 
maxThreshold <- max(readWindows)
minThreshold <- min(readWindows)
threshold <- (maxThreshold + minThreshold)/2
thresholdVector <- rep(0,length(readWindows))
thresholdVector[readWindows>threshold] <- 1
plot(1:length(thresholdVector), thresholdVector, type = "l",xlab = "points", ylab = "magnitude")

#Part 2
# Import the annotation file
rawData <- readLines("sequenceAnnotation.txt")
refinedData <- str_match(rawData,"location=[A-Za-z0-9$&+,:;=?@#|'<>-^*()%!]*\\.\\.[A-Za-z0-9$&+,:;=?@#|'<>-^*()%!]*") # regex used to find the locations in the annotation file
refinedData <- refinedData[!is.na(refinedData)]

# extract start location
furtherRefinedData <- str_match(refinedData,"\\d+\\.\\.") #Fetches the actual start locations
startfinalRefined <- gsub("\\.\\.", "", furtherRefinedData)
startfinalRefined <- as.numeric(startfinalRefined)
startfinalRefined <- startfinalRefined[!is.na(startfinalRefined)]

# extract end location
furtherRefinedData2 <- str_match(refinedData,"\\.\\.>*\\d+") #Fetches the actual end locations
endfinalRefined <- gsub("\\.\\.>*", "", furtherRefinedData2)
endfinalRefined <- as.numeric(endfinalRefined)
endfinalRefined <- endfinalRefined[!is.na(endfinalRefined)]

ecoliSeq <- read.fasta(file="sequenceFasta.fasta")
ecoliSeq <- getSequence(ecoliSeq[[1]], as.string = FALSE)

# generate indicator sequence
indicatorGSeqEcoli <- character()
dnaCopy2 <- character()

dnaCopy2 <- ecoliSeq[1:length(ecoliSeq)]
dnaCopy2 <- tolower(dnaCopy2)
indicatorGSeqEcoli <- dnaCopy2

# Replaces g's with ones and the rest with zeroes
indicatorGSeqEcoli[dnaCopy2!='g'] <- 0
indicatorGSeqEcoli[dnaCopy2=='g'] <- 1
indicatorGSeqEcoli <- as.numeric(indicatorGSeqEcoli)
#print(indicatorGSeq)

shiftValsEcoli <- seq(1,length(indicatorGSeqEcoli),by=3)
readWindowsEcoli <- sapply(shiftValsEcoli,myfunc)
readWindowsEcoli <- readWindows[!is.na(readWindowsEcoli)]

thresholdEcoli <- (max(readWindowsEcoli)+min(readWindowsEcoli))/2
thresholdVectorEcoli <- rep(0,length(readWindowsEcoli))
thresholdVectorEcoli[readWindows>thresholdEcoli] <- 1
plot(1:length(thresholdVectorEcoli), thresholdVectorEcoli, type = "l",xlab = "points", ylab = "magnitude")

# Find out where the exons in the indicator fft start and end
starterInd <- function(y) indicatorGSeqEcoli[y-1]==0 && indicatorGSeqEcoli[y] == 1  # Starting exon condition
starter <- sapply(1:length(indicatorGSeqEcoli), starterInd)

endInd <- function(y) indicatorGSeqEcoli[y]==1 && indicatorGSeqEcoli[y+1] == 0  # Ending of exon condition
ender <- sapply(1:length(indicatorGSeqEcoli), endInd)

starter <- as.integer(starter)
ender <- as.integer(ender)

# Calculates the true positives
truePos <- function(y) any(starter[(startfinalRefined[y]/3):(endfinalRefined[y]/3)] == 1 | ender[(startfinalRefined[y]/3): (endfinalRefined[y]/3)] == 1)
truePositive <- sapply(1:length(startfinalRefined), truePos)
truePositive <-sum(as.integer(truePositive))

#Question 2.3
indicatorGSeqEcoliGC <- character()
# dnaCopy2 <- character()
# 
# dnaCopy2 <- ecoliSeq[1:length(ecoliSeq)]
# dnaCopy2 <- tolower(dnaCopy2)
indicatorGSeqEcoliGC <- dnaCopy2

# Replaces g's with ones and the rest with zeroes
indicatorGSeqEcoliGC[dnaCopy2!='g' || dnaCopy2!='c'] <- 0
indicatorGSeqEcoliGC[dnaCopy2=='g' || dnaCopy2=='c'] <- 1
indicatorGSeqEcoli <- as.numeric(indicatorGSeqEcoli)

myfunc<- function(y) {
  (Mod(fft((indicatorGSeqGC[y:(y+window-1)]))))[117]
}
shiftVals <- seq(1,length(dnaSequence),by=3) # shift 3 after every read for the length of the sequence
# shiftVals

readWindows <- sapply(shiftVals,myfunc)
readWindows <- readWindows[!is.na(readWindows)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

shiftValsEcoli <- seq(1,length(indicatorGSeqEcoliGC),by=3)
readWindowsEcoli <- sapply(shiftValsEcoli,myfunc)
readWindowsEcoli <- readWindows[!is.na(readWindowsEcoli)]

#thresholdEcoli <- (max(readWindowsEcoli)+min(readWindowsEcoli))/2
thresholdVectorEcoli <- rep(0,length(readWindowsEcoli))
thresholdVectorEcoli[readWindows>thresholdEcoli] <- 1
plot(1:length(thresholdVectorEcoli), thresholdVectorEcoli, type = "l",xlab = "points", ylab = "magnitude")

# Find out where the exons in the indicator fft start and end
starterInd <- function(y) indicatorGSeqEcoliGC[y-1]==0 && indicatorGSeqEcoliGC[y] == 1  # Starting exon condition
starter <- sapply(1:length(indicatorGSeqEcoliGC), starterInd)

endInd <- function(y) indicatorGSeqEcoliGC[y]==1 && indicatorGSeqEcoliGC[y+1] == 0  # Ending of exon condition
ender <- sapply(1:length(indicatorGSeqEcoli), endInd)

starter <- as.integer(starter)
ender <- as.integer(ender)

# Calculates the true positives
truePos <- function(y) any(starter[(startfinalRefined[y]/3):(endfinalRefined[y]/3)] == 1 | ender[(startfinalRefined[y]/3): (endfinalRefined[y]/3)] == 1)
truePositive <- sapply(1:length(startfinalRefined), truePos)
truePositive <-sum(as.integer(truePositive))