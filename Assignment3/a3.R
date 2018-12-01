library('seqinr')
library(stringr)

# Question 1
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

differenceOverFiveThousandExons <- integer()
differenceOverFiveThousandIntrons <- integer()

i <- 1
j <- 1

# Four values of exons greater than or equal to 5000
while(i<=length(endfinalRefined)){
  if((endfinalRefined[i]-startfinalRefined[i])>=5000){
    differenceOverFiveThousandExons[length(differenceOverFiveThousandExons)+1] <- i
  }
  i <- i+1
}

exonSequence <- ecoliSeq[startfinalRefined[differenceOverFiveThousandExons[1]]:endfinalRefined[differenceOverFiveThousandExons[1]]]

while(j<=length(endfinalRefined-1)){
  if((startfinalRefined[j+1]-endfinalRefined[j])>=5000){
    differenceOverFiveThousandIntrons[length(differenceOverFiveThousandIntrons)+1] <- j
  }
  j <- j+1
}

intronSequence <- ecoliSeq[endfinalRefined[differenceOverFiveThousandIntrons[1]]:startfinalRefined[differenceOverFiveThousandIntrons[1]+1]]


indicatorGExon <- character()
dnaCopyExon <- character()

dnaCopyExon <- exonSequence[1:length(exonSequence)]
dnaCopyExon <- tolower(dnaCopyExon)
indicatorGExon <- dnaCopyExon

# Replaces g's with ones and the rest with zeroes
indicatorGExon[dnaCopyExon !='g'] <- 0
indicatorGExon[dnaCopyExon =='g'] <- 1
indicatorGExon <- as.numeric(indicatorGExon)


indicatorGIntron <- character()
dnaCopyIntron <- character()

dnaCopyIntron <- intronSequence[1:length(intronSequence)]
dnaCopyIntron <- tolower(dnaCopyIntron)
indicatorGIntron <- dnaCopyIntron

# Replaces g's with ones and the rest with zeroes
indicatorGIntron[dnaCopyIntron !='g'] <- 0
indicatorGIntron[dnaCopyIntron =='g'] <- 1
indicatorGIntron <- as.numeric(indicatorGIntron)


# Arg for fft exons - Used for questions 2.1 and 2.2
window <- 351 #characters displayed in every read
exonShift<- function(y) {
  (Arg(fft((indicatorGExon[y:(y+window-1)]))))[118] # 118th value produces better results than 117th
}
shiftValsExon <- seq(1,length(indicatorGExon),by=3) # shift 3 after every read for the length of the sequence

readWindowsExon <- sapply(shiftValsExon,exonShift)
readWindowsExon <- readWindowsExon[!is.na(readWindowsExon)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

intronShift<- function(y) {
  (Arg(fft((indicatorGIntron[y:(y+window-1)]))))[118]
}
shiftValsIntron <- seq(1,length(indicatorGIntron),by=3) # shift 3 after every read for the length of the sequence

readWindowsIntron <- sapply(shiftValsIntron,intronShift)
readWindowsIntron <- readWindowsIntron[!is.na(readWindowsIntron)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

plot(1:length(readWindowsIntron), readWindowsIntron, type = "l", main ="Introns Sequence", xlab = "Window", ylab = "Phase")
plot(1:length(readWindowsExon), readWindowsExon, type = "l", main ="Exon Sequence", xlab = "Window", ylab = "Phase")


#Question 2
hist(readWindowsIntron,breaks = 200) # diplays uniformity

histOriginal <- hist(readWindowsExon, breaks=200) # displays values that clump together


# Deletes the middle value for 1 deletion
indicatorGExon1Del <- indicatorGExon[-round(length(indicatorGExon)/2)]

exonShift1Del<- function(y) {
  (Arg(fft((indicatorGExon1Del[y:(y+window-1)]))))[118]
}
shiftValsExon1Del <- seq(1,length(indicatorGExon1Del),by=3) # shift 3 after every read for the length of the sequence
readWindowsExon1Del <- sapply(shiftValsExon1Del,exonShift1Del)
readWindowsExon1Del <- readWindowsExon1Del[!is.na(readWindowsExon1Del)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

hist1Del <- hist(readWindowsExon1Del,breaks = 200) # Histogram for 1 deletion shows 2 clusters


# Deletes the 900th element from the already-deleted sequence
indicatorGExon2Del <- indicatorGExon1Del[-900]

exonShift2Del<- function(y) {
  (Arg(fft((indicatorGExon2Del[y:(y+window-1)]))))[118]
}
shiftValsExon2Del <- seq(1,length(indicatorGExon2Del),by=3) # shift 3 after every read for the length of the sequence
readWindowsExon2Del <- sapply(shiftValsExon2Del,exonShift2Del)
readWindowsExon2Del <- readWindowsExon2Del[!is.na(readWindowsExon2Del)] # deletes all the NA values from the data. Easier than figuring out the exact shiftVal sequence

hist2Del=hist(readWindowsExon2Del, breaks =85)


# Question 2.4
# Finding the three peaks
cPhiOriginal = hist1Del$breaks[which.max(hist1Del$counts)] # phi is the value that corresponds to the highest peak on the unedited sequence
cPhiLeftOriginal = cPhiOriginal-(2*pi/3)
cPhiRightOriginal = cPhiOriginal+(2*pi/3)

# Classifying the original sequence based on these peaks
leftHalfOriginal <-(cPhiOriginal+cPhiLeftOriginal)/2
rightHalfOriginal <- (cPhiOriginal+cPhiRightOriginal)/2

counts <- 1
counts1D <- 1
counts2D <- 1

classifierOriginal <- integer()
classifier1Del <- integer()
classifier2Del <- integer()
rollingSumOriginal <- integer()
rollingSum1Del <- integer()
rollingSum2Del <- integer()

countSum <- 2
sums<-0

# Questions 2.5 and 2.6
# Maps these sequences to values based on their distances to the closest phi values
while (counts<=length(readWindowsExon)) {
  val <- readWindowsExon[counts]
  if(val <= leftHalfOriginal){
    classifierOriginal[counts] <- -1
  }
  else if (val >= leftHalfOriginal && val <= rightHalfOriginal){
    classifierOriginal[counts] <- 0
  }
  else{
    classifierOriginal[counts] <- 1
  }
  counts <- counts +1
}

# Calculates a running total and stores these values in a vector for question 2.6
rollingSumOriginal[1] <- classifierOriginal[1]
while (countSum<=length(classifierOriginal)) {
  sums <- classifierOriginal[countSum]
  sums <- sums + rollingSumOriginal[countSum-1]
  rollingSumOriginal[countSum] <- sums
  countSum <- countSum + 1
}

rollingOriginalSeq <- seq(1,length(rollingSumOriginal),1)
plot(rollingOriginalSeq,rollingSumOriginal,main = "Exon sequence without any dels")

countSum <- 2

while (counts1D<=length(readWindowsExon1Del)) {
  val <- readWindowsExon1Del[counts1D]
  if(val <= leftHalfOriginal){
    classifier1Del[counts1D] <- -1
  }
  else if (val >= leftHalfOriginal && val <= rightHalfOriginal){
    classifier1Del[counts1D] <- 0
  }
  else{
    classifier1Del[counts1D] <- 1
  }
  counts1D <- counts1D +1
}

rollingSum1Del[1] <- classifier1Del[1]
while (countSum<=length(classifier1Del)) {
  sums <- classifier1Del[countSum]
  sums <- sums + rollingSum1Del[countSum-1]
  rollingSum1Del[countSum] <- sums
  countSum <- countSum + 1
}

countSum <- 2

rolling1DelSeq <- seq(1,length(rollingSum1Del),1)
plot(rolling1DelSeq,rollingSum1Del,main = "f(j) for the Exon sequence one del")


while (counts2D<=length(readWindowsExon2Del)) {
  val <- readWindowsExon2Del[counts2D]
  if(val <= leftHalfOriginal){
    classifier2Del[counts2D] <- -1
  }
  else if (val >= leftHalfOriginal && val <= rightHalfOriginal){
    classifier2Del[counts2D] <- 0
  }
  else{
    classifier2Del[counts2D] <- 1
  }
  counts2D <- counts2D + 1
}

countSum <- 2

rollingSum2Del[1] <- classifier2Del[1]
while (countSum<=length(classifier2Del)) {
  sums <- classifier2Del[countSum]
  sums <- sums + rollingSum2Del[countSum-1]
  rollingSum2Del[countSum] <- sums
  countSum <- countSum + 1
}

rolling2DelSeq <- seq(1,length(rollingSum2Del),1)
plot(rolling2DelSeq,rollingSum2Del,main = "f(j) Exon sequence two del")

# This analysis results in very unreliable data.
# Hence, the prediction is also incorrect.