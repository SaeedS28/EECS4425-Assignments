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
while(i<length(endfinalRefined)){
  if((endfinalRefined[i]-startfinalRefined[i])>=5000){
    differenceOverFiveThousandExons[length(differenceOverFiveThousandExons)+1] <- i
  }
  i <- i+1
}

exonSequence <- ecoliSeq[startfinalRefined[differenceOverFiveThousandExons[1]]:endfinalRefined[differenceOverFiveThousandExons[1]]]

while(j<length(endfinalRefined-1)){
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


#Arg for fft exons
window <- 351 #characters displayed in every read
exonShift<- function(y) {
  (Arg(fft((indicatorGExon[y:(y+window-1)]))))[118]
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

hist(readWindowsIntron)
