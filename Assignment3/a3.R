library('seqinr')
library(stringr)

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
j <- 2

# Four values of exons greater than or equal to 5000
while(i<length(endfinalRefined)){
  if((endfinalRefined[i]-startfinalRefined[i])>=5000){
    differenceOverFiveThousandExons[length(differenceOverFiveThousandExons)+1] <- i
  }
  i <- i+1
}

while(j<length(endfinalRefined)){
  if((startfinalRefined[j]-endfinalRefined[j-1])>=5000){
    differenceOverFiveThousandIntrons[length(differenceOverFiveThousandIntrons)+1] <- j
  }
  j <- j+1
}

