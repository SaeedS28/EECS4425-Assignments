# Question 2

library(seqinr)

# Part a
dnaRaw1 <- read.fasta(file="sequence1.fasta") # Import fasta file
seq1 <- getSequence(dnaRaw1[[1]], as.string = FALSE) # Extracts the nucleotides from fasta file
count(seq1,1, freq = TRUE) # outputs the percentage values of the four bases that make up the sequence
count(seq1,1, freq = FALSE) # outputs the freqency of bases