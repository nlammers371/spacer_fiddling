rm(list = ls())
setwd(paste0(getwd(),"./GitHub/spacer_fiddling/code/utilities/"))
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

wtEve2 <- as.character(unlist(read.csv(paste0(ReadPath,'eve2depace.txt'), header = FALSE)))
spc1Eve2 <- as.character(unlist(read.csv(paste0(ReadPath,'eve2spc1depace.txt'), header = FALSE)))
wtDNA <- DNAString(wtEve2)
spc1DNA <- DNAString(spc1Eve2)

sp1Breaks <- append(as.vector(unlist(gregexpr('[a-z][A-Z]|[A-Z][a-z]',spc1Eve2))),nchar(spc1Eve2))

pInd <- 1
cons <- vector()
spacer1 <- vector()
for (i in 1:(length(sp1Breaks))){
  Ind <- sp1Breaks[i]
  sp1Seq <- substr(spc1Eve2,pInd,Ind)
  conSeq <- substr(wtEve2,pInd,Ind)
  spacer1 <- append(spacer1,sp1Seq)
  cons <- append(cons,conSeq)
  pInd <- Ind + 1
}

write(spacer1[],paste0(AnalyzePath,'spacer1.txt'))


