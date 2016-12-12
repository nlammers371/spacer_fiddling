rm(list = ls())
setwd(getwd())
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

wtEve2 <- as.character(unlist(read.csv(paste0(ReadPath,'eve2depace.txt'), header = FALSE)))
spc1Eve2 <- as.character(unlist(read.csv(paste0(ReadPath,'eve2spc1depace.txt'), header = FALSE)))

sp1Breaks <- as.vector(unlist(gregexpr('[a-z][A-Z]|[A-Z][a-z]',spc1Eve2)))

pInd <- 1

for (i in 1:length(sp1Breaks))