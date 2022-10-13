#!usr/bin/Rscript

.libPaths("/data/kronos/kronos/acarrasco/R_libs")
library(data.table)
library(dplyr)
library(tidyverse)

ChunkFileCreator <- function(nSamples, ChunkSize) {
  numberSamples = nSamples
  totalChunks = ceiling(numberSamples/ChunkSize)
  Skipnumber = 0
  fileNumber = 1
  ChunkSize = ChunkSize
  columnnames = colnames(fread("SNPsDF_recodeA.raw", nrow = 1, drop=c(1,3:6)))
  for (chunk in seq(1, totalChunks)) {
    
    
    if (chunk == 1) {
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = ChunkSize - 1,
                        skip = Skipnumber, drop=c(1,3:6))
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    else if ( Skipnumber + ChunkSize > numberSamples) {
      
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = numberSamples - Skipnumber,
                        skip = Skipnumber, drop=c(1,3:6), col.names = columnnames)
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    else {
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = ChunkSize,
                        skip = Skipnumber, drop=c(1,3:6), col.names = columnnames)
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    Skipnumber = Skipnumber + ChunkSize
    fileNumber = fileNumber + 1
  }
  
}


# Calling the function to split the big genetic file and read it
# Make sure the number of samples is equal to the number of samples + 1. This is because on the .raw file, one of the rows is for the colnames)

ChunkFileCreator(1699, 150)