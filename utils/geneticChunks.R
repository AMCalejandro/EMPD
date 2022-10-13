#!usr/bin/Rscript

# Loading libraries
.libPaths("/data/kronos/kronos/acarrasco/R_libs")
library(dplyr)
library(data.table)
library(tidyverse)

# Loading the sample chunks of the genetic file
print("loading sample chunks of the genetic data")
partone <- readRDS("part1.rds")
parttwo <- readRDS("part2.rds")
partthree <- readRDS("part3.rds")
partfour <- readRDS("part4.rds")
partfive <- readRDS("part5.rds")
partsix <- readRDS("part6.rds")
partseven <- readRDS("part7.rds")
parteight <- readRDS("part8.rds")
partnine <- readRDS("part9.rds")
partten <- readRDS("part10.rds")
parteleven <- readRDS("part11.rds")
parttwelve <- readRDS("part12.rds")

print("rbinding the sample chunks of the genetic data")
genetic_data <- rbind(partone, parttwo, partthree, partfour, partfive,
                      partsix, partseven, parteight, partnine, partten, parteleven, parttwelve)

print("dimensions of the genetic dataset")
dim(genetic_data)

# Generation of the genetic chunks conataining 50K SNPs each
print("starting the for loop to generate all the 50000 SNPs  chunks")
numberSNPs = dim(genetic_data)[2]
n=1
SNPchunk = 50000
for (i in seq(2,numberSNPs, by = SNPchunk)) {
  start_index = as.integer(i)
  end_index = as.integer(i + 49999)
  if (end_index > numberSNPs) {end_index = numberSNPs}
  
  chunkDF = genetic_data[ , start_index:end_index]
  chunkDF = cbind(IID = genetic_data$IID, chunkDF)
  write.table(chunkDF, file = paste0("GeneticChunksFiles/genetic_",n,".txt"), col.names = TRUE, row.names = FALSE, sep = " ", quote = FALSE)
  n = n+1
}