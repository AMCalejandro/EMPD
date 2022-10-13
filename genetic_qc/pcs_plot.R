#!/bin/Rscript

#--- Setting up the libpath, loading necesary libraries, and setting up some variable names ----####
args<-commandArgs(trailingOnly=TRUE)
library.path <- .libPaths()
library("tidyverse", lib.loc = library.path)

args = commandArgs(trailingOnly = TRUE)
cohort = args[1]


#---Load principal components generated from GCTA as well as the HapMap populations file---####
#This is the PCA from out data  merged with HapMap
#Using linkage-pruned SNPs, MAF > 5%, excluding palindromic SNPs and flipping missnps


#Read in eigenvectors
geneticPCA_hapmap.eigenvec <- as_tibble(read.table(args[2], sep = ""))

#Change column names
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>%
  dplyr::rename(FID = V1,
                IID = V2,
                PC1 = V3,
                PC2 = V4,
                PC3 = V5,
                PC4 = V6,
                PC5 = V7,
                PC6 = V8,
                PC7 = V9,
                PC8 = V10,
                PC9 = V11,
                PC10 = V12)

#Read in hapmap populations file
HapMap_pops <- as_tibble(read.table("/data/kronos/kronos/acarrasco/FILES/CEU.CHB.JPT.YRI_hapmap.txt",
                                    header = TRUE))

#---Merge PCA data with HapMap population data---####
geneticPCA_hapmap.eigenvec$FID <- as.factor(geneticPCA_hapmap.eigenvec$FID)
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>%
  left_join(HapMap_pops)

# If the population is missing, these should be the samples we are studying
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>%
  dplyr::mutate(group = ifelse(!is.na(population), population, cohort))

geneticPCA.eigenvec.CEU <- geneticPCA_hapmap.eigenvec %>%
  filter(group == "CEU")

#---Plot first two Principal Components with HapMap samples---####
ggplot(data = geneticPCA_hapmap.eigenvec, aes(x = PC1, y = PC2, colour = group, size = group)) +
  geom_point() +
  geom_vline(xintercept = mean(geneticPCA.eigenvec.CEU$PC1) + 7*(sd(geneticPCA.eigenvec.CEU$PC1))) +
  #geom_hline(yintercept = mean(geneticPCA.eigenvec.CEU$PC2) + 6*(sd(geneticPCA.eigenvec.CEU$PC2))) +
  #geom_vline(xintercept = mean(geneticPCA.eigenvec.CEU$PC1) - 6*(sd(geneticPCA.eigenvec.CEU$PC1))) +
  geom_hline(yintercept = mean(geneticPCA.eigenvec.CEU$PC2) - 7*(sd(geneticPCA.eigenvec.CEU$PC2))) +
  scale_color_manual(values = rainbow(5)) +
  scale_size_manual(values=c(5,3,3,1,3)) +
  theme_bw() +
  ggtitle(paste0("Hapmap reference panel populations and ", cohort, " cohort")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(paste0("PCA_HapMaPops_",cohort,".png"))