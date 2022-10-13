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
HapMap_pops <- as_tibble(read.table("/mnt/rreal/RDS/REFERENCE_GENOME/HAPMAP/relationships_w_pops_121708.txt",
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

geneticPCA.eigenvec.MyData <- geneticPCA_hapmap.eigenvec %>%
  filter(group == cohort)


## Remove mySample individuals who are outliers on any PC
#Make results table, with one row for each individual, including only the individuals of the current study
PC.outlierResults <- as.data.frame(matrix(ncol = 11,  nrow = length(which(geneticPCA_hapmap.eigenvec$group == cohort))))
PC.outlierResults[, 1] <- geneticPCA_hapmap.eigenvec$IID[which(geneticPCA_hapmap.eigenvec$group == cohort)] #Put your IIDs in the first column


#For loop to calculate the SDs of each Principal Component (first 10 PCs only)
#This outputs into a results table
for (i in 3:12) {
  mean <- mean(geneticPCA.eigenvec.CEU[[i]])
  sd <- sd(geneticPCA.eigenvec.CEU[[i]])
  
  PC.outlierResults[, i-1] <- geneticPCA.eigenvec.MyData %>%
    mutate(outlier = ifelse(geneticPCA.eigenvec.MyData[[i]] > mean + 6*sd, "outlier",
                            ifelse(geneticPCA.eigenvec.MyData[[i]] < mean - 6*sd, "outlier", "keep"))) %>%
    dplyr::select(outlier)
}

#Rename column names
PC.outlierResults <- PC.outlierResults %>%
  dplyr::rename(ID = V1,
                PC1_result = V2,
                PC2_result = V3,
                PC3_result = V4,
                PC4_result = V5,
                PC5_result = V6,
                PC6_result = V7,
                PC7_result = V8,
                PC8_result = V9,
                PC9_result = V10,
                PC10_result = V11)

#Now merge the outlier results with the main dataset
geneticPCA.eigenvec.MyData <- geneticPCA.eigenvec.MyData %>%
  left_join(PC.outlierResults, by = c("IID" = "ID"))

#If any of the PC results are outliers, flag as outlier
geneticPCA.eigenvec.MyData <- geneticPCA.eigenvec.MyData %>%
  mutate(PCA_outlier = ifelse(PC1_result == "outlier" |
                                PC2_result == "outlier" |
                                PC3_result == "outlier" |
                                PC4_result == "outlier" |
                                PC5_result == "outlier" |
                                PC6_result == "outlier" |
                                PC7_result == "outlier" |
                                PC8_result == "outlier" |
                                PC9_result == "outlier" |
                                PC10_result == "outlier", "outlier", "valid"))

#Plot first 2 PCs by outlier status (this includes both the HapMap samples and your samples).
ggplot(data = geneticPCA.eigenvec.MyData, mapping = aes(x = PC1, y = PC2, color = PCA_outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  ggtitle(paste0("Outliers on the ", cohort, " population")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(paste0("PCA_CEU_", cohort,"outliers.png"))


#Count how many outliers
geneticPCA.eigenvec.MyData %>%
  group_by(PCA_outlier) %>%
  dplyr::summarise(count = n())

#---Write list of samples that are population outliers---####
#To remove in PLINK

PCA_outliers <- geneticPCA.eigenvec.MyData %>%
  filter(PCA_outlier == "outlier") %>%
  select(FID, IID)

#Write text file of FID and IID
write.table(PCA_outliers, "PCA_outliersAncestry_4SD_PSP_CEUonly.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)