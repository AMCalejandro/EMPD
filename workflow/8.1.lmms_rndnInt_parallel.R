#!usr/bin/Rscript

.libPaths("/data/kronos/kronos/acarrasco/R_libs")
library(plyr)           # To split up long data and apply a function to each group
library(tidyverse)      # For tidy manipulation of data
library(here)           # For file path construction
library(DT)             # To display data tables
library(lme4)           # For mixed models generation
library(MuMIn)          # To get marginal and conditional R2 values
library(lmerTest)       # To get pvalues based on Satterthwaite approach
library(doParallel)     # To setup parallel backend
library(foreach)        # To run a parallelizable for loop
library(data.table)

# Loading data
OPDC_PCs <- read.table("../Rsq0.4_3PCs/PCA.eigenvec", quote="\"", comment.char="")
OPDC_PCs <- OPDC_PCs[, 2:12]
colnames(OPDC_PCs) <- c("ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# Clinical_data + Outcome variables
OPDC_clinical <- readRDS("../Rsq0.4_3PCs/OPDC_final.rds")

# Coming up with a function so that I can parallelize this process
GWAS_earlyPD <- function(clinical, PCs, chunk_path) {
  #Loading libraries
  library(data.table)
  library(tidyverse)
  library(survival)
  
  #geneticChunk <- read.csv(paste0("GeneticChunksFiles/genetic_",i,".txt"), sep = "",
  #                         header = TRUE)
  geneticChunk <- read.csv(chunk_path, sep = "",  header = TRUE)
  
  PCs <- PCs
  clinical <- clinical
  
  TABLE<- PCs %>%
    inner_join(clinical) %>%
    inner_join(geneticChunk, by = c("ID" = "IID")) %>%
    arrange(ID, visit) %>%
    dplyr::mutate(visit = plyr::mapvalues(visit, from = 2, to = 1.5))
  
  
  columns_df <-colnames(TABLE)
  column_number = NULL
  for (column in 1:length(columns_df))  {
    to_check = columns_df[column]
    if (startsWith(columns_df[column], "chr") == "TRUE") {
      column_number = column
      break
    }
  }
  
  coefficients_1 <-as.data.frame(matrix(ncol= 13))
  names(coefficients_1) <- c("SNP","Coeff_slope", "se_slope", "pVal_slope","tval_slope", 
                             "coeff_incpt", "se_intcpt", "pVal_intcpt", "BIC", "MargR2", 
                             "CondR2", "N", "OutcomeVar")
  
  index_1 = 0
  for (z in column_number:ncol(TABLE)) {
    snp <- TABLE[,c(z,1:(column_number - 1 ))]
    
    model_total <- lmerTest::lmer(UPDRSIIItotal_imputed ~ 1 +
                                    + (1 + visit | ID) +
                                    snp[,1] +
                                    visit*snp[,1] +
                                    visit +
                                    age_diag.std + gender +
                                    PC1 + PC2 + PC3 + PC4 + PC5,
                                  data = snp,
                                  REML = TRUE)
    
    model_limb <- lmerTest::lmer(UPDRSIIIlimb_imputed ~ 1 +
                                   + (1 + visit | ID) +
                                   snp[,1] +
                                   visit*snp[,1] +
                                   visit +
                                   age_diag.std + gender +
                                   PC1 + PC2 + PC3 + PC4 + PC5,
                                 data = snp,
                                 REML = TRUE)
    
    model_axial <- lmerTest::lmer(UPDRSIIIaxial_imputed ~ 1 +
                                    + (1 + visit | ID) +
                                    snp[,1] +
                                    visit*snp[,1] +
                                    visit +
                                    age_diag.std + gender +
                                    PC1 + PC2 + PC3 + PC4 + PC5,
                                  data = snp,
                                  REML = TRUE)
    
    index_1 =  index_1 + 1
    coefficients_1[index_1, 1] = colnames(TABLE)[z]
    coeffs <- coef(summary(model_total))
    coefficients_1[index_1, 2] = coeffs[11,1]
    coefficients_1[index_1, 3] = coeffs[11,2]
    coefficients_1[index_1, 4] = coeffs[11,5]
    coefficients_1[index_1, 5] = coeffs[11,4]
    coefficients_1[index_1, 6] = coeffs[2,1]
    coefficients_1[index_1, 7] = coeffs[2,2]
    coefficients_1[index_1, 8] = coeffs[2,5]
    coefficients_1[index_1, 9] = BIC(model_total)
    r2 = MuMIn::r.squaredGLMM(model_total)
    coefficients_1[index_1, 10] = r2[1]
    coefficients_1[index_1, 11] = r2[2]
    coefficients_1[index_1, 12] = model_total@devcomp$dims[[1]] # Number of data points inc>
    coefficients_1[index_1, 13] = "model_total"
    
    index_1 =  index_1 + 1
    coefficients_1[index_1, 1] = colnames(TABLE)[z]
    coeffs <- coef(summary(model_limb))
    coefficients_1[index_1, 2] = coeffs[11,1]
    coefficients_1[index_1, 3] = coeffs[11,2]
    coefficients_1[index_1, 4] = coeffs[11,5]
    coefficients_1[index_1, 5] = coeffs[11,4]
    coefficients_1[index_1, 6] = coeffs[2,1]
    coefficients_1[index_1, 7] = coeffs[2,2]
    coefficients_1[index_1, 8] = coeffs[2,5]
    coefficients_1[index_1, 9] = BIC(model_limb)
    r2 = MuMIn::r.squaredGLMM(model_limb)
    coefficients_1[index_1, 10] = r2[1]
    coefficients_1[index_1, 11] = r2[2]
    coefficients_1[index_1, 12] = model_limb@devcomp$dims[[1]] # Number of data points inc>
    coefficients_1[index_1, 13] = "model_limb"
    
    
    index_1 =  index_1 + 1
    coefficients_1[index_1, 1] = colnames(TABLE)[z]
    coeffs <- coef(summary(model_axial))
    coefficients_1[index_1, 2] = coeffs[11,1]
    coefficients_1[index_1, 3] = coeffs[11,2]
    coefficients_1[index_1, 4] = coeffs[11,5]
    coefficients_1[index_1, 5] = coeffs[11,4]
    coefficients_1[index_1, 6] = coeffs[2,1]
    coefficients_1[index_1, 7] = coeffs[2,2]
    coefficients_1[index_1, 8] = coeffs[2,5]
    coefficients_1[index_1, 9] = BIC(model_axial)
    r2 = MuMIn::r.squaredGLMM(model_axial)
    coefficients_1[index_1, 10] = r2[1]
    coefficients_1[index_1, 11] = r2[2]
    coefficients_1[index_1, 12] = model_axial@devcomp$dims[[1]] # Number of data points inc>
    coefficients_1[index_1, 13] = "model_axial"
  }
  coefficients_1
}


files = list.files(path = "../Rsq0.4_3PCs/GeneticChunksFiles", full.names=T)
total_chunks = length(files)

cores=detectCores()
cl <- makeCluster(20) #not to overload your computer
clusterEvalQ(cl, .libPaths("/data/kronos/kronos/acarrasco/R_libs"))
registerDoParallel(cl)

finalMatrix <- foreach(i=1:total_chunks, .combine=rbind) %dopar% {
  myfile = files[i]
  chunkMatrix = GWAS_earlyPD(clinical = OPDC_clinical, PCs = OPDC_PCs, chunk = myfile)
  chunkMatrix #Equivalent to finalMatrix = cbind(finalMatrix, chunkMatrix)
}

stopCluster(cl)
fwrite(finalMatrix, file = "OPDC.GWAS_earlyPD_RStime_SNPTimeinteraction.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)