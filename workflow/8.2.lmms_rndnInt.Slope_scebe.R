#!/bin/Rscript
# Loading packages
.libPaths("/mnt/rreal/RDS/RDS/acarrasco/R_libs")
library(MASS)
library(nlme)
library(Matrix)
library(lme4)
library(memoise)
library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(tidyverse)
library(data.table)
# Loading data   - #Modify this so that I can take data from stdinn
# Loading data
PCs <- read.table("/mnt/rreal/RDS/RDS/acarrasco/ANALYSES_WORKSPACE/EARLY_PD/OPDC/GENETICCHUNKS_NAMESUPDATED_2/PCA.eigenvec")[2:7]
colnames(PCs) = c("ID","PC1","PC2","PC3","PC4","PC5")
clinical_data = readRDS("/mnt/rreal/RDS/RDS/acarrasco/ANALYSES_WORKSPACE/EARLY_PD/OPDC/OPDC_final.rds")

# PhenoData N x 3 df ( ID, response, time_var)
phenoData = data.frame(ID = clinical_data$ID,
                       OUTCOME_total = clinical_data$UPDRSIIItotal_imputed,
                       OUTCOME_limb = clinical_data$UPDRSIIIlimb_imputed,
                       OUTCOME_axial = clinical_data$UPDRSIIIaxial_imputed,
                       YEAR = clinical_data$visit,
                       GENDER = clinical_data$gender,
                       AGE = clinical_data$age_diag.std) %>%
  inner_join(PCs) %>%
  arrange(ID, YEAR) %>%
  dplyr::mutate(YEAR = plyr::mapvalues(YEAR, from = c(1,2), to = c(0,1.5)))


# We get the phenotype data frame and the basic model fitted for each data framestring_combs = combn(c("total","limb","axial"), m = 2, simplify=F)
string_combs = combn(c("total","limb","axial"), m = 2, simplify=F)
names(string_combs) = c("axial", "limb", "total")
list_dfs = map(string_combs, function(x) dplyr::select(phenoData, !matches(x)))

fit0 = function(df) {
  df = as.data.frame(df)
  lme4::lmer(df[,2] ~ 1 + (1 + YEAR | ID) +
               YEAR + GENDER + AGE +
               PC1 + PC2 + PC3 + PC4 + PC5,
             REML = T,
             data = df)
}
base_models = map(list_dfs, fit0)


# We load in MEM the SCEBE algorithm
source(here::here("fun.dyngwas-final.R"))


# Get all the chunks files
data_path = "/mnt/rreal/RDS/RDS/acarrasco/ANALYSES_WORKSPACE/EARLY_PD/OPDC/GENETICCHUNKS_NAMESUPDATED_2/GeneticChunksFiles"
chunks = list.files(data_path)
#chunks_path = paste(data_path, chunks, sep = "/")[1:2]
chunks_path = paste(data_path, chunks, sep = "/")


# Custom function that maps genetic chunks to SCEBE
scebe_chunks = function(path, fit0, base_models, list_dfs) {
  genoData = fread(path) %>% dplyr::select(-IID)
  my_list <- vector(mode = "list", length = 3)
  
  for (i in seq_along(base_models)) {
    outcome_var = paste("OUTCOME", names(base_models)[i], sep = "_")
    
    res = scebe_sim(phenoData = list_dfs[[i]],
                    genoData = genoData,
                    fit0 = base_models[[i]],
                    Time = "YEAR",
                    pheno = outcome_var,
                    method = "scebe") %>%
      as.data.frame() %>%
      rownames_to_column(var = "SNP") %>%
      as_tibble() %>%
      mutate(OUTCOME = names(base_models)[i])
    
    #my_list[[i]] <- rbind(my_list[[i]], as_tibble(res) )
    my_list[[i]] <- res
  }
  my_list
}


start_time = Sys.time()

res_merged = map_dfr(chunks_path, ~ scebe_chunks(path = ., fit0 = fit0, base_models = base_models, list_dfs = list_dfs))
end_time <- Sys.time()
cat("\n Time to run ",end_time - start_time)

fwrite(res_merged, "OPDC.SCEBE.modelRandomSlopesTime_SCEBE.0.8.txt", quote=F, sep = "\t", row.names=F, col.names = T)