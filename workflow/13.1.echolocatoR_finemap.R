#!/bin/Rscript

## DESCRIPTION ##
# This is the workflow to automate the finemapping from echocolocatoR

# This script takes a txt file with the following informations
# RS_PATH
# SS_path
## Example usage
#Rscript FINEMAPPING.R [path_to_metadata_file.txt]



# Activate echoR conda env
#system("conda activate echoR")

# Load libraries
library(data.table) # Efficient multiprocessor data reading/writing/processing
library(echolocatoR) # To fo the finemapping
library(tidyverse) # For data wrangling / munging
library(tools) # To perform format changes in the SS file
library(rlang) # Data masking and more
library(vautils) # Use util to find nearest gene from an input SNP
library(here) # Efficient path management within the R project
# Read file with metadata to perform finemapping
library(colochelpR, lib.loc = "/mnt/rreal/RDS/acarrasco/R_libs") # Necessary to do som QC
library(SNPlocs.Hsapiens.dbSNP144.GRCh37, lib.loc = "/mnt/rreal/RDS/acarrasco/R_libs") # SNP reference panel to do some QC

# We load some functions that I created to automate running fienmapping as much as possible
# Functions can be found here -> https://github.com/AMCalejandro/echor_wrapperPipes
source("/mnt/rreal/RDS/acarrasco/TOOLS/echor_wrapperPipes/R/finemap_plotting.R")  



args <- commandArgs(trailingOnly = TRUE)
metadata = readLines(args[1])
# Making sure I remove any problematic white spaces
metadata = metadata[which(metadata!="")]
metadata = trimws(gsub("\\s+", "", metadata))

# Storing metadata information
# This script should take the data harmonised. Perform harmonization outside this script
fullSS_path = metadata[1]
fullRS_path = metadata[2]
newSS_name = metadata[3]


source("/mnt/rreal/RDS/acarrasco/TOOLS/echor_wrapperPipes/R/utils.R")

# Create RS if not exist
make_results_dir(fullRS_path)

# Read GWAS and do some wrangling
data = fread(fullSS_path)


# TODO
# Use function from postGWAS tool to get to know the SNP format and make decisions
# So, use source to load postGWAS QC tool. Then, harmonise input data


if ("MarkerName" %in% colnames(data)) {
  data = data %>% dplyr::rename(SNP = MarkerName, Pval = `P-value`) %>%
    select(SNP, everything())
  
  ### Data wrangling
  # I am assuming for now that the input data has rsids
  # This will need some extra functions as I receive new data from users
  dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37
  
  
  data =  colochelpR::convert_rs_to_loc(df = data, SNP_column = "SNP", dbSNP = dbSNP) %>%
    separate(loc, into = c("CHR", "POS")) %>%
    mutate(CHR = as.numeric(as.character(CHR)),
           POS = as.numeric(POS))
  
  # Perform some processing
  # We get CHRPOS for each SNP, and perform some further renaming
  # Run here the SS_CHECK FUNCTION
  metadata = SS_check(data, metadata)
} else {
  metadata = SS_check(data, metadata)
}

# Update metada file.
metadata = c(metadata[2], metadata[4])

# Get the top_SNPs file
# data_qc = fread(metadata[2])
# lead_variants = gwas_lead_snps(data_qc, pval_thresh)
# top_SNPs = make_topSNPs(lead_variants, build = "hg19",
#              write.out = T,
#              .metadata_file = metadata,
#              custom_gene = "ACP6")

# Test for non significant loci
data_qc = fread(metadata[2])
lead_variants = gwas_lead_snps(data_qc, pval_thres = 5e-7)
top_SNPs = make_topSNPs(lead_variants, build = "hg19",
                        write.out = T,
                        .metadata_file = metadata,
                        custom_gene = c("LINC00511"))

## At this point
# We have the SS GWAS updates if needed, RS updated if needed, top_SNPs file created
# We are ready to run finemap_loci
source("/mnt/rreal/RDS/acarrasco/TOOLS/echor_wrapperPipes/R/finemap_plotting.R")
finemapping_wrapper(top_SNPs = top_SNPs,
                    study_name = "PD_GWAS",
                    study_type = "motor_progression",
                    build = "hg19",
                    finemap_tools = c("ABF", "FINEMAP", "SUSIE", "POLYFUN_SUSIE"),
                    ld_ref = c("UKB"), #c("1KGphase3"),
                    mean_SS = 3500)