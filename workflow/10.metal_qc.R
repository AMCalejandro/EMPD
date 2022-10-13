#!/usr/bin/Rscript

# Loading data from stding
args = commandArgs(trailingOnly=TRUE)
metal = fread(args[1])
N = as.numeric(args[2])
## Library loading
.libPaths("/mnt/rreal/RDS/acarrasco/R_libs/")
library(data.table)
library(stringr)
library(tidyverse)


# To perfrom the harmonization QC, we used the code from -> https://github.com/AMCalejandro/postgwasQC
source("/mnt/rreal/RDS/acarrasco/TOOLS/postgwasQC/R/Utils.R")


harmonise_metal(metal, N)