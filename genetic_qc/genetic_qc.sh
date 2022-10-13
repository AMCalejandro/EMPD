#!/bin/bash




############################################################################################################################
# This is just my standard QC pipeline for genetic data. This pipeline has to be readapted depending on the aim of the study.
############################################################################################################################

# Storing the genetic filename and the new name I will give to the output files.
FILENAME=$1
NEW_NAME=$2
SUBSET_TO_EXTRACT=$3
CLINICAL_GENDERS=$4

##########################################################################
### 1. Filter out individuals I am not interested in for the analysis ####
##########################################################################

#The first step is to extract the individuals I am interested for the analysis.
#This is because different groups may come from slightly different populations, e.g. if I am just doing a GWAS in PD cases, I don't want to include any controls (even if my study recruited controls)
#because the cases and controls may be sampled from slightly different populations. This will change the variants and individuals that you remove during your QC,
#as the distributions and means etc. will be different depending on who you incude.

plink1.9 --bfile ../RAW_DATA/$FILENAME \
--keep $SUBSET_TO_EXTRACT \
--make-bed \
--out $NEW_NAME.filtered



# Prior doing any specific processing it is good to exclude variants with high missing rate
plink1.9 --bfile OPDC.filtered.updatedsex \
--geno 0.02 \
--make-bed \
--out $NEW_NAME.filtered.geno_0.02




##################################
######### 2. SampleQC ############
##################################

#In this second step we will remove samples that do not meet specific quality checkings

#2.1 Sex checking ##############
################################
################################

# First update cinical gender in fam file

plink1.9 --bfile $NEW_NAME.filtered \
--update-sex $CLINICAL_GENDERS \
--make-bed \
--out $NEW_NAME.filtered.updatedsex

# Then I check the sex for males and famales is around the expected values
plink1.9 --bfile $NEW_NAME.filtered.updatedsex \
--check-sex 0.2 0.7 \
--out $NEW_NAME.filtered.updatedsex.sexcheck_results

# I keep then all the samples whose sex status is OK
cat $NEW_NAME.filtered.updatedsex.sexcheck_results.sexcheck | awk '($5=="OK") {print $1 "\t"$2}' > sex_samples_to_keep.txt

plink1.9 --bfile $NEW_NAME.filtered.updatedsex \
--keep sex_samples_to_keep.txt \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass


# 2.2 Missingness and heterozygosity ##
#######################################
#######################################

# With respect to heterozygosity, we remove samples that are 2SD away from the mean
# With respect to the missingness per individual, we remove samples that do not have a call rate >98%
# This cutoffs are pretty standard - Can be adjusted depending on my data/study

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass \
--missing \
--het \
--out $NEW_NAME.filtered.updatedsex.sexpass.mishet

# We generate the files having the outliers and also we plot the sample
R --vanilla --slave \
--args $NEW_NAME.filtered.updatedsex.sexpass.mishet.het $NEW_NAME.filtered.updatedsex.sexpass.mishet.imiss < ./het_callrate.R

#Then we just remove the samples
plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass \
--remove samples_to_remove.txt \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD


# 2.3 Relatedness ###############
#################################
#################################

#In this step we look at pairs of samples that are highly related (grm-cutoff 0.125)

# As a first step, we prun the data
plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD \
--geno 0.01 \
--maf 0.05 \
--indep-pairwise 50 5 0.5 --out pruning

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD \
--extract pruning.prune.in \
--make-bed \
--out pruned_data

# We generate the grm matrix which records the paired relationships of our sample
../TOOLS/gcta64 --bfile pruned_data \
--make-grm \
--out GRM_matrix \
--autosome \
--maf 0.05

# And then we select those pairs of inidividuals distantly related.
# A values on the matrix of 1 for two individuals would mean that it is the same individuals
# A value of 0.5 for a pair of individuals of 0.5 indicates parent/child relationship
../TOOLS/gcta64 --grm-cutoff 0.125 \
--grm GRM_matrix \
--out GRM_matrix_0125 \
--make-grm

# We keep those inidividuals that did not overpass the threshold
plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD \
--keep GRM_matrix_0125.grm.id \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125

# To see individuals before and after...
cut -f 1,2 $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.fam > IDs_before_relatedness_filter.txt
cut -f 1,2 $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.fam > IDs_after_relatedness_filter.txt

# To print the individuals removed, I could do something like , sort both files, count the number of times each line appeared
# to finally print only the lines that appeared once
# sort file1 file2 | uniq -c | awk '$1 == 1 {print $2}'


# 2.4 Ancestry checking #########
#################################
#################################

# As a first step, we plot the PC's cohort along with the ones from the Reference panel
# We mainly want to decide which is a good threshold with this exploration of the PC distribution

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125 \
--extract /data/kronos/kronos/NGS_Reference/HapMap_Reference/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps>
--make-bed \
--out NeuroChip.hapmap_SNPs

# Then we get the CEU, CHB, JPT, and YRI from HapMap.
plink1.9 --bfile /data/kronos/kronos/NGS_Reference/HapMap_Reference/hapmap3r2_CEU.CHB.JPT.YRI.founde>
--keep /data/kronos/kronos/acarrasco/FILES/CEU.CHB.JPT.YRI_hapmap.txt  \
--make-bed \
--out HapMapIII_CGRCh38.CHB.JPT.YRI.CEU
# Note that sometimes I will have to exclude extra SNPs because they do not match
# even after the flip


# Merge my Chip_hapmap SNPs with the hapmap CEU_CHB_YRI_JPT data and extract the pruned SNPs
# Note that the pruning.prune.in was already got on the previous step
plink1.9 --bfile NeuroChip.hapmap_SNPs \
--bmerge HapMapIII_CGRCh38.CHB.JPT.YRI.CEU \
--extract pruning.prune.in \
--make-bed \
--out NeuroChip.hapmap_SNPs.merged-pruned

# Need to flip missnps, beacuase some alleles may be swapped and that may be the reason
# why alleled did not match in my data to the hapmap data
plink1.9 --bfile NeuroChip.hapmap_SNPs \
--flip NeuroChip.hapmap_SNPs.merged-pruned-merge.missnp \
--make-bed \
--out NeuroChip.hapmap_SNPs.flipped_missnps

# Now we can remerge again the my SNPs hapmap file that now contains swapped the alleles that
#did not match with the CEU hapmap allele

plink1.9 --bfile NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps \
--bmerge hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only \
--extract pruning.prune.in \
--make-bed \
--out NeuroChip.hapmap_SNPs.flipped_missnps.merged-pruned


# Now I am ready to get the PCA plot for the HapMap populations together with my cohort
# Once I have checked I have at least 10K variants (the minimum for having a powered analysis)
# I can calculate the PCs from the merged ( my data + hapmap CEU) and pruned dataset

../TOOLS/gcta64 --bfile NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned \
--make-grm \
--autosome \
--out NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned.matrix

#Run PCA to generate 10 principal components
../../TOOLS/gcta64 --grm NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned.matrix \
--pca 10 \
--out PCA


# Then we call the R script
# Note that this script must be modifed each time
R --vanilla --slave --args PPMI PCA_Chip.HapMapPops.eigenvec < ./

# Then, once we have a visualization, we should decide the SD threshold and this time,
R --vanilla --slave --args PPMI PCA.ceu_chip.eigenvec < ./ancestry_outliers.R


plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125 \
--remove PCA_outliers.txt \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep

####################################
######### 2. VariantQC  ############
####################################

# Now we have to perform som variant QC steps.


# Missingness ########
######################
######################

# We remove those variants that have a missing rate > 0.05

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep \
--geno 0.05 \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05

# MAF ##################
########################
########################

# markers with a call rate less than 95% are removed
plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05 \
--maf 0.05 \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05

# Missingness by case control ####
##################################
##################################

#Note that in analysis that anly onvolves cases, no controls, this step will be skipped

#plink1.9 --bfile XXXXXX \
#--test-missing \
#--out missing_snps
#awk {if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
#plink1.9 --bfile XXXXXXXXX \
#--exclude missing_snps_1E4.txt \
#sort -u missing_snps_1E4.txt > VARIANT_TEST_MISSING_SNPS.txt
#--make-bed \
#--out XXXXXXXXXXX


# Missingness byhaplotype ########
##################################
##################################

# Now we test whether genotype calls at the two adjacent variants could be used to predict
# missingness status of the current variant. This can help judge the safety of assuming missing
# missing calls are randomly distributed.

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05 \
--test-mishap \
--out missing_hap

awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt

sort -u missing_haps_1E4_final.txt > HAPLOTYPE_TEST_MISSING_SNPS.txt
plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05 \
--exclude missing_haps_1E4_final.txt \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05.mishap


### HWE ##########################
##################################
##################################

# This is to check whether the variants of Control samples meet the HWE.

# Dependging on the data we work with, we will set the threshold lower or higher
# In this case, as I am only working with PD cases, I will set the threshold quite low
# just to make sure I am only removing serioud genotyping errors and not genuine SNP-trait
# assocations which deviate slightly from the HWE

plink1.9 --bfile $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05.mishap \
--hwe 1E-10 \
--make-bed \
--out $NEW_NAME.filtered.updatedsex.sexpass.sample_0.98.het_4SD.related_0.125.PCA_keep.geno_0.05.maf_0.05.mishap.hwe_1E10