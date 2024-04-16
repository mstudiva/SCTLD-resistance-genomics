#### packages ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   BiocManager::install(version = "3.18")
# BiocManager::install("SNPRelate")
library(SNPRelate)
library(gdsfmt)

# Install VariantAnnotation package if not installed
# BiocManager::install("VariantAnnotation")
library(VariantAnnotation)

# Install vcfR package if not installed
# BiocManager::install("vcfR")
library(vcfR)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("dendextend", "ggdendro", "tidyverse")

# install.packages("poppr", repos = "https://cloud.r-project.org/")
library(poppr)

# install,packages("ape")
library(ape)


#### Vcf conversion ####

# Reading in vcf files
vcf_wgs_noclones <- read.vcfR("ofav_wgs_snp_passing_noclones.vcf.gz") 

# Convert vcf to gds
geno_wgs_noclones <- extract.gt(vcf_wgs_noclones, as.numeric = T) 
snpgdsCreateGeno("ofav_wgs_snp_passing_noclones.gds", geno_wgs_noclones)

# Creating a genind object for poppr
genlight_wgs_noclones <- vcfR2genlight(vcf_wgs_noclones) # WGS 

# Converting to snpclone object
snpclone_wgs_noclones <- poppr::as.snpclone(genlight_wgs_noclones) 


#### Dissimilarity matrix ####

dist_wgs_noclones <- bitwise.dist(snpclone_wgs_noclones)
write.csv(as.matrix(dist_wgs_noclones), file = "ofav_wgs_snp_passing_noclones_dist.csv")