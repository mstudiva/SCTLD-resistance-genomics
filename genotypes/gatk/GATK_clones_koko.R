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
vcf_wgs <- read.vcfR("ofav_wgs_snp_passing.vcf.gz") 

# Reading in list of bam files and metadata
cloneMeta_wgs = read.csv("bams_wgs.csv", header = T) 

# Convert vcf to gds
geno_wgs <- extract.gt(vcf_wgs, as.numeric = T) 
snpgdsCreateGeno("ofav_wgs_snp_passing.gds", geno_wgs)

# Creating a genind object for poppr
genlight_wgs <- vcfR2genlight(vcf_wgs) # WGS 

# Converting to snpclone object
snpclone_wgs <- poppr::as.snpclone(genlight_wgs) 


#### Multi-locus genotyping ####

# Creating a dissimilarity matrix
dist_wgs <- bitwise.dist(snpclone_wgs)
write.csv(as.matrix(dist_wgs), file = "ofav_wgs_snp_passing_dist.csv")
# Look for the highest genetic distance value between technical replicates, and plug that number into the mlg.filter code below

# Tests for best threshold of clonal detection
# Commenting out as this method set too high of a threshold and lumped together a bunch of unrelated genotypes
# pdf(file = "mlgFilter_gatk_wgs.pdf") # Saves resulting plots
# threshtest_wgs <- filter_stats(
#   snpclone_wgs,
#   distance = bitwise.dist,
#   plot = TRUE,
#   threads = 1L)
# dev.off()
# Prints the best threshold value
# print(thresh_wgs <- cutoff_predictor(threshtest_wgs$farthest$THRESHOLDS)) # 0.05288876

# mlg.filter(snpclone_wgs, threads = 1L) <- thresh_wgs # Applies the filter based on the threshold determined above
mlg.filter(snpclone_wgs, threads = 1L) <- 0.021821604 # Applies the filter based on the dissimilarity matrix

# Multi-locus genotype assignments
mlg_wgs <- slot(snpclone_wgs, "mlg")
MLG_wgs <- slot(mlg_wgs, "mlg")
MLGeno_wgs <- as.data.frame(MLG_wgs$contracted) # Convert the slot data to a data.frame
cloneMeta_wgs_mlg <- cbind(cloneMeta_wgs, MLGeno_wgs) # Adding MLGs to metadata file
write.csv(cloneMeta_wgs_mlg, file = "ofav_gatk_wgs_mlg.csv") # Writing to file
