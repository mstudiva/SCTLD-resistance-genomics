#### packages ####

# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("LEA")
library(LEA)

#### data conversion ####

# Run once then comment out
ped2lfmm("ofav_wgs_gwas.ped", force = T)
lfmm2geno("ofav_wgs_gwas.lfmm", force = T)


#### Missing genotypes ####

# running the imputation of missing genotypes, then comment out
impute(structure, "ofav_wgs_gwas.lfmm", method='mode', K=5, run=best)
