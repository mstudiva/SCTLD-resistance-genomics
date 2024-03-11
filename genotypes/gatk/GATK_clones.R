#### packages ####

# install.packages("BiocManager")
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
vcf_rad <- read.vcfR("ofav_2brad_snp_passing.vcf.gz") # 2bRAD 

# Convert vcf to gds
geno_rad <- extract.gt(vcf_rad, as.numeric = T) # 2bRAD 
snpgdsCreateGeno("ofav_2brad_snp_passing.gds", geno_rad) 

# Creating a genind object for poppr
genlight_rad <- vcfR2genlight(vcf_rad) # 2bRAD 

# Converting to snpclone object
snpclone_rad <- poppr::as.snpclone(genlight_rad) 


#### Multi-locus genotyping ####

# Creating a dissimilarity matrix
dist_rad <- bitwise.dist(snpclone_rad)
write.csv(as.matrix(dist_rad), file = "ofav_2brad_snp_passing_dist.csv")
# Look for the highest genetic distance value between technical replicates, and plug that number into the mlg.filter code below

# Tests for best threshold of clonal detection
# Commenting out as this method set too high of a threshold and lumped together a bunch of unrelated genotypes
# pdf(file = "mlgFilter_gatk_2brad.pdf") # Saves resulting plots
# threshtest_rad <- filter_stats(
#   snpclone_rad,
#   distance = bitwise.dist,
#   plot = TRUE,
#   threads = 1L)
# dev.off()
# Prints the best threshold value
# print(thresh_rad <- cutoff_predictor(threshtest_rad$farthest$THRESHOLDS)) # 0.009723445

# mlg.filter(snpclone_rad, threads = 1L) <- thresh_rad # Applies the filter based on the threshold determined above
mlg.filter(snpclone_rad, threads = 1L) <- thresh_rad # Applies the filter based on the dissimilarity matrix

# Multi-locus genotype assignments
mlg_rad <- slot(snpclone_rad, "mlg") # Pulling the 'mlg' slot from the snpclone object
MLG_rad <- slot(mlg_rad, "mlg") # Pulling the 'mlg' slot from the mlg object
MLGeno_rad <- as.data.frame(MLG_rad$contracted) # Convert the slot data to a data.frame
cloneMeta_rad_mlg <- cbind(cloneMeta_rad, MLGeno_rad) # Adding MLGs to metadata file
write.csv(cloneMeta_rad_mlg, file = "ofav_gatk_2brad_mlg.csv") # Writing to file


#### Dendrogram ####

# Reading in dissimilarity matrices
cloneMa_rad <- t(as.matrix(read.csv("ofav_2brad_snp_passing_dist.csv", row.names = 1, head = T)))

cloneMa_wgs <- t(as.matrix(read.csv("ofav_wgs_snp_passing_dist.csv", row.names = 1, head = T)))

# Reading in metadata
cloneMeta_rad = read.csv("../bams_2brad.csv", header = T) # 2bRAD 

cloneMeta_wgs = read.csv("../bams_wgs.csv", header = T) # WGS 

# Creating row and column names for dissimilarity matrix
dimnames(cloneMa_rad) = list(cloneMeta_rad[,1],cloneMeta_rad[,1]) # 2bRAD 
clonesHc_rad = hclust(as.dist(cloneMa_rad),"ave")

dimnames(cloneMa_wgs) = list(cloneMeta_wgs[,1],cloneMeta_wgs[,1]) # WGS 
clonesHc_wgs = hclust(as.dist(cloneMa_wgs),"ave")

# Coding factors for tree plotting
cloneGeno_rad = cloneMeta_rad$PutGeno # 2bRAD 
cloneReps_rad = cloneMeta_rad$GenoRep
cloneSeq_rad = cloneMeta_rad$Seq

cloneGeno_wgs = cloneMeta_wgs$PutGeno # WGS 
cloneReps_wgs = cloneMeta_wgs$GenoRep
cloneSeq_wgs = cloneMeta_wgs$Seq

# Creating dendogram
cloneDend_rad = cloneMa_rad %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram() # 2bRAD 
cloneDData_rad = cloneDend_rad %>% dendro_data()

cloneDend_wgs = cloneMa_wgs %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram() # WGS 
cloneDData_wgs = cloneDend_wgs %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
cloneDData_rad$segments$yend2 = cloneDData_rad$segments$yend # 2bRAD 
for(i in 1:nrow(cloneDData_rad$segments)) {
  if (cloneDData_rad$segments$yend2[i] == 0) {
    cloneDData_rad$segments$yend2[i] = (cloneDData_rad$segments$y[i] - 0.01)}}

cloneDendPoints_rad = cloneDData_rad$labels
cloneDendPoints_rad$geno = cloneGeno_rad[order.dendrogram(cloneDend_rad)]
cloneDendPoints_rad$reps=cloneReps_rad[order.dendrogram(cloneDend_rad)]
rownames(cloneDendPoints_rad) = cloneDendPoints_rad$label

cloneDData_wgs$segments$yend2 = cloneDData_wgs$segments$yend # WGS 
for(i in 1:nrow(cloneDData_wgs$segments)) {
  if (cloneDData_wgs$segments$yend2[i] == 0) {
    cloneDData_wgs$segments$yend2[i] = (cloneDData_wgs$segments$y[i] - 0.01)}}

cloneDendPoints_wgs = cloneDData_wgs$labels
cloneDendPoints_wgs$geno = cloneGeno_wgs[order.dendrogram(cloneDend_wgs)]
cloneDendPoints_wgs$reps=cloneReps_wgs[order.dendrogram(cloneDend_wgs)]
rownames(cloneDendPoints_wgs) = cloneDendPoints_wgs$label

# Making points at the ends of branches so sample IDs line up
point_rad = as.vector(NA) # 2bRAD
for(i in 1:nrow(cloneDData_rad$segments)) {
  if (cloneDData_rad$segments$yend[i] == 0) {
    point_rad[i] = cloneDData_rad$segments$y[i] - 0.01
  } else {
    point_rad[i] = NA}}
cloneDendPoints_rad$y = point_rad[!is.na(point_rad)]

point_wgs = as.vector(NA)  # WGS
for(i in 1:nrow(cloneDData_wgs$segments)) {
  if (cloneDData_wgs$segments$yend[i] == 0) {
    point_wgs[i] = cloneDData_wgs$segments$y[i] - 0.01
  } else {
    point_wgs[i] = NA}}
cloneDendPoints_wgs$y = point_wgs[!is.na(point_wgs)]

# Reading in metadata for technical replicates (sequenced duplicates)
techReps_rad <- c(read.csv("../techReps_rad.csv", head = F)) # 2bRAD

techReps_wgs <- c(read.csv("../techReps_wgs.csv", head = F)) # WGS

# Plotting dendrogram based on genetic dissimilarity
cloneDendA_rad = ggplot() + # 2bRAD 
  geom_segment(data = segment(cloneDData_rad), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints_rad, aes(x = x, y = y, fill = geno), size = 4, stroke = 0.25) +
  # scale_shape_manual(values = c(21, 22), name = "Sequencing Pipeline")+
  # geom_hline(yintercept = 0.75, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold based on maximum genetic dissimilarity between technical replicates
  geom_text(data = subset(cloneDendPoints_rad, subset = reps %in% techReps_rad$V1), aes(x = x, y = (y - .0275), label = reps), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints_rad, subset = !reps %in% techReps_rad$V1), aes(x = x, y = (y - .0075), label = geno), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

# Applying formatting
cloneDend_rad = cloneDendA_rad + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 12, color = "black", angle = 90),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.line.y = element_line(),
  axis.ticks.y = element_line(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
cloneDend_rad
ggsave("cloneDend_gatk_2brad.pdf", plot = cloneDend_rad, height = 8, width = 28, units = "in", dpi = 300) # 2bRAD 

cloneDendA_wgs = ggplot() +  # WGS 
  geom_segment(data = segment(cloneDData_wgs), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints_wgs, aes(x = x, y = y, fill = geno), size = 4, stroke = 0.25) +
  # scale_shape_manual(values = c(21, 22), name = "Sequencing Pipeline")+
  geom_hline(yintercept = 0.014201934, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold based on maximum genetic dissimilarity between technical replicates
  geom_text(data = subset(cloneDendPoints_wgs, subset = reps %in% techReps_wgs$V1), aes(x = x, y = (y - .0275), label = reps), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints_wgs, subset = !reps %in% techReps_wgs$V1), aes(x = x, y = (y - .0075), label = geno), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

cloneDend_wgs = cloneDendA_wgs + theme(  
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 12, color = "black", angle = 90),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.line.y = element_line(),
  axis.ticks.y = element_line(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
cloneDend_wgs
ggsave("cloneDend_gatk_wgs.pdf", plot = cloneDend_wgs, height = 8, width = 28, units = "in", dpi = 300) # WGS 
