#### packages ####

# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("LEA")
library(LEA)

# install.packages('multidplyr')
library(multidplyr)
library(parallel)

# setting up a parallelized work environment to use multiple CPU cores
cluster <- new_cluster(parallel::detectCores()-1)
cluster_library(cluster, c('LEA'))

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
vcf_2brad_noclones <- read.vcfR("ofav_2brad_snp_passing_noclones.vcf.gz") 

# Convert vcf to gds
geno_2brad_noclones <- extract.gt(vcf_2brad_noclones, as.numeric = T) 
snpgdsCreateGeno("ofav_2brad_snp_passing_noclones.gds", geno_2brad_noclones)

# Creating a genind object for poppr
genlight_2brad_noclones <- vcfR2genlight(vcf_2brad_noclones) 

# Converting to snpclone object
snpclone_2brad_noclones <- poppr::as.snpclone(genlight_2brad_noclones) 


#### Dissimilarity matrix ####

dist_2brad_noclones <- bitwise.dist(snpclone_2brad_noclones)
write.csv(as.matrix(dist_2brad_noclones), file = "ofav_2brad_snp_passing_noclones_dist.csv")


#### data conversion ####

# Run once then comment out
# ped2lfmm("ofav_wgs_gwas.ped", force = T)
# lfmm2geno("ofav_wgs_gwas.lfmm", force = T)

# ped2lfmm("ofav_2brad_snp_passing_noclones.ped", force = T)
# lfmm2geno("ofav_2brad_snp_passing_noclones.lfmm", force = T)


#### PCA ####

# running the PCA, then comment out
# pc = pca("ofav_wgs_gwas.lfmm", scale = TRUE)
pc_2brad = pca("ofav_2brad_snp_passing_noclones.lfmm", scale = TRUE)

# loading in a PCA that has already been run
pc = load.pcaProject("ofav_wgs_gwas.pcaProject")

# perform Tracy-Widom tests for eigenvalues
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component
pdf(file = "PCA_variance.pdf")
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
dev.off()

# write PC variance to csv
write.csv(tw$percentage, file = "PCA_variance.csv")

pdf(file = "PCA.pdf", width = 14)
par(mfrow=c(1,2))
# PC1-PC2 plot
plot(pc$projections, pch=19, xlab="PC 1 (14.0% variation)",ylab="PC 1 (7.3% variation)")
# PC3-PC4 plot 
plot(pc$projections[,3:4], pch=19, xlab="PC 3 (5.5% variation)",ylab="PC 4 (4.3% variation)")
dev.off()

# looks like 5 major sample clusters

par(mfrow=c(1,1))
# Plot standard deviations.
pdf(file = "PCA_stdev.pdf")
plot(pc$sdev)
dev.off()


#### Structure ####

# running the structure analysis, then comment out
# structure = snmf("ofav_wgs_gwas.geno",
#                K = 1:10,
#                entropy = TRUE,
#                repetitions = 10,
#                project = "new",
#                CPU = 12)

# loading in a structure project that has already been run
structure = load.snmfProject("ofav_wgs_gwas.snmfProject")

# plot cross-entropy criterion for all runs in the snmf project
pdf(file = "Structure_Kselection.pdf")
plot(structure, col = "blue", pch = 19, cex = 1.2)
dev.off()
# there's no clear break to indicate the best K to use, but there is an 'elbow' at K=5, and 5 sample clusters in the PCA

# select the run with the lowest cross-entropy
best = which.min(cross.entropy(structure, K = 5))
# run 6

# set a custom color palette
my.colors <- c("#762a83", "#af8dc3","#d9f0d3", "#7fbf7b", "#1b7837")

# creating structure barplot
pdf(file = "Structure.pdf", width = 16)
barchart(structure, K = 5, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)
dev.off()

# this is the estimated membership proportion of each sample to each of the 5 ancestral lineages
res = Q(structure, K = 5, run = best)
# cluster assignment for each individual
cluster <- apply(res, 1, which.max)
res <- cbind(res,cluster)

write.csv(res, file = "ofav_wgs_gwas_structure.csv")
# now we can go back and recolor the points of the PCA to correspond to the dominant lineage for each sample

# Re-running PCA with Structure-based lineages
# load in metadata based on ancestral lineages below
metadata <- read.csv("ofav_wgs_gwas_metadata.csv", head=T)
metadata$source <- as.factor(metadata$source)
levels(metadata$source)

pdf(file = "PCA_lineages.pdf", width = 14)
par(mfrow=c(1,2))
# PC1-PC2 plot
plot(pc$projections, col=metadata$color, pch=c(15,16,17,18,13,20), xlab="PC 1 (14.0% variation)",ylab="PC 1 (7.3% variation)")
legend("topright", legend=c("Birthday", "Horseshoe", "Jaap", "Key West","Wonderland","unknown"), pch=c(15,16,17,18,20,13), bty="n")
# PC3-PC4 plot 
plot(pc$projections[,3:4], col=metadata$color, pch=c(15,16,17,18,13,20), xlab="PC 3 (5.5% variation)",ylab="PC 4 (4.3% variation)")
legend("topleft", legend=c("Birthday", "Horseshoe", "Jaap", "Key West","Wonderland","unknown"), pch=c(15,16,17,18,20,13), bty="n")
dev.off()


#### Population differentiation ####

# Genome scan for selection: population differentiation tests
pop = snmf.pvalues(structure,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 5)
pvalues = pop$pvalues

# plotting
pdf(file = "Population_differentiation.pdf", width = 10)
par(mfrow = c(1,2))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
dev.off()


#### Missing genotypes ####

# running the imputation of missing genotypes, then comment out
# impute(structure, "ofav_wgs_gwas.lfmm", method='mode', K=5, run=best) 


#### GWAS ####

# creating trait data matrix based on metadata, then comment out
# write.env(metadata[7:14], "ofav_wgs_gwas_traits.env")
write.env(metadata[14], "ofav_wgs_gwas_traits_resistance.env")
write.env(metadata[13], "ofav_wgs_gwas_traits_susceptibility.env")
write.env(metadata[10], "ofav_wgs_gwas_traits_D.env")

# Resistance
lfmm_resistance <- lfmm2(input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                         env = "ofav_wgs_gwas_traits_resistance.env", 
                         K = 5,
                         effect.sizes = TRUE)

pvalues_resistance <- lfmm2.test(object = lfmm_resistance, 
                                 input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                                 env = "ofav_wgs_gwas_traits_resistance.env", 
                                 genomic.control = TRUE,
                                 full = FALSE,
                                 linear = TRUE)

write.csv(pvalues_resistance$pvalues, file = "ofav_wgs_gwas_pvalues_resistance.csv")

# plot(-log10(lfmm_pvalues$pvalues), col = "grey", cex = .4, pch = 19)
# points(target, -log10(pv$pvalues[target]), col = "red")

pdf(file = "GWAS_Manhattan_resistance.pdf", width = 10)
par(mfrow = c(2,1))
hist(pvalues_resistance$pvalues, col = "lightblue")
plot(-log10(pvalues_resistance$pvalues), pch = 19, col = "blue", cex = .7)
dev.off()

# Susceptibility
lfmm_susceptibility <- lfmm2(input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                             env = "ofav_wgs_gwas_traits_susceptibility.env", 
                             K = 5,
                             effect.sizes = TRUE)

pvalues_susceptibility <- lfmm2.test(object = lfmm_susceptibility, 
                                     input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                                     env = "ofav_wgs_gwas_traits_susceptibility.env", 
                                     genomic.control = TRUE,
                                     full = FALSE,
                                     linear = TRUE)

write.csv(pvalues_susceptibility$pvalues, file = "ofav_wgs_gwas_pvalues_susceptibility.csv")

# plot(-log10(lfmm_pvalues$pvalues), col = "grey", cex = .4, pch = 19)
# points(target, -log10(pv$pvalues[target]), col = "red")

pdf(file = "GWAS_Manhattan_susceptibility.pdf", width = 10)
par(mfrow = c(2,1))
hist(pvalues_susceptibility$pvalues, col = "lightblue")
plot(-log10(pvalues_susceptibility$pvalues), pch = 19, col = "blue", cex = .7)
dev.off()

# Durusdinium
lfmm_D <- lfmm2(input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                env = "ofav_wgs_gwas_traits_D.env", 
                K = 5,
                effect.sizes = TRUE)

pvalues_D <- lfmm2.test(object = lfmm_D, 
                        input = "ofav_wgs_gwas.lfmm_imputed.lfmm", 
                        env = "ofav_wgs_gwas_traits_D.env", 
                        genomic.control = TRUE,
                        full = FALSE,
                        linear = TRUE)

write.csv(pvalues_D$pvalues, file = "ofav_wgs_gwas_pvalues_D.csv")

# plot(-log10(lfmm_pvalues$pvalues), col = "grey", cex = .4, pch = 19)
# points(target, -log10(pv$pvalues[target]), col = "red")

pdf(file = "GWAS_Manhattan_D.pdf", width = 10)
par(mfrow = c(2,1))
hist(pvalues_D$pvalues, col = "lightblue")
plot(-log10(pvalues_D$pvalues), pch = 19, col = "blue", cex = .7)
dev.off()


#### Dendrograms ####

if (!require("pacman")) install.packages("pacman")
pacman::p_load("dendextend", "ggdendro", "tidyverse")

# Reading in dissimilarity matrices
dissim_wgs_noclones <- t(as.matrix(read.csv("ofav_wgs_snp_passing_noclones_dist.csv", row.names = 1, head = T)))
dissim_2brad_noclones <- t(as.matrix(read.csv("ofav_2brad_snp_passing_noclones_dist.csv", row.names = 1, head = T)))

# Reading in metadata
wgs_noclones_metadata = read.csv("bams_wgs_noclones_metadata.csv", header = T) # WGS 
rad_noclones_metadata = read.csv("bams_2brad_noclones_metadata.csv", header = T) # 2bRAD 

# Creating row and column names for dissimilarity matrix
dimnames(dissim_wgs_noclones) = list(wgs_noclones_metadata[,1],wgs_noclones_metadata[,1]) # WGS 
wgs_noclones_hclust = hclust(as.dist(dissim_wgs_noclones),"ave")

dimnames(dissim_2brad_noclones) = list(rad_noclones_metadata[,1],rad_noclones_metadata[,1]) # WGS 
rad_noclones_hclust = hclust(as.dist(dissim_2brad_noclones),"ave")

# Coding factors for tree plotting
Rgroup = wgs_noclones_metadata$Rgroup # WGS 
Rgroup5clusters = wgs_noclones_metadata$Rgroup5clusters
DomSym = wgs_noclones_metadata$domSym
Lineage = as.factor(wgs_noclones_metadata$K)
Type = as.factor(wgs_noclones_metadata$type)
Source = as.factor(wgs_noclones_metadata$source)
Hybrid = as.factor(wgs_noclones_metadata$hybrid)

# Creating dendogram
wgs_noclones_dend = dissim_wgs_noclones %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram() # WGS 
wgs_noclones_data = wgs_noclones_dend %>% dendro_data()

rad_noclones_dend = dissim_2brad_noclones %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram() # WGS 
rad_noclones_data = rad_noclones_dend %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
wgs_noclones_data$segments$yend2 = wgs_noclones_data$segments$yend # WGS 
for(i in 1:nrow(wgs_noclones_data$segments)) {
  if (wgs_noclones_data$segments$yend2[i] == 0) {
    wgs_noclones_data$segments$yend2[i] = (wgs_noclones_data$segments$y[i] - 0.01)}}

rad_noclones_data$segments$yend2 = rad_noclones_data$segments$yend # WGS 
for(i in 1:nrow(rad_noclones_data$segments)) {
  if (rad_noclones_data$segments$yend2[i] == 0) {
    rad_noclones_data$segments$yend2[i] = (rad_noclones_data$segments$y[i] - 0.01)}}

wgs_noclones_dendpoints = wgs_noclones_data$labels
wgs_noclones_dendpoints$Rgroup = Rgroup[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$Rgroup5clusters=Rgroup5clusters[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$domSym=DomSym[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$lineage=Lineage[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$type=Type[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$source=Source[order.dendrogram(wgs_noclones_dend)]
wgs_noclones_dendpoints$hybrid=Hybrid[order.dendrogram(wgs_noclones_dend)]
rownames(wgs_noclones_dendpoints) = wgs_noclones_dendpoints$label

rad_noclones_dendpoints = rad_noclones_data$labels
rownames(rad_noclones_dendpoints) = rad_noclones_dendpoints$label

# Making points at the ends of branches so sample IDs line up
point_wgs = as.vector(NA)  # WGS
for(i in 1:nrow(wgs_noclones_data$segments)) {
  if (wgs_noclones_data$segments$yend[i] == 0) {
    point_wgs[i] = wgs_noclones_data$segments$y[i] - 0.01
  } else {
    point_wgs[i] = NA}}
wgs_noclones_dendpoints$y = point_wgs[!is.na(point_wgs)]

point_2brad = as.vector(NA)  # WGS
for(i in 1:nrow(rad_noclones_data$segments)) {
  if (rad_noclones_data$segments$yend[i] == 0) {
    point_2brad[i] = rad_noclones_data$segments$y[i] - 0.01
  } else {
    point_2brad[i] = NA}}
rad_noclones_dendpoints$y = point_2brad[!is.na(point_2brad)]

# plotting dendrogram
wgs_noclones_dendrogram = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y), size = 4, stroke = 0.25) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram = wgs_noclones_dendrogram + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram
ggsave("wgs_noclones_dendrogram.pdf", plot = wgs_noclones_dendrogram, height = 6, width = 24, units = "in", dpi = 300) # WGS 

rad_noclones_dendrogram = ggplot() +  # 2bRAD 
  geom_segment(data = segment(rad_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = rad_noclones_dendpoints, aes(x = x, y = y), size = 4, stroke = 0.25) +
  geom_text(data = rad_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

rad_noclones_dendrogram = rad_noclones_dendrogram + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
rad_noclones_dendrogram
ggsave("2brad_noclones_dendrogram.pdf", plot = rad_noclones_dendrogram, height = 6, width = 24, units = "in", dpi = 300) # 2bRAD 

# 4 susceptibility groups
wgs_noclones_dendrogram_4clusters = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=Rgroup), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("#d7191c","#2b83ba","#abdda4", "#fdae61")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_4clusters = wgs_noclones_dendrogram_4clusters + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_4clusters
ggsave("wgs_noclones_dendrogram_4clusters.pdf", plot = wgs_noclones_dendrogram_4clusters, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_4clusters_nolegend = wgs_noclones_dendrogram_4clusters + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_4clusters_nolegend
ggsave("wgs_noclones_dendrogram_4clusters_nolegend.pdf", plot = wgs_noclones_dendrogram_4clusters_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# 5 susceptibility groups
wgs_noclones_dendrogram_5clusters = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=Rgroup5clusters), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("#abdda4","#d7191c","#2b83ba","#1a9641", "#fdae61")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_5clusters = wgs_noclones_dendrogram_5clusters + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_5clusters
ggsave("wgs_noclones_dendrogram_5clusters.pdf", plot = wgs_noclones_dendrogram_5clusters, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_5clusters_nolegend = wgs_noclones_dendrogram_5clusters + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_5clusters_nolegend
ggsave("wgs_noclones_dendrogram_5clusters_nolegend.pdf", plot = wgs_noclones_dendrogram_5clusters_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# dominant algal symbionts
wgs_noclones_dendrogram_domSym = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=domSym), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("yellow","black")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_domSym = wgs_noclones_dendrogram_domSym + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_domSym
ggsave("wgs_noclones_dendrogram_domSym.pdf", plot = wgs_noclones_dendrogram_domSym, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_domSym_nolegend = wgs_noclones_dendrogram_domSym + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_domSym_nolegend
ggsave("wgs_noclones_dendrogram_domSym_nolegend.pdf", plot = wgs_noclones_dendrogram_domSym_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# ancestral lineages (K)
wgs_noclones_dendrogram_K = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=lineage), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("#762a83", "#af8dc3","#d9f0d3", "#7fbf7b", "#1b7837","black")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_K = wgs_noclones_dendrogram_K + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_K
ggsave("wgs_noclones_dendrogram_K.pdf", plot = wgs_noclones_dendrogram_K, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_K_nolegend = wgs_noclones_dendrogram_K + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_K_nolegend
ggsave("wgs_noclones_dendrogram_K_nolegend.pdf", plot = wgs_noclones_dendrogram_K_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# genotype type
wgs_noclones_dendrogram_type = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=type), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("cornflowerblue", "orange2","black")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_type = wgs_noclones_dendrogram_type + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_type
ggsave("wgs_noclones_dendrogram_type.pdf", plot = wgs_noclones_dendrogram_type, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_type_nolegend = wgs_noclones_dendrogram_type + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_type_nolegend
ggsave("wgs_noclones_dendrogram_type_nolegend.pdf", plot = wgs_noclones_dendrogram_type_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# genotype source
wgs_noclones_dendrogram_source = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=source), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("#01665e", "#5ab4ac","#f6e8c3","#d8b365","black","#8c510a")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_source = wgs_noclones_dendrogram_source + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_source
ggsave("wgs_noclones_dendrogram_source.pdf", plot = wgs_noclones_dendrogram_source, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_source_nolegend = wgs_noclones_dendrogram_source + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_source_nolegend
ggsave("wgs_noclones_dendrogram_source_nolegend.pdf", plot = wgs_noclones_dendrogram_source_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 

# genotype hybrid
wgs_noclones_dendrogram_hybrid = ggplot() +  # WGS 
  geom_segment(data = segment(wgs_noclones_data), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = wgs_noclones_dendpoints, aes(x = x, y = y, color=hybrid), size = 4, stroke = 0.25) +
  scale_color_manual(values=c("red","black","blue")) +
  geom_text(data = wgs_noclones_dendpoints, aes(x = x, y = (y - .0065), label = label), angle = 90) +
  labs(y = "Genetic distance") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

wgs_noclones_dendrogram_hybrid = wgs_noclones_dendrogram_hybrid + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "right")
wgs_noclones_dendrogram_hybrid
ggsave("wgs_noclones_dendrogram_hybrid.pdf", plot = wgs_noclones_dendrogram_hybrid, height = 6, width = 24, units = "in", dpi = 300) # WGS 

wgs_noclones_dendrogram_hybrid_nolegend = wgs_noclones_dendrogram_hybrid + theme(  
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
  # legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "none")
wgs_noclones_dendrogram_hybrid_nolegend
ggsave("wgs_noclones_dendrogram_hybrid_nolegend.pdf", plot = wgs_noclones_dendrogram_hybrid_nolegend, height = 6, width = 24, units = "in", dpi = 300) # WGS 
