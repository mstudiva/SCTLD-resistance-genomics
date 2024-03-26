#### packages ####

# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("LEA")
library(LEA)

#### data conversion ####

# Run once then comment out
# ped2lfmm("ofav_wgs_gwas.ped", force = T)
# lfmm2geno("ofav_wgs_gwas.lfmm", force = T)


#### PCA ####

# running the PCA, then comment out
# pc = pca("ofav_wgs_gwas.lfmm", scale = TRUE)

# loading in a PCA that has already been run
pc = load.pcaProject("ofav_wgs_gwas.pcaProject")

# perform Tracy-Widom tests for eigenvalues
tw = tracy.widom(pc)

# load in metadata based on ancestral lineages below
metadata <- read.csv("ofav_wgs_gwas_metadata.csv", head=T)

# plot the percentage of variance explained by each component
pdf(file = "PCA_variance.pdf")
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
dev.off()

pdf(file = "PCA.pdf", width = 7)
par(mfrow=c(1,1))
# PC1-PC2 plot
plot(pc$projections, col=metadata$color, pch=19, xlab="PC 1 (13.8% variation)",ylab="PC 1 (0.1% variation)")
# PC3-PC4 plot (not plotting since the variance explained is so low)
# plot(pc$projections[,3:4], col=metadata$color, pch=19, xlab="PC 3 (0.05% variation)",ylab="PC 4 (0.04% variation)")
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
plot(structure, col = "blue", pch = 19, cex = 1.2)
# there's no clear break to indicate the best K to use, but there is an 'elbow' at K=5, and 5 sample clusters in the PCA

# select the run with the lowest cross-entropy
best = which.min(cross.entropy(structure, K = 5))

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
res = Q(structure, K = 5, run = 8)
# cluster assignment for each individual
cluster <- apply(res, 1, which.max)
res <- cbind(res,cluster)

write.csv(res, file = "ofav_wgs_gwas_structure.csv")
# now we can go back and recolor the points of the PCA to correspond to the dominant lineage for each sample


#### Population differentiation ####

# Genome scan for selection: population differentiation tests
pop = snmf.pvalues(structure,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 5)
pvalues = pop$pvalues

# plotting
pdf(file = "Population differentiation.pdf", height = 10)
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
dev.off()


#### Missing genotypes ####

# running the imputation of missing genotypes, then comment out
# impute(structure, "ofav_wgs_gwas.lfmm", method='mode', K=5, run=best) 


#### GWAS ####

# creating trait data matrix based on metadata
traits <- metadata[]

mod.lfmm2 <- lfmm2(input = "ofav_wgs_gwas.lfmm_imputed.lfmm", env = traits, K = 5)

pv <- lfmm2.test(object = mod.lfmm2, input = "ofav_wgs_gwas.lfmm_imputed.lfmm", env = traits, linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")
