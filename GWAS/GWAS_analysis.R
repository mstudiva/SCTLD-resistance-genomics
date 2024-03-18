#### packages ####

# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("LEA")
library(LEA)

#### data conversion ####

ped2lfmm("ofav_wgs_gwas.ped", force = T)
lfmm2geno("ofav_wgs_gwas.lfmm", force = T)


#### PCA ####

pc = pca("ofav_wgs_gwas.lfmm", scale = TRUE)
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component
pdf(file = "PCA_variance.pdf")
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
dev.off()

pdf(file = "PCA.pdf", width = 15)
par(mfrow=c(1,2))
# PC1-PC2 plot.
plot(pc$projections)
# PC3-PC4 plot.
plot(pc$projections[,3:4])
dev.off()

par(mfrow=c(1,1))
# Plot standard deviations.
pdf(file = "PCA_stdev.pdf")
plot(pc$sdev)
dev.off()


#### Structure ####

structure = snmf("ofav_wgs_gwas.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new",
               CPU = 8)

project = load.snmfProject("ofav_wgs_gwas.snmfProject")

# plot cross-entropy criterion for all runs in the snmf project
plot(structure, col = "blue", pch = 19, cex = 1.2)

# get the cross-entropy of each run for K = 4
ce = cross.entropy(project, K = 4)

# select the run with the lowest cross-entropy
best = which.min(ce)

# plot the best run for K = 4 (ancestry coefficients).
barplot(t(Q(project, K = 4, run = best)))

# get the cross-entropy for the 2nd run for K = 4
ce = cross.entropy(project, K = 4, run = 2)

# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 4. 
res = G(project, K = 4, run = 2)

my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")

barchart(project, K = 4, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)