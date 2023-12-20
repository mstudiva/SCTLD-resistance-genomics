#### IBS matrix ####

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

# Replace 'your_data.vcf.gz' with your actual VCF file name
vcf <- read.vcfR("ofav_2brad_snp_passing.vcf.gz") # 2bRAD only
# vcf <- read.vcfR("ofav_snp_passing.vcf.gz") # WGS and 2bRAD

# Convert vcf to GDS
geno <- extract.gt(vcf, as.numeric = T)
snpgdsCreateGeno("ofav_2brad_snp_passing.gds", geno) # 2bRAD only
# snpgdsCreateGeno("ofav_snp_passing.gds", geno) # WGS and 2bRAD

# Creating IBS matrix
# genofile<-snpgdsOpen("ofav_2brad_snp_passing.gds") # 2bRAD only
genofile<-snpgdsOpen("ofav_snp_passing.gds") # WGS and 2bRAD
set.seed(100)
ibs <- snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))
write.csv(ibs$dist, file = "ofav_2brad_snp_passing_ibs.csv") # 2bRAD only
# write.csv(ibs$dist, file = "ofav_snp_passing_ibs.csv") # WGS and 2bRAD


#### Dendrogram ####

if (!require("pacman")) install.packages("pacman")

pacman::p_load("dendextend", "ggdendro", "tidyverse")

# cloneMeta = read.csv("bams_2brad.csv", header = T) # 2bRAD only
cloneMeta = read.csv("bams.csv", header = T) # WGS and 2bRAD

cloneMa = ibs$dist

dimnames(cloneMa) = list(cloneMeta[,1],cloneMeta[,1])
clonesHc = hclust(as.dist(cloneMa),"ave")

cloneGeno = cloneMeta$PutGeno
cloneReps = cloneMeta$GenoRep
cloneSeq = cloneMeta$Seq

cloneDend = cloneMa %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram()
cloneDData = cloneDend %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
cloneDData$segments$yend2 = cloneDData$segments$yend
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend2[i] == 0) {
    cloneDData$segments$yend2[i] = (cloneDData$segments$y[i] - 0.01)}}

cloneDendPoints = cloneDData$labels
cloneDendPoints$geno = cloneGeno[order.dendrogram(cloneDend)]
cloneDendPoints$reps=cloneReps[order.dendrogram(cloneDend)]
rownames(cloneDendPoints) = cloneDendPoints$label

# Making points at the leaves to place symbols for sequencing pipeline
point = as.vector(NA)
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend[i] == 0) {
    point[i] = cloneDData$segments$y[i] - 0.01
  } else {
    point[i] = NA}}

cloneDendPoints$y = point[!is.na(point)]

techReps <- c(read.csv("techReps.csv", head = F))

cloneDendA = ggplot() +
  geom_segment(data = segment(cloneDData), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints, aes(x = x, y = y, fill = geno), size = 4, stroke = 0.25) +
  #scale_fill_brewer(palette = "Dark2", name = "Population") +
  # scale_fill_manual(values = flPal, name= "Population")+
  scale_shape_manual(values = c(21, 22), name = "Sequencing Pipeline")+
  # geom_hline(yintercept = 0.75, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold
  geom_text(data = subset(cloneDendPoints, subset = reps %in% techReps$V1), aes(x = x, y = (y - .025), label = reps), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints, subset = !reps %in% techReps$V1), aes(x = x, y = (y - .010), label = geno), angle = 90) +
  labs(y = "Genetic distance (1 - IBS)") +
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  theme_classic()

cloneDend = cloneDendA + theme(
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

cloneDend

ggsave("cloneDend_2brad.pdf", plot = cloneDend, height = 7, width = 28, units = "in", dpi = 300) # 2bRAD only
# ggsave("cloneDend_wgs_2brad.pdf", plot = cloneDend, height = 7, width = 49, units = "in", dpi = 300) # WGS and 2bRAD
