if (!require("pacman")) install.packages("pacman")

pacman::p_load("dendextend", "ggdendro", "tidyverse")

cloneBams = read.csv("bamsClones.csv", head=F) # list of bam files

cloneMa = as.matrix(read.table("ofavClones.ibsMat")) # reads in IBS matrix produced by ANGSD 

dimnames(cloneMa) = list(cloneBams[,1],cloneBams[,1])
clonesHc = hclust(as.dist(cloneMa),"ave")

# clonePops = cloneBams$region
# cloneDepth = cloneBams$depthZone

cloneDend = cloneMa %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram()
cloneDData = cloneDend %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
cloneDData$segments$yend2 = cloneDData$segments$yend
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend2[i] == 0) {
    cloneDData$segments$yend2[i] = (cloneDData$segments$y[i] - 0.01)}}

cloneDendPoints = cloneDData$labels
# cloneDendPoints$pop = clonePops[order.dendrogram(cloneDend)]
# cloneDendPoints$depth=cloneDepth[order.dendrogram(cloneDend)]
rownames(cloneDendPoints) = cloneDendPoints$label

# Making points at the leaves to place symbols for populations
# point = as.vector(NA)
# for(i in 1:nrow(cloneDData$segments)) {
#   if (cloneDData$segments$yend[i] == 0) {
#     point[i] = cloneDData$segments$y[i] - 0.01
#   } else {
#     point[i] = NA}}

# cloneDendPoints$y = point[!is.na(point)]

techReps = c("CPR_210_2.trim.bt2.bam", "CPR_210_3.trim.bt2.bam", "CPR_210.trim.bt2.bam", "CPR_226_2.trim.bt2.bam", "CPR_226_3.trim.bt2.bam", "CPR_226.trim.bt2.bam", "CPR_375_2.trim.bt2.bam", "CPR_375_3.trim.bt2.bam", "CPR_375.trim.bt2.bam", "CPR_376_2.trim.bt2.bam", "CPR_376_3.trim.bt2.bam", "CPR_376.trim.bt2.bam", "CPR_sperm_2.trim.bt2.bam", "CPR_sperm_3.trim.bt2.bam", "CPR_sperm.trim.bt2.bam")
# cloneDendPoints$depth = factor(cloneDendPoints$depth,levels(cloneDendPoints$depth)[c(2,1)])

# cloneDendPoints$pop = factor(cloneDendPoints$pop,levels(cloneDendPoints$pop)[c(4, 1, 3, 2)])

# flPal = paletteer_d("vapoRwave::jazzCup")[c(2:5)]

cloneDendA = ggplot() +
  geom_segment(data = segment(cloneDData), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  # geom_point(data = cloneDendPoints, aes(x = x, y = y), size = 4, stroke = 0.25) +
  #scale_fill_brewer(palette = "Dark2", name = "Population") +
  # scale_fill_manual(values = flPal, name= "Population")+
  # scale_shape_manual(values = c(24, 25), name = "Depth Zone")+
  geom_hline(yintercept = 0.12, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold
  geom_text(data = subset(cloneDendPoints, subset = label %in% techReps), aes(x = x, y = (y - .015), label = label), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints, subset = !label %in% techReps), aes(x = x, y = (y - .010), label = label), angle = 90) +
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
  legend.position = "bottom")

cloneDend

ggsave("cloneDend.pdf", plot = cloneDend, height = 18, width = 35, units = "in", dpi = 300)
# ggsave("../figures/cloneDend.eps", plot = cloneDend, height = 8, width = 35, units = "in", dpi = 300)