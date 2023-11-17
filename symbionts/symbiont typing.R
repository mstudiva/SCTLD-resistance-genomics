#### packages ####

library(tidyverse)
library(reshape2)
library(ggplot2)
library(RColorBrewer)


#### WGS data import ####

wgs = read.delim("../WGS/ZooxReads", header = FALSE, check.names = FALSE)

head(wgs)

#Reconstruct read mapping output into dataframe usable for analysis
wgs$V2[is.na(wgs$V2)] <- as.character(wgs$V1[is.na(wgs$V2)])
wgs = wgs %>% filter(wgs$V1 != "*")
wgsLst = split(wgs$V2, as.integer(gl(length(wgs$V2), 20, length(wgs$V2))))

wgsMaps = NULL

for(i in wgsLst){
  wgsMaps = rbind(wgsMaps, data.frame(t(i)))
}

#convert characters to numeric
str(wgsMaps)

for(i in c(2:20)){
  wgsMaps[,i] = as.numeric(wgsMaps[,i])
}

str(wgsMaps)

#collapse fake chromosomes into representativ genera
wgsMaps$Symbiodinium = rowSums(wgsMaps[2:6])
wgsMaps$Breviolum = rowSums(wgsMaps[7:10])
wgsMaps$Cladocopium = rowSums(wgsMaps[11:16])
wgsMaps$Durusdinium = rowSums(wgsMaps[17:20])

#keep genera totals and turn into proportions for barplot
wgsMaps = wgsMaps[,c(1, 21:24)]
wgsProp = wgsMaps
wgsProp$sum = apply(wgsProp[, c(2:length(wgsProp[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})
wgsProp = cbind(wgsProp$X1, (wgsProp[, c(2:(ncol(wgsProp)-1))]
                                   / wgsProp$sum))

colnames(wgsProp)[1] = "Sample"

head(wgsProp)

#Check that all samples total to 1
apply(wgsProp[, c(2:(ncol(wgsProp)))], 1, function(x) {
  sum(x, na.rm = T)
})

# exporting proportion data
write.csv(wgsProp, file="WGS_zoox.csv")

# transposing and reformatting dataframes to make abundance a single column
wgsPerc <- reshape2::melt(wgsProp, id = "Sample")
wgsPerc$variable=factor(wgsPerc$variable, levels=c("Symbiodinium","Breviolum","Cladocopium","Durusdinium")) 
wgsPerc %>% mutate(across('Sample', str_replace, '.un.bt2.bam','')) -> wgsPerc

# percent stacked barplot
wgsPlot <- ggplot(wgsPerc, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Sample",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "WGS (Experimental Genotypes") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  # geom_text(x = 1.8, y=1.035, label = "Site: F2,21 = 0.9, R2 = 0.085, p = ns") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
wgsPlot 

ggsave("WGS symbionts.pdf", plot= wgsPlot, width=64, height=4, units="in", dpi=300, limitsize = F)


#### 2bRAD data import ####

rad = read.delim("../2bRAD/ZooxReads", header = FALSE, check.names = FALSE)

head(rad)

#Reconstruct read mapping output into dataframe usable for analysis
rad$V2[is.na(rad$V2)] <- as.character(rad$V1[is.na(rad$V2)])
rad = rad %>% filter(rad$V1 != "*")
radLst = split(rad$V2, as.integer(gl(length(rad$V2), 20, length(rad$V2))))

radMaps = NULL

for(i in radLst){
  radMaps = rbind(radMaps, data.frame(t(i)))
}

#convert characters to numeric
str(radMaps)

for(i in c(2:20)){
  radMaps[,i] = as.numeric(radMaps[,i])
}

str(radMaps)

#collapse fake chromosomes into representativ genera
radMaps$Symbiodinium = rowSums(radMaps[2:6])
radMaps$Breviolum = rowSums(radMaps[7:10])
radMaps$Cladocopium = rowSums(radMaps[11:16])
radMaps$Durusdinium = rowSums(radMaps[17:20])

#keep genera totals and turn into proportions for barplot
radMaps = radMaps[,c(1, 21:24)]
radProp = radMaps
radProp$sum = apply(radProp[, c(2:length(radProp[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})
radProp = cbind(radProp$X1, (radProp[, c(2:(ncol(radProp)-1))]
                             / radProp$sum))

colnames(radProp)[1] = "Sample"

head(radProp)

#Check that all samples total to 1
apply(radProp[, c(2:(ncol(radProp)))], 1, function(x) {
  sum(x, na.rm = T)
})

# exporting proportion data
write.csv(radProp, file="2bRAD_zoox.csv")

# transposing and reformatting dataframes to make abundance a single column
radPerc <- reshape2::melt(radProp, id = "Sample")
radPerc$variable=factor(radPerc$variable, levels=c("Symbiodinium","Breviolum","Cladocopium","Durusdinium")) 
radPerc %>% mutate(across('Sample', str_replace, '.trim.un.bt2.bam','')) %>%
  mutate(across('Sample', str_replace, 'CPR_','')) -> radPerc

# percent stacked barplot
radPlot <- ggplot(radPerc, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Sample",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "2bRAD (Experimental Genotypes") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  # geom_text(x = 1.8, y=1.035, label = "Site: F2,21 = 0.9, R2 = 0.085, p = ns") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
radPlot 

ggsave("2bRAD symbionts.pdf", plot= radPlot, width=64, height=4, units="in", dpi=300, limitsize = F)
