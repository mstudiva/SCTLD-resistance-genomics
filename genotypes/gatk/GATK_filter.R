#### packages ####

library('gridExtra')
library('ggplot2')
library('tidyverse')


#### Data import ####

rad_snp_filtered <- read.csv('ofav_2brad_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") # 2bRAD SNPs only (2bRAD samples don't contain indels)
wgs_snp_filtered <- read.csv('ofav_wgs_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") # WGS SNPs
wgs_indel_filtered <- read.csv('ofav_wgs_indel_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") # WGS indels 

rad_snp_filtered$Pass <- factor(rad_snp_filtered$FILTER == "PASS")
wgs_snp_filtered$Pass <- factor(wgs_snp_filtered$FILTER == "PASS")
wgs_indel_filtered$Pass <- factor(wgs_indel_filtered$FILTER == "PASS")

rad_snp_filtered %>%
  filter(FILTER == "PASS") -> rad_snp_passing
wgs_snp_filtered %>%
  filter(FILTER == "PASS") -> wgs_snp_passing
wgs_indel_filtered %>%
  filter(FILTER == "PASS") -> wgs_indel_passing


#### Pass/fail 2bRAD SNPs ####

snps <- '#A9E2E4'

# combined depth per variant across samples
DP_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
  xlim(0,5000)

# variant confidence standardized by depth
QD_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=1)

# Phred-scaled probability that there is strand bias at the site
FS_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60), linewidth=1, colour = c(snps)) + scale_x_continuous(trans='log10') 

# sequencing bias in which one DNA strand is favored over the other
SOR_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=SOR, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(3), linewidth=1, colour = c(snps))

# mapping quality of a variant
MQ_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=MQ, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=1, colour = c(snps))

# rank sum test for mapping qualities
MQRankSum_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=MQRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-12.5, linewidth=1, colour = c(snps))

# do all the reads support a SNP call tend to be near the end of a read?
ReadPosRankSum_rad_snp_passfail <- ggplot(rad_snp_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-8), linewidth=1, colour = c(snps))

# same as above, just with an x axis limit
ReadPosRankSum_rad_snp_passfail_xlim <- ggplot(rad_snp_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-8), linewidth=1, colour = c(snps)) + xlim(-5,5)

pdf("ofav_rad_snp_passfail.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(DP_rad_snp_passfail, QD_rad_snp_passfail, FS_rad_snp_passfail, SOR_rad_snp_passfail, MQ_rad_snp_passfail, MQRankSum_rad_snp_passfail, ReadPosRankSum_rad_snp_passfail, ReadPosRankSum_rad_snp_passfail_xlim, nrow=4)
dev.off()


#### Pass/fail WGS SNPs ####

snps <- '#A9E2E4'

# combined depth per variant across samples
DP_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
  xlim(0,5000)

# variant confidence standardized by depth
QD_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=1)

# Phred-scaled probability that there is strand bias at the site
FS_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60), linewidth=1, colour = c(snps)) + scale_x_continuous(trans='log10') 

# sequencing bias in which one DNA strand is favored over the other
SOR_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=SOR, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(3), linewidth=1, colour = c(snps))

# mapping quality of a variant
MQ_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=MQ, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=1, colour = c(snps))

# rank sum test for mapping qualities
MQRankSum_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=MQRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-12.5, linewidth=1, colour = c(snps))

# do all the reads support a SNP call tend to be near the end of a read?
ReadPosRankSum_wgs_snp_passfail <- ggplot(wgs_snp_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-8), linewidth=1, colour = c(snps))

# same as above, just with an x axis limit
ReadPosRankSum_wgs_snp_passfail_xlim <- ggplot(wgs_snp_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-8), linewidth=1, colour = c(snps)) + xlim(-5,5)

pdf("ofav_wgs_snp_passfail.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(DP_wgs_snp_passfail, QD_wgs_snp_passfail, FS_wgs_snp_passfail, SOR_wgs_snp_passfail, MQ_wgs_snp_passfail, MQRankSum_wgs_snp_passfail, ReadPosRankSum_wgs_snp_passfail, ReadPosRankSum_wgs_snp_passfail_xlim, nrow=4)
dev.off()


#### Pass/fail WGS indels ####

indels <- '#F4CCCA'

# combined depth per variant across samples
DP_wgs_indel_passfail <- ggplot(wgs_indel_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
  xlim(0,5000)

# variant confidence standardized by depth
QD_wgs_indel_passfail <- ggplot(wgs_indel_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=1)

# Phred-scaled probability that there is strand bias at the site
FS_wgs_indel_passfail <- ggplot(wgs_indel_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(200), linewidth=1, colour = c(indels)) + scale_x_continuous(trans='log10') 

# do all the reads support a indel call tend to be near the end of a read?
ReadPosRankSum_wgs_indel_passfail <- ggplot(wgs_indel_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-20), linewidth=1, colour = c(indels))

# same as above, just with an x axis limit
# ReadPosRankSum_wgs_indel_passfail_xlim <- ggplot(wgs_indel_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
   # geom_vline(xintercept=c(-20), linewidth=1, colour = c(indels)) + xlim(-5,5)

pdf("ofav_wgs_indel_passfail.pdf", height=10, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(DP_wgs_indel_passfail, QD_wgs_indel_passfail, FS_wgs_indel_passfail, ReadPosRankSum_wgs_indel_passfail, nrow=2)
dev.off()


#### Filter validation ####

# Just a sanity check to make sure all the filters performed as intended
# Each line below should return a value of 0
# If not, rerun the filtering with the filters in a different order

# 2bRAD SNPs
sum(na.omit(rad_snp_passing$QD) < 2) # 0
sum(na.omit(rad_snp_passing$QUAL) < 30) # 0
sum(na.omit(rad_snp_passing$SOR) > 3) # 0
sum(na.omit(rad_snp_passing$FS) > 60) # 0
sum(na.omit(rad_snp_passing$MQ) < 40) # 0
sum(na.omit(rad_snp_passing$MQRankSum) < -12.5) # 0
sum(na.omit(rad_snp_passing$ReadPosRankSum) < -8) # 0

# WGS SNPs
sum(na.omit(wgs_snp_passing$QD) < 2) # 0
sum(na.omit(wgs_snp_passing$QUAL) < 30) # 0
sum(na.omit(wgs_snp_passing$SOR) > 3) # 0
sum(na.omit(wgs_snp_passing$FS) > 60) # 0
sum(na.omit(wgs_snp_passing$MQ) < 40) # 0
sum(na.omit(wgs_snp_passing$MQRankSum) < -12.5) # 0
sum(na.omit(wgs_snp_passing$ReadPosRankSum) < -8) # 0

# WGS indels
sum(na.omit(wgs_indel_passing$QD) < 2) # 0
sum(na.omit(wgs_indel_passing$QUAL) < 30) # 0
sum(na.omit(wgs_indel_passing$FS) > 200) # 0
sum(na.omit(wgs_indel_passing$ReadPosRankSum) < -20) # 0


#### Filter modification ####

# All things considered, the default filtering thresholds recommended by GATK best practices (https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2)
# seem to do an adequate job removing low-confidence variants, while still retaining high-quality
# SNPs. Return to GATK_processing_README.txt for further processing.
