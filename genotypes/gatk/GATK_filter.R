library('gridExtra')
library('ggplot2')
library('tidyverse')


#### Data import ####

# snp_filtered <- read.csv('ofav_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") 
# indel_filtered <- read.csv('ofav_indel_filtered.table', header = T, na.strings=c("","NA"), sep = "\t")

# VCF_filtered <- rbind(snp_filtered, indel_filtered)
# VCF_filtered$Variant <- factor(c(rep("SNPs", dim(snp_filtered)[1]), rep("Indels", dim(indel_filtered)[1])))

rad_snp_filtered <- read.csv('ofav_2brad_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") 

# VCF_filtered %>%
#   filter(FILTER == "PASS") -> VCF_passing

# snp_filtered$Pass <- factor(snp_filtered$FILTER == "PASS")
# indel_filtered$Pass <- factor(indel_filtered$FILTER == "PASS")

rad_snp_filtered$Pass <- factor(rad_snp_filtered$FILTER == "PASS")

# snp_filtered %>%
#   filter(FILTER == "PASS") -> snp_passing

rad_snp_filtered %>%
  filter(FILTER == "PASS") -> rad_snp_passing

# indel_filtered %>%
#   filter(FILTER == "PASS") -> indel_passing


#### Raw variants ####

snps <- '#A9E2E4'
# indels <- '#F4CCCA'

# DP <- ggplot(VCF_filtered, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) +
#   xlim(0,5000)
# 
# QD <- ggplot(VCF_filtered, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=2, linewidth=1) 
# 
# FS <- ggplot(VCF_filtered, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(60, 200), linewidth=1, colour = c(snps,indels)) + scale_x_continuous(trans='log10') 
# 
# SOR <- ggplot(VCF_filtered, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(3), linewidth=1, colour = c(snps)) 
# 
# MQ <- ggplot(VCF_filtered, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=40, linewidth=1, colour = c(snps)) 
# 
# MQRankSum <- ggplot(VCF_filtered, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=-12.5, linewidth=1, colour = c(snps))
# 
# ReadPosRankSum <- ggplot(VCF_filtered, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-8,-20), linewidth=1, colour = c(snps,indels)) 
# 
# ReadPosRankSum_xlim <- ggplot(VCF_filtered, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-8,-20), linewidth=1, colour = c(snps,indels)) + xlim(-5,5)
# 
# pdf("Ofav_variants.pdf", height=20, width=15)
# theme_set(theme_gray(base_size = 18))
# grid.arrange(DP, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, ReadPosRankSum_xlim, nrow=4)
# dev.off()


#### Filtered variants ####

# DP_filtered <- ggplot(VCF_passing, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) + 
#   xlim(0,5000)
# 
# QD_filtered <- ggplot(VCF_passing, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=2, linewidth=1) 
# 
# FS_filtered <- ggplot(VCF_passing, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(60, 200), linewidth=1, colour = c(snps,indels)) + scale_x_continuous(trans='log10') 
# 
# SOR_filtered <- ggplot(VCF_passing, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(3), linewidth=1, colour = c(snps)) 
# 
# MQ_filtered <- ggplot(VCF_passing, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=40, linewidth=1, colour = c(snps)) 
# 
# MQRankSum_filtered <- ggplot(VCF_passing, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=-12.5, linewidth=1, colour = c(snps))
# 
# ReadPosRankSum_filtered <- ggplot(VCF_passing, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-8,-20), linewidth=1, colour = c(snps,indels)) 
# 
# ReadPosRankSum_filtered_xlim <- ggplot(VCF_passing, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-8,-20), linewidth=1, colour = c(snps,indels)) + xlim(-5,5)
# 
# pdf("Ofav_variants_filtered.pdf", height=20, width=15)
# theme_set(theme_gray(base_size = 18))
# grid.arrange(DP_filtered, QD_filtered, FS_filtered, SOR_filtered, MQ_filtered, MQRankSum_filtered, ReadPosRankSum_filtered, ReadPosRankSum_filtered_xlim, nrow=4)
# dev.off()

# side-by-side of raw and filtered variants

# pdf("Ofav_variants_raw_vs_filtered.pdf", height=30, width=15)
# theme_set(theme_gray(base_size = 18))
# grid.arrange(QD, QD_filtered, FS, FS_filtered, SOR, SOR_filtered, MQ, MQ_filtered, MQRankSum, MQRankSum_filtered, ReadPosRankSum, ReadPosRankSum_filtered, nrow=6)
# dev.off()


#### Pass/fail 2bRAD SNPs ####

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

pdf("Ofav_rad_snp_passfail.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(DP_rad_snp_passfail, QD_rad_snp_passfail, FS_rad_snp_passfail, SOR_rad_snp_passfail, MQ_rad_snp_passfail, MQRankSum_rad_snp_passfail, ReadPosRankSum_rad_snp_passfail, ReadPosRankSum_rad_snp_passfail_xlim, nrow=4)
dev.off()


# #### Pass/fail indels ####
# 
# DP_indel_passfail <- ggplot(indel_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
#   xlim(0,5000)
# 
# QD_indel_passfail <- ggplot(indel_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=2, linewidth=1)
# 
# QD_indel_passfail_ylim <- ggplot(indel_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=2, linewidth=1) +ylim(0,0.05)
# 
# FS_indel_passfail <- ggplot(indel_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(200), linewidth=1, colour = c(indels)) + scale_x_continuous(trans='log10') 
# 
# # the SOR filter is not used for indels, but including for reference
# SOR_indel_passfail <- ggplot(indel_filtered, aes(x=SOR, fill=Pass)) + geom_density(alpha=.3) 
# 
# # the MQ filter is not used for indels, but including for reference
# MQ_indel_passfail <- ggplot(indel_filtered, aes(x=MQ, fill=Pass)) + geom_density(alpha=.3)
# 
# # the MQRankSum filter is not used for indels, but including for reference
# MQRankSum_indel_passfail <- ggplot(indel_filtered, aes(x=MQRankSum, fill=Pass)) + geom_density(alpha=.3) 
# 
# ReadPosRankSum_indel_passfail <- ggplot(indel_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-20), linewidth=1, colour = c(indels)) 
# 
# ReadPosRankSum_indel_passfail_xlim <- ggplot(indel_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
#   geom_vline(xintercept=c(-20), linewidth=1, colour = c(indels)) + xlim(-5,5)
# 
# pdf("Ofav_indel_passfail.pdf", height=15, width=15)
# theme_set(theme_gray(base_size = 18))
# grid.arrange(DP_indel_passfail, FS_indel_passfail, QD_indel_passfail, QD_indel_passfail_ylim, ReadPosRankSum_indel_passfail, ReadPosRankSum_indel_passfail_xlim, nrow=3)
# dev.off()


#### Filter validation ####

# Just a sanity check to make sure all the filters performed as intended
# Each line below should return a value of 0
# If not, rerun the filtering with the filters in a different order

# SNPs
sum(na.omit(rad_snp_passing$QD) < 2) # 0
sum(na.omit(rad_snp_passing$QUAL) < 30) # 0
sum(na.omit(rad_snp_passing$SOR) > 3) # 0
sum(na.omit(rad_snp_passing$FS) > 60) # 0
sum(na.omit(rad_snp_passing$MQ) < 40) # 0
sum(na.omit(rad_snp_passing$MQRankSum) < -12.5) # 0
sum(na.omit(rad_snp_passing$ReadPosRankSum) < -8) # 0

# Indels
# sum(na.omit(indel_passing$QD) < 2) # 0
# sum(na.omit(indel_passing$QUAL) < 30) # 0
# sum(na.omit(indel_passing$FS) > 200) # 0
# sum(na.omit(indel_passing$ReadPosRankSum) < -20) # 0


#### Filter modification ####

# All things considered, the default filtering thresholds recommended by GATK best practices (https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2)
# seem to do an adequate job removing low-confidence variants, while still retaining high-quality
# SNPs. Return to GATK_processing_README.txt for further processing.
