library('gridExtra')
library('ggplot2')


#### Data import ####

snp_filtered <- read.csv('ofav_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") 
indel_filtered <- read.csv('ofav_indel_filtered.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(snp_filtered)
dim(indel_filtered)
# VCF_filtered <- rbind(snp_filtered, indel_filtered)
# VCF_filtered$Variant <- factor(c(rep("SNPs", dim(snp_filtered)[1]), rep("Indels", dim(indel_filtered)[1])))


#### Raw variants ####

snp <- read.csv('ofav_snp.table', header = T, na.strings=c("","NA"), sep = "\t") 
indel <- read.csv('ofav_indel.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(snp)
dim(indel)
VCF <- rbind(snp, indel)
VCF$Variant <- factor(c(rep("SNPs", dim(snp)[1]), rep("Indels", dim(indel)[1])))

snps <- '#A9E2E4'
indels <- '#F4CCCA'

DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200))

QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=0.7)

FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60, 200), linewidth=0.7) + ylim(0,0.1)

MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=0.7)

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-20, linewidth=0.7)

SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), linewidth=1, colour = c(snps,indels))

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10,-20,20), linewidth=1, colour = c(snps,snps,indels,indels)) + xlim(-30, 30)

pdf("Ofav_variants_raw.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()


#### Filtered variants ####

snp_filtered <- read.csv('ofav_snp_filtered.table', header = T, na.strings=c("","NA"), sep = "\t") 
indel_filtered <- read.csv('ofav_indel_filtered.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(snp_filtered)
dim(indel_filtered)
VCF_filtered <- rbind(snp_filtered, indel_filtered)
VCF_filtered$Variant <- factor(c(rep("SNPs", dim(snp_filtered)[1]), rep("Indels", dim(indel_filtered)[1])))

snps <- '#A9E2E4'
indels <- '#F4CCCA'

DP_filtered <- ggplot(VCF_filtered, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200))

QD_filtered <- ggplot(VCF_filtered, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=0.7)

FS_filtered <- ggplot(VCF_filtered, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60, 200), linewidth=0.7) + ylim(0,0.1)

MQ_filtered <- ggplot(VCF_filtered, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=0.7)

MQRankSum_filtered <- ggplot(VCF_filtered, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-20, linewidth=0.7)

SOR_filtered <- ggplot(VCF_filtered, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), linewidth=1, colour = c(snps,indels))

ReadPosRankSum_filtered <- ggplot(VCF_filtered, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10,-20,20), linewidth=1, colour = c(snps,snps,indels,indels)) + xlim(-30, 30)

pdf("Ofav_variants_filtered.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD_filtered, DP_filtered, FS_filtered, MQ_filtered, MQRankSum_filtered, SOR_filtered, ReadPosRankSum_filtered, nrow=4)
dev.off()


#### Pass/fail SNPs ####

snp_filtered$Pass <- factor(snp_filtered$FILTER == "PASS")

DP_snp_passfail <- ggplot(snp_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200))

QD_snp_passfail <- ggplot(snp_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=0.7)

FS_snp_passfail <- ggplot(snp_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60), linewidth=0.7) + ylim(0,0.1)

MQ_snp_passfail <- ggplot(snp_filtered, aes(x=MQ, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=0.7)

MQRankSum_snp_passfail <- ggplot(snp_filtered, aes(x=MQRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-20, linewidth=0.7)

SOR_snp_passfail <- ggplot(snp_filtered, aes(x=SOR, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), linewidth=1)

ReadPosRankSum_snp_passfail <- ggplot(snp_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10), linewidth=1) + xlim(-30, 30)

pdf("Ofav_snp_passfail.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD_snp_passfail, DP_snp_passfail, FS_snp_passfail, MQ_snp_passfail, MQRankSum_snp_passfail, SOR_snp_passfail, ReadPosRankSum_snp_passfail, nrow=4)
dev.off()


#### Pass/fail indels ####

indel_filtered$Pass <- factor(indel_filtered$FILTER == "PASS")

DP_indel_passfail <- ggplot(indel_filtered, aes(x=DP, fill=Pass)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200))

QD_indel_passfail <- ggplot(indel_filtered, aes(x=QD, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, linewidth=0.7)

FS_indel_passfail <- ggplot(indel_filtered, aes(x=FS, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(200), linewidth=0.7) + ylim(0,0.1)

MQ_indel_passfail <- ggplot(indel_filtered, aes(x=MQ, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, linewidth=0.7)

MQRankSum_indel_passfail <- ggplot(indel_filtered, aes(x=MQRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-20, linewidth=0.7)

SOR_indel_passfail <- ggplot(indel_filtered, aes(x=SOR, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), linewidth=1)

ReadPosRankSum_indel_passfail <- ggplot(indel_filtered, aes(x=ReadPosRankSum, fill=Pass)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-20,20), linewidth=1) + xlim(-30, 30)

pdf("Ofav_indel_passfail.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD_indel_passfail, DP_indel_passfail, FS_indel_passfail, MQ_indel_passfail, MQRankSum_indel_passfail, SOR_indel_passfail, ReadPosRankSum_indel_passfail, nrow=4)
dev.off()