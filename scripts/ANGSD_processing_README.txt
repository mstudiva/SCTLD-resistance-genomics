
## Analysis of Next-Generation Sequencing Data (ANGSD) pipeline, version May 21, 2024
# Created by Michael Studivan (studivanms@gmail.com) based on 2bRAD pipeline by Ryan Eckert (reckert2017@fau.edu)
https://ryaneckert.github.io/Stephanocoenia_FKNMS_PopGen/code/
WGS pipeline developed from GATK best practices
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-


#------------------------------
## Downloading scripts

cd ~/bin
svn checkout https://github.com/RyanEckert/Stephanocoenia_FKNMS_PopGen/trunk/scripts/
mv scripts/* .
rm -rf scripts

wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
unzip PGDSpider_2.0.7.1.zip
rm PGDSpider_2.0.7.1.zip

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py


#------------------------------
## Fuzzy Genotyping (ANGSD)

copy all *.bam* files to a new directory 'ANGSD'
mkdir project/directory/ANGSD
mv project/directory/mappedReads/*.bam* .

ls *bam >bamsClones

module load angsd-0.933-gcc-9.2.0-65d64pp

# Initial run ANGSD settings:
# -minMapQ 20: only highly unique mappings (prob of erroneous mapping =< 1%)
# -maxDepth: highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd: the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.
export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1850 -minInd 92"
export TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

echo '#!/bin/bash' >radDD.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out dd -nThreads 20 >>radDD.sh

sbatch --mem=200GB -o radDD.o%j -e radDD.e%j --mail-user=studivanms@gmail.com --mail-type=ALL --partition=shortq7 radDD.sh

module load R/3.6.1
echo '#!/bin/bash' >RQC.sh
echo Rscript ~/bin/plotQC.R prefix=dd >>RQC.sh
echo gzip -9 dd.counts >>RQC.sh
sbatch -e RQC.e%j -o RQC.o%j --dependency=afterok:460550 --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com RQC.sh

cat quality.txt # proportion of sites covered at >5X

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ, -minIndDepth and -minInd filters in subsequent ANGSD runs
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/\*.pdf .

# Second ANGSD run settings:
module load angsd-0.933-gcc-9.2.0-65d64pp
# change the -minInd flag to 80% of the number of samples you have
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 148 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

echo '#!/bin/bash' > radClones.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out radClones >>radClones.sh

sbatch --mem=200GB -o radClones.o%j -e radClones.e%j --mail-type=ALL --mail-user=studivanms@gmail.com --partition=shortq7 radClones.sh

# scp ibs matrix and bcf file to local machine
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/radClones.ibsMat .
scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/radClones.bcf .

# Convert bcf to vcf with bcftools
bcftools convert radClones.bcf -O z -o radClones.vcf.gz

#------------------------------
# WGS-specific settings

module load angsd-0.933-gcc-9.2.0-65d64pp

export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1800 -minInd 90"
export TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo '#!/bin/bash' >wgsDD.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out dd -nThreads 20 >>wgsDD.sh
sbatch --mem=200GB -o wgsDD.o%j -e wgsDD.e%j --mail-user=studivanms@gmail.com --mail-type=ALL --partition=longq7 wgsDD.sh

module load R/3.6.1
echo '#!/bin/bash' >RQC.sh
echo 'module load R/3.6.1' >>RQC.sh
echo Rscript ~/bin/plotQC.R prefix=dd >>RQC.sh
sbatch --partition=shortq7 -o RQC.o%j -e RQC.e%j --dependency=afterok:460550 --constraint="epyc7702" --mem=0 --mail-type=ALL --mail-user=studivanms@gmail.com RQC.sh
# --constraint="epyc7702" --mem=0 specifies a node with 1Tb memory, and allows use of all the memory

echo '#!/bin/bash' >zip.sh
echo gzip -9 dd.counts >>zip.sh
sbatch --partition=longq7 -o zip.o%j -e zip.e%j --mail-type=ALL --mail-user=studivanms@gmail.com zip.sh

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 144 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
echo '#!/bin/bash' > wgsClones.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out wgsClones >>wgsClones.sh
sbatch --mem=200GB -o wgsClones.o%j -e wgsClones.e%j --mail-type=ALL --mail-user=studivanms@gmail.com --partition=longq7 wgsClones.sh

scp mstudiva@koko-login.hpc.fau.edu:~/project/directory/ANGSD/wgsClones.ibsMat .


#------------------------------
# Next, run ANGSD_clones.R script
