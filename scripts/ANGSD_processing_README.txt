
## Analysis of Next-Generation Sequencing Data (ANGSD) pipeline, version November 29, 2023
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
## Fuzzy Genotyping (ANGSD; 2bRAD and WGS together)

# If working with WGS and 2bRAD samples together, copy all *.bam* files to a master directory 'ANGSD'
mkdir project/directory/ANGSD
mv project/directory/2bRAD/mappedReads/*.bam* .
mv project/directory/WGS/mappedReads/*.bam* .

ls *bam >bamsClones

module load angsd-0.933-gcc-9.2.0-65d64pp

# ANGSD settings:
# -minMapQ 20: only highly unique mappings (prob of erroneous mapping =< 1%)
# -maxDepth: highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd: the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.
export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 3650 -minInd 182"
export TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

echo '#!/bin/bash' >ofavDD.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out dd -nThreads 20 >>ofavDD.sh

sbatch --mem=200GB -o ofavDD.o%j -e ofavDD.e%j --mail-user=studivanms@gmail.com --mail-type=ALL --partition=longq7 ofavDD.sh

module load R/3.6.1
echo '#!/bin/bash' >RQC.sh
echo Rscript ~/bin/plotQC.R prefix=dd >>RQC.sh
echo gzip -9 dd.counts >>RQC.sh
sbatch -e RQC.e%j -o RQC.o%j --dependency=afterok:460550 --mem=200GB RQC.sh

cat quality.txt

# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ, -minIndDepth and -minInd filters in subsequent ANGSD runs
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/resist/ANGSD/\*.pdf .

# Clones

module load angsd-0.933-gcc-9.2.0-65d64pp

# change the -minInd flag to 80% of the number of samples you have
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 292 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

echo '#!/bin/bash' > ofavClones.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out ofavClones >>ofavClones.sh

sbatch --mem=200GB -o ofavClones.o%j -e ofavClones.e%j --mail-type=ALL --mail-user=studivanms@gmail.com --partition=longq7 ofavClones.sh

# scp ibs matrix and list of bam files to local machine
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/resist/ANGSD/ofavClones.ibsMat .
scp mstudiva@koko-login.hpc.fau.edu:~/resist/ANGSD/bamsClones .


#------------------------------
# Next, run 2bRAD clones.R script
