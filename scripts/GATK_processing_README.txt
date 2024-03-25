## Genome Analysis ToolKit (GATK) pipeline, version March 25, 2024
# Created by Michael Studivan (studivanms@gmail.com) based on GATK best practices
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

# NOTE: I have had issues getting some of the GATK tools to work on KoKo due to java memory errors, which could not be fixed with high-memory (1 Tb) nodes
# I transitioned to running GATK jobs on a local machine, which may take a while, but at least ran successfully


#------------------------------
## Hard Genotyping (GATK)

# First, some housekeeping, package installation, and data formatting

mkdir project/directory/GATK
mv project/directory/2bRAD/mappedReads/*.bam* .
mv project/directory/WGS/mappedReads/*.bam* .

# Need to add a read group header to each file before GATK can run
# Read groups specify which sample a read is assigned to in multiplexed samples
# on local machine, download bam files to local directory
git clone https://github.com/djhshih/rgsam.git
cd rgsam
make
make install

# Extracts read group information from bam files using rgsam and puts it in a text file
cd ~/where/bams/are/
for F in *.bam; do
samtools view $F | rgsam collect -s $F -o $F.rg.txt;
done

# Now adding header into each respective bam file
for F in *.bam; do
samtools view -h $F |
  rgsam tag -r $F.rg.txt |
  samtools view -b - > $F.rg;
done

# Double check that read groups are now in your bams
samtools view -H Sample.bt2.bam.rg | grep '@RG'
# All good? Now scp back to KoKo
scp *bam* mstudiva@koko-login.hpc.fau.edu:~/resist/GATK/

# Move all original .bam and .bam.bai files elsewhere
mv *.bam ../
mv *.bai ../

# Now remaking index files for each .rg
module load samtools-1.10-gcc-8.3.0-khgksad
echo '#!/bin/bash' > index.sh
for F in *.rg; do
echo "samtools index $F" >> index.sh;
done
launcher_creator.py -j index.sh -n index -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch index.slurm


#------------------------------
## Individual variant calls

conda activate GATKenv
export GENOME_REF=~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa

echo '#!/bin/bash' > genos.sh
echo 'conda activate GATKenv' >> genos.sh
for F in *.rg; do
echo "gatk --java-options "-Xmx4g" HaplotypeCaller \
   -R $GENOME_REF \
   -ERC GVCF \
   -I $F\
   -O $F.vcf" >> genos.sh;
   done
launcher_creator.py -j genos.sh -n genos -q mediumq7 -t 24:00:00 -N 16 -e studivanms@gmail.com
sbatch genos.slurm

# Some files may not be done after 24h - if so:
# on a local machine
gsplit -l 3 -d --additional-suffix=.sh genos2.sh genos2

# does not work with launcher_creator, consider breaking up script and running multiple jobs
chmod +x *.sh
sbatch --partition=longq7 -o genos2.o%j -e genos2.e%j --mail-type=ALL --mail-user=studivanms@gmail.com genos2.sh # run sbatch command with all the other versions of your script
# some large files took several days to complete on KoKo

# Some housekeeping
mkdir mappedReads
mv *.rg *.bai *.txt mappedReads/
cd mappedReads/
zipper.py -a -9 -f rg --launcher -e studivanms@gmail.com
sbatch zip.slurm


#------------------------------
## Creating genomics database

# Create a sample/vcf lookup tab-delimited file 'GenomicsDBImport' on your local machine and scp it to KoKo
scp vcfs.list mstudiva@koko-login.hpc.fau.edu:~/resist/GATK/

# Create a lookup table of genome scaffolds
cd ~/db/ofavgenome/
grep -e ">" Orbicella_faveolata_gen_17.scaffolds.fa | awk 'sub(/^>/, "")' > intervals.list
mv intervals.list ~/resist/GATK/intervals.list
cd ~/resist/GATK/

conda activate GATKenv
# To determine max heap size (i.e., memory) allocation for java
java -XX:+PrintFlagsFinal -version | grep HeapSize
# Rule of thumb to specify ~80% of max value in -Xmx flag below

# Combining single vcf files into a genomics database
echo '#!/bin/bash' > vcfs.sh
echo 'conda activate GATKenv' >> vcfs.sh
echo "gatk --java-options "-Xmx12g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path wgs_database \
       --batch-size 50 \
       -L intervals.list \
       --sample-name-map vcfs_wgs.list \
       --tmp-dir /mnt/beegfs/home/mstudiva/scratch/tmp" >> vcfs.sh
sbatch --partition=longq7 -o vcfs.o%j -e vcfs.e%j -c epyc7702 --mem=0 --mail-type=ALL --mail-user=studivanms@gmail.com vcfs.sh
# -c epyc7702 --mem=0 specifies a node with 1Tb memory, and allows use of all the memory

# For some reason, I cannot get gatk to run on KoKo; keep getting out of memory errors
# But it is working on my local machine
# Follow the instructions to install GATK4 and all dependencies locally in ~/bin/: https://github.com/broadinstitute/gatk

# first, scp genome assembly and index files to local machine
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa .
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa.fai .
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.dict .

# Creating genomicsdb on local machine for 2bRAD samples
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path 2brad_database \
       --batch-size 50 \
       -L ../intervals.list \
       --sample-name-map ../vcfs_2brad.list \
       -R ../Orbicella_faveolata_gen_17.scaffolds.fa \
       --tmp-dir /Volumes/tmp

       # Creating genomicsdb on local machine for WGS samples
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
      GenomicsDBImport \
      --genomicsdb-workspace-path wgs_database \
      --batch-size 50 \
      -L ../intervals.list \
      --sample-name-map ../vcfs_wgs.list \
      -R ../Orbicella_faveolata_gen_17.scaffolds.fa \
      --tmp-dir /Volumes/tmp
# This took about 3 weeks on my local machine

# To add additional samples to the genomics db
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
       GenomicsDBImport \
       --genomicsdb-update-workspace-path wgs_database \
       --batch-size 50 \
       -L intervals.list \
       --sample-name-map vcfs2_wgs.list \
       -R ../Orbicella_faveolata_gen_17.scaffolds.fa \
       --tmp-dir /Volumes/tmp # empty hard drive for temp files


#------------------------------
## Joint genotyping

# joint SNP calling for 2bRAD genomics db
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
      GenotypeGVCFs \
      -R ../Orbicella_faveolata_gen_17.scaffolds.fa \
      -V gendb://2brad_database \
      -O ../ofav_2brad.vcf.gz

# joint SNP calling for WGS genomics db
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
      GenotypeGVCFs \
      -R ../Orbicella_faveolata_gen_17.scaffolds.fa \
      -V gendb://wgs_database \
      -O ../ofav_wgs.vcf.gz


#------------------------------
## Filtering variants (finding the best-quality SNPs and indels for your dataset)

# GATK recommends using Variant Quality Score Recalibration (VQSR), but this tool requires 'truth' and 'training' datasets for the machine learning algorithm
# These resources do not exist for non-model organisms, so we will use hard-filtering instead (applying a standardized threshold for good SNPs)

# Selecting SNPs from 2bRAD samples
~/bin/gatk-4.4.0.0/gatk SelectVariants \
      --variant ofav_2brad.vcf.gz \
      --select-type SNP \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      --output ofav_2brad_snp.vcf.gz

# SNP filtering for 2bRAD samples
# Filters based on GATK recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
~/bin/gatk-4.4.0.0/gatk VariantFiltration \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V ofav_2brad_snp.vcf.gz \
      -O ofav_2brad_snp_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# Make sure the number of variants is the same between the original and filtered vcf files (should be the same)
zgrep -v "^#" ofav_2brad_snp.vcf.gz | wc -l                                                  # 61812
zgrep -v "^#" ofav_2brad_snp_filtered.vcf.gz | wc -l                                         # 61812
# But when we sort by variants passing filter
zgrep -v "^#" ofav_2brad_snp_filtered.vcf.gz | cut -f 7 | sort | uniq -c > 2brad_snp_filtered.txt  # 29767 PASS

# Separating SNPs from indels for WGS samples
~/bin/gatk-4.4.0.0/gatk SelectVariants \
      --variant ofav_wgs.vcf.gz \
      --select-type SNP \
      --output ofav_wgs_snp.vcf.gz

~/bin/gatk-4.4.0.0/gatk SelectVariants \
      --variant ofav_wgs.vcf.gz \
      --select-type INDEL \
      --output ofav_wgs_indel.vcf.gz

# SNP filtering for WGS samples
~/bin/gatk-4.4.0.0/gatk VariantFiltration \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V ofav_wgs_snp.vcf.gz \
      -O ofav_wgs_snp_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# Make sure the number of variants is the same between the original and filtered vcf files (should be the same)
zgrep -v "^#" ofav_wgs_snp.vcf.gz | wc -l                                                  # 24097804
zgrep -v "^#" ofav_wgs_snp_filtered.vcf.gz | wc -l                                         # 24097804
# But when we sort by variants passing filter
zgrep -v "^#" ofav_wgs_snp_filtered.vcf.gz | cut -f 7 | sort | uniq -c > wgs_snp_filtered.txt  # 10366466 PASS

# Indel filtering for WGS samples only (2bRAD samples don't contain indels)
# Filters based on GATK recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
~/bin/gatk-4.4.0.0/gatk VariantFiltration \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V ofav_wgs_indel.vcf.gz \
      -O ofav_wgs_indel_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "FS > 200.0" --filter-name "FS200" \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

# Make sure the number of variants is the same between the original and filtered vcf files (should be the same)
zgrep -v "^#" ofav_wgs_indel.vcf.gz | wc -l                                                     # 3238598
zgrep -v "^#" ofav_wgs_indel_filtered.vcf.gz | wc -l                                            # 3238598
# But when we sort by variants passing filter
zgrep -v "^#" ofav_wgs_indel_filtered.vcf.gz | cut -f 7 | sort | uniq -c > wgs_indel_filtered.txt   # 3145092 PASS


#------------------------------
## Assessing raw vs filtered variants

# filtered SNPs from 2bRAD samples (includes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
    -R Orbicella_faveolata_gen_17.scaffolds.fa \
    -V ofav_2brad_snp_filtered.vcf.gz \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
    -O ofav_2brad_snp_filtered.table \
    --show-filtered \
    -F FILTER

# passing SNPs from 2bRAD samples (excludes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
    -R Orbicella_faveolata_gen_17.scaffolds.fa \
    -V ofav_2brad_snp_filtered.vcf.gz \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
    -O ofav_2brad_snp_passing.table

# filtered SNPs from WGS samples (includes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_wgs_snp_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_wgs_snp_filtered.table \
     --show-filtered \
     -F FILTER

# passing SNPs from WGS samples (excludes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
    -R Orbicella_faveolata_gen_17.scaffolds.fa \
    -V ofav_snp_filtered.vcf.gz \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
    -O ofav_snp_passing.table

# filtered indels from WGS samples (includes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_wgs_indel_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_wgs_indel_filtered.table \
     --show-filtered \
     -F FILTER

# passing indels from WGS samples (excludes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_wgs_indel_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_wgs_indel_passing.table

# Now, use the R script GATK_filter.R to visualize the results of the filtering
# Once done, either adjust the filtering thresholds above based on your data and repeat filtering/QAQC, or proceed below


#------------------------------
## Exporting passing SNPs and indels

# Creating a vcf of 2bRAD SNPs passing filter
~/bin/gatk-4.4.0.0/gatk SelectVariants \
     --variant ofav_2brad_snp_filtered.vcf.gz \
     --exclude-filtered \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     --output ofav_2brad_snp_passing.vcf.gz

# Creating a vcf of WGS SNPs passing filter
~/bin/gatk-4.4.0.0/gatk SelectVariants \
     --variant ofav_wgs_snp_filtered.vcf.gz \
     --exclude-filtered \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     --output ofav_wgs_snp_passing.vcf.gz

# Creating a vcf of WGS indels passing filter
~/bin/gatk-4.4.0.0/gatk SelectVariants \
    --variant ofav_wgs_indel_filtered.vcf.gz \
    --exclude-filtered \
    -R Orbicella_faveolata_gen_17.scaffolds.fa \
    --output ofav_wgs_indel_passing.vcf.gz

# Now, use the R script GATK_clones.R to visualize genotype relatedness among samples
# and identify multi-locus genotypes


#------------------------------
## GATK_clones.R on KoKo

# The R package vcfR in GATK_clones.R is running out of RAM processing the full WGS SNP vcf
# Let's try on KoKo

# Start by creating a conda environment for R
conda create -n R r r-essentials
conda activate R
R # to activate R

# Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("vcfR")
BiocManager::install("VariantAnnotation")
if (!require("pacman")) install.packages("pacman")
pacman::p_load("dendextend", "ggdendro", "tidyverse")
quit() # to exit R session, answer 'no' to not save workspace

# scp GATK_clones_koko.R script, ofav_wgs_snp_passing.vcf.gz, bams_wgs.csv, and techReps_wgs.csv to KoKo

# Now create job script to run R scipt
echo '#!/bin/bash' >R.sh
echo 'conda activate R' >>R.sh
echo Rscript GATK_clones_koko.R >>R.sh
sbatch --partition=shortq7 -o R.o%j -e R.e%j --constraint="epyc7702" --mem=0 --mail-type=ALL --mail-user=studivanms@gmail.com R.sh
# --constraint="epyc7702" --mem=0 specifies a node with 1Tb memory, and allows use of all the memory

# Once completed, scp the dissimilarity matrix (ofav_wgs_snp_passing_dist.csv) to your local machine and determine the highest dissimilarity value between technical replicates
# Then edit the mlg.filter line to add in a rounded-up estimate of the dissimilarity as the threshold, and rerun the job script
# Once completed a second time, scp files (ofav_wgs_snp_passing.gds, mlgFilter_gatk_wgs.pdf, ofav_gatk_wgs_mlg.csv) to your local machine and continue with GATK_clones.R to generate dendrogram
