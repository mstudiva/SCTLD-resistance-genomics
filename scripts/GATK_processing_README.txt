## Genome Analysis ToolKit (GATK) pipeline, version December 6, 2023
# Created by Michael Studivan (studivanms@gmail.com) based on GATK best practices
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

# NOTE: I have had issues getting some of the GATK tools to work on KoKo due to java memory errors, which could not be fixed with high-memory (1 Tb) nodes
# I transitioned to running GATK jobs on a local machine, which may take a while, but at least ran successfully


#------------------------------
## Hard Genotyping (GATK; 2bRAD and WGS together)

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
sbatch --partition=longq7 -o genos2.o%j -e genos2.e%j genos2.sh # run sbatch command with all the other versions of your script
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

# Create a temp directory in your scratch directory
mkdir ~/scratch/tmp

conda activate GATKenv
# To determine max heap size (i.e., memory) allocation for java
java -XX:+PrintFlagsFinal -version | grep HeapSize
# Rule of thumb to specify ~80% of max value in -Xmx flag below

# Combining single vcf files into a genomics database
echo '#!/bin/bash' > vcfs.sh
echo 'conda activate GATKenv' >> vcfs.sh
echo "gatk --java-options "-Xmx12g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ofav_database \
       --batch-size 50 \
       -L intervals.list \
       --sample-name-map vcfs.list \
       --tmp-dir /mnt/beegfs/home/mstudiva/scratch/tmp" >> vcfs.sh
sbatch --partition=longq7 -o vcfs.o%j -e vcfs.e%j vcfs.sh -c epyc7702 --mem=0 # -c epyc7702 --mem=0 specifies a node with 1Tb memory, and allows use of all the memory

# For some reason, I cannot get gatk to run on KoKo; keep getting out of memory errors
# But it is working on my local machine
# Follow the instructions to install GATK4 and all dependencies locally in ~/bin/: https://github.com/broadinstitute/gatk
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ofav_database \
       --batch-size 50 \
       -L intervals.list \
       --sample-name-map vcfs.list \
       --tmp-dir /Volumes/tmp # empty hard drive for temp files
# This took about 2 weeks on my local machine

# To add additional samples to the genomicsdb
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
       GenomicsDBImport \
       --genomicsdb-update-workspace-path ofav_database \
       --batch-size 50 \
       -L intervals.list \
       --sample-name-map vcfs2.list \
       --tmp-dir /Volumes/tmp # empty hard drive for temp files


#------------------------------
## Joint genotyping

# scp genome assembly and index files to local machine
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa .
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa.fai .
scp mstudiva@koko-login.hpc.fau.edu:~/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.dict .

# joint SNP calling across all samples in genomics db
~/bin/gatk-4.4.0.0/gatk --java-options "-Xmx32g" \
      GenotypeGVCFs \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V gendb://ofav_database \
      -O ofav.vcf.gz
# This took a week on my local machine


#------------------------------
## Filtering variants (finding the best-quality SNPs and indels for your dataset)

# GATK recommends using Variant Quality Score Recalibration (VQSR), but this tool requires 'truth' and 'training' datasets for the machine learning algorithm
# These resources do not exist for non-model organisms, so we will use hard-filtering instead (applying a standardized threshold for good SNPs)

# First separating SNPs from indels
~/bin/gatk-4.4.0.0/gatk SelectVariants \
      --variant ofav.vcf.gz \
      --select-type SNP \
      --output ofav_snp.vcf.gz

# for use on KoKo
echo 'conda activate GATKenv' > snps.sh
echo "gatk SelectVariants \
      --variant ofav.vcf.gz \
      --select-type SNP \
      --output ofav_snp.vcf.gz" >> snps.sh
launcher_creator.py -j snps.sh -n snps -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch snps.slurm
# This took nearly 3x as long on KoKo as it did on my local machine

# ~/bin/gatk-4.4.0.0/gatk SelectVariants \
      --variant ofav.vcf.gz \
      --select-type INDEL \
      --output ofav_indel.vcf.gz

# SNP filtering
# Filters based on GATK recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
~/bin/gatk-4.4.0.0/gatk VariantFiltration \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V ofav_snp.vcf.gz \
      -O ofav_snp_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# Make sure the number of variants is the same between the original and filtered vcf files (should be the same)
zgrep -v "^#" ofav_snp.vcf.gz | wc -l                                                  # 24847138
zgrep -v "^#" ofav_snp_filtered.vcf.gz | wc -l                                         # 24847138
# But when we sort by variants passing filter
zgrep -v "^#" ofav_snp_filtered.vcf.gz | cut -f 7 | sort | uniq -c > snp_filtered.txt  # 10653428 PASS

# Indel filtering
# Filters based on GATK recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
# ~/bin/gatk-4.4.0.0/gatk VariantFiltration \
      -R Orbicella_faveolata_gen_17.scaffolds.fa \
      -V ofav_indel.vcf.gz \
      -O ofav_indel_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "FS > 200.0" --filter-name "FS200" \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

# Make sure the number of variants is the same between the original and filtered vcf files (should be the same)
# zgrep -v "^#" ofav_indel.vcf.gz | wc -l                                                     # 3340789
# zgrep -v "^#" ofav_indel_filtered.vcf.gz | wc -l                                            # 3340789
# But when we sort by variants passing filter
# zgrep -v "^#" ofav_indel_filtered.vcf.gz | cut -f 7 | sort | uniq -c > indel_filtered.txt   # 3249070 PASS


#------------------------------
## Assessing raw vs filtered variants

# filtered SNPs (includes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_snp_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_snp_filtered.table \
     --show-filtered \
     -F FILTER

# filtered indels (includes those that failed filters)
# ~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_indel_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_indel_filtered.table \
     --show-filtered \
     -F FILTER

# passing SNPs (excludes those that failed filters)
~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_snp_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_snp_passing.table

# passing indels (excludes those that failed filters)
# ~/bin/gatk-4.4.0.0/gatk VariantsToTable \
     -R Orbicella_faveolata_gen_17.scaffolds.fa \
     -V ofav_indel_filtered.vcf.gz \
     -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
     -O ofav_indel_passing.table

# Now, use the R script GATK_filter.R to visualize the results of the filtering
# Once done, either adjust the filtering thresholds above based on your data, or proceed below


#------------------------------
## Merging SNPs and indels

# java -jar ~/bin/picard/build/libs/picard.jar SortVcf \
     I=ofav_snp_filtered.vcf.gz \
     I=ofav_indel_filtered.vcf.gz \
     O=ofav_filtered.vcf.gz

# Creating a vcf of just variants passing filter
# ~/bin/gatk-4.4.0.0/gatk SelectVariants \
     --variant ofav_filtered.vcf.gz \
     --exclude-filtered \
     --output ofav_passing.vcf.gz

# Creating a vcf of just SNPs passing filter
~/bin/gatk-4.4.0.0/gatk SelectVariants \
    --variant ofav_snp_filtered.vcf.gz \
    --exclude-filtered \
    --output ofav_snp_passing.vcf.gz


#------------------------------
## Creating identical by state (IBS) matrix for phylogenetic trees and further analysis
# Thanks to ChatGPT for this bit!

# Install dependencies
brew install bcftools
brew install vcftools
# Download plink from https://www.cog-genomics.org/plink/1.9/, put it in your working directory, then right click and open to allow run access

# converting any multiallelic variants into biallelic variants
# bcftools norm -m-any ofav_passing.vcf.gz -o ofav_passing_split.vcf.gz

# convert vcf to plink data format with vcftools
# ./plink --vcf ofav_passing_split.vcf.gz --make-bed --out ofav_passing --allow-extra-chr
./plink --vcf ofav_snp_passing.vcf.gz --make-bed --out ofav_snp_passing --allow-extra-chr --double-id

# confirm that sample IDs match between the vcf and plink files
bcftools query -l ofav_passing.vcf.gz # prints all names
bcftools query -l ofav_passing.vcf.gz | wc -l # counts the number, in this case, 365
# run plink_sampleIDs.py for plink sample IDs
python plink_sampleIDs.py # prints all names
python plink_sampleIDs.py | wc -l # 365

# calculate IBS matrix using plink2
# ./plink --bfile ofav_passing --distance ibs square --out ofav_ibsmatrix --allow-extra-chr
./plink --bfile ofav_snp_passing --distance ibs square --out ofav_snp_ibsmatrix --allow-extra-chr
