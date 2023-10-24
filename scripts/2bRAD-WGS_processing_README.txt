## WGS reads processing pipeline, version October 20, 2023
# Created by Ryan Eckert (reckert2017@fau.edu), modified by Michael Studivan (studivanms@gmail.com) for this project
https://ryaneckert.github.io/Stephanocoenia_FKNMS_PopGen/code/


#------------------------------
## Modules

# Add to your .bashrc
module load angsd-0.933-gcc-9.2.0-65d64pp
module load bayescan-2.1-gcc-8.3.0-7gakqmd
module load qt-5.15.2-gcc-9.2.0-zi7wcem BayeScEnv/1.1
module load bcftools-1.9-gcc-8.3.0-il4d373
module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load cdhit-4.8.1-gcc-8.3.0-bcay75d
module load htslib-1.9-gcc-8.3.0-jn7ehrc
module load kraken2-2.1.1-gcc-9.2.0-ocivj3u
module load python-3.7.4-gcc-8.3.0-3tniqr5
module load launcher
module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
module load ncbi-toolkit-22_0_0-gcc-9.2.0-jjhd2wa
module load ngsadmix-32-gcc-8.3.0-qbnwmpq
module load ngsRelate/v2
module load R/3.6.1
module load samtools-1.10-gcc-8.3.0-khgksad
module load vcftools-0.1.14-gcc-8.3.0-safy5vc


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

# Trim Galore (a wrapper of cutadapt and fastQC)
git clone https://github.com/FelixKrueger/TrimGalore.git
mv TrimGalore/* .
rm -rf TrimGalore/

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py


#------------------------------
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n cutadaptenv cutadapt
conda activate cutadaptenv # this is specifically for cutadapt, which doesn't play well with other KoKo modules

conda create -n GATKenv gatk4 picard openjdk
conda activate GATKenv # this is specifically for cutadapt, which requires a different version of java than KoKo


#------------------------------
## 2bRAD reads

# transfer directly from Illumina BaseSpace)
echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n JA23220 -o .' >> downloadReads.sh
# -n is the project name and -o is the output directory

echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'rmdir SA*' >>downloadReads.sh
#echo 'mkdir ../concatReads' >> downloadReads.sh
#echo 'cp *.gz ../concatReads' >> downloadReads.sh
#echo 'cd ../concatReads' >> downloadReads.sh
echo 'mergeReads.sh -o mergeTemp' >> downloadReads.sh
# -o is the directory to put output files in

echo 'rm *L00*' >> downloadReads.sh
echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'gunzip *.gz' >> downloadReads.sh
echo 'rmdir mergeTemp' >> downloadReads.sh

chmod +x downloadReads.sh

launcher_creator.py -b 'srun downloadReads.sh' -n downloadReads -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch downloadReads.slurm

#------------------------------
## WGS reads

# scp sequences from local machine to project directory on KoKo
scp *.gz mstudiva@koko-login.hpc.fau.edu:~/path/to/project/directory/

# FastQC of raw reads
module load fastqc/v0.11.9
> rawqc.sh
for F in *.gz; do
echo "fastqc $F" >> rawqc.sh;
done

launcher_creator.py -j rawqc.sh -n fastqc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch --mem=200GB fastqc.slurm
# when done, scp .html files to local machine

# Count raw reads (doesn't work when files are gzipped [.gz])
# echo '#!/bin/bash' >rawReads.sh
# echo readCounts.sh -e gz -o Raw >>rawReads.sh
# cat RawReadCounts
# sbatch -o rawReads.o%j -e rawReads.e%j rawReads.sh --mail-type=ALL --mail-user=studivanms@gmail.com


#------------------------------
## 2bRAD deduplication and trimming

# Deduplicates row pools into separate 3ill-BC's (1-12), using reverse complement as the ID
2bRAD_trim_launch_dedup.pl fastq > trims.sh
launcher_creator.py -j trims.sh -n trims -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB trims.slurm

# Do we have the correct number of files?
ls -l *.tr0 | wc -l

mkdir trimmedReads
srun mv *.tr0 trimmedReads/ &

# Rezips the raw fastq's for storage
zipper.py -f fastq -a -9 --launcher -e studivanms@gmail.com
sbatch --mem=200GB zip.slurm

cd ../trimmedReads

# Renames files based on two column lookup table (sampleID.csv): filename, then sample ID
srun sampleRename.py -i sampleID -f tr0

# For loop to generate a list of commands for each file
echo '#!/bin/bash' > trimse.sh
echo 'module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj' >> trimse.sh
echo 'conda activate cutadaptenv' >> trimse.sh
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file/.tr0/}.trimlog.txt" >> trimse.sh;
done

# Old-school way of submitting jobs
sbatch -o trimse.o%j -e trimse.e%j --mem=200GB trimse.sh
sbatch -o trimse2.o%j -e trimse2.e%j --mem=200GB trimse2.sh

conda deactivate

# Do we have the correct number of files?
ls -l *.trim | wc -l

# Counting the trimmed reads
echo '#!/bin/bash' >cleanReads
echo readCounts.sh -e trim -o Filt >>cleanReads
sbatch --mem=200GB cleanReads

mkdir ../filteredReads
mv *.trim ../filteredReads

# Rezips the row pools for storage
zipper.py -f tr0 -a -9 --launcher -e studivanms@gmail.com
sbatch zip.slurm

cat FiltReadCounts


#------------------------------
## WGS trimming

# on a local machine
gsplit -l 4 -d --additional-suffix=.sh trim.sh trim

# back on KoKo
mkdir filteredReads

conda activate cutadaptenv
module load fastqc/v0.11.9

# does not work with launcher_creator, consider breaking up script and running multiple jobs
chmod +x *.sh
sbatch -o trim.o%j -e trim.e%j trim.sh # run sbatch command with all the other versions of your script

conda deactivate

# Do we have the correct number of files?
cd trimmedReads/
ls -l *.gz | wc -l

# FastQC of trimmed reads
module load fastqc/v0.11.9
> trimqc.sh
for F in *.gz; do
echo "fastqc $F" >> trimqc.sh;
done

launcher_creator.py -j trimqc.sh -n fastqc -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch --mem=200GB fastqc.slurm
# when done, scp .html files to local machine

# Counting the trimmed reads (doesn't work when files are gzipped [.gz])
# echo '#!/bin/bash' >cleanReads
# echo readCounts.sh -e gz -o Filt >>cleanReads
# sbatch --mem=200GB cleanReads
# cat FiltReadCounts


#------------------------------
## Genome formatting

cd ~/db/

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load samtools-1.10-gcc-8.3.0-khgksad

echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build Orbicella_faveolata_gen_17.scaffolds.fa OfaveolataGenome >>genomeBuild.sh
echo samtools faidx Orbicella_faveolata_gen_17.scaffolds.fa >>genomeBuild.sh

sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mem=200GB genomeBuild.sh

# For GATK (Hard call genotyping) only
conda activate GATKenv
cd ~/bin/
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar

cd ~/db/ofavgenome/
java -jar ~/bin/picard/build/libs/picard.jar CreateSequenceDictionary R=Orbicella_faveolata_gen_17.scaffolds.fa  O=Orbicella_faveolata_gen_17.scaffolds.dict


#------------------------------
## 2bRAD alignment

cd ~/project/directory/2bRAD/filteredReads
mkdir symbionts

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
2bRAD_bowtie2_launcher.py -g ~/db/ofavgenome/OfaveolataGenome -f .trim --split -u un -a al --undir symbionts --launcher -e studivanms@gmail.com
sbatch --mem=200GB maps.slurm

ls *.sam | wc -l

ls *al | cut -d '.' -f 1 >align1
grep "% overall" maps.e* | cut -d ' ' -f 1 >align2
paste <(awk -F' ' '{print $1}' align1) <(awk -F' ' '{print $1}' align2) >alignmentRates
rm align1 align2
less alignmentRates

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

launcher_creator.py -j s2b -n s2b -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB s2b.slurm

ls *bam | wc -l

zipper.py -a -9 -f sam --launcher -e studivanms@gmail.com
sbatch zip.slurm

zipper.py -a -9 -f trim --launcher -e studivanms@gmail.com
sbatch zip.slurm

mkdir ../mappedReads
mv *.sam.gz ../mappedReads
mv *.al ../mappedReads


#------------------------------
## WGS alignment

cd ~/project/directory/WGS/filteredReads
mkdir symbionts

# on a local machine, create a bowtie job script like this:
# bowtie2 --score-min L,16,1 --local -L 16 -x /mnt/beegfs/home/mstudiva/db/ofavgenome/OfaveolataGenome -1 s001_R1_001_val_1.fq.gz -2 s001_R2_001_val_2.fq.gz -S s001.bt2.sam --no-unal --al-conc ./s001.al --un-conc symbionts/s001.un
# scp to KoKo
chmod +x bowtie.sh

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
launcher_creator.py -j bowtie.sh -n maps -q mediumq7 -t 24:00:00 -e studivanms@gmail.com -N 24
sbatch maps.slurm

ls *.sam | wc -l

ls *.1.al | cut -d '.' -f 1 >align1
grep "% overall" maps.e* | cut -d ' ' -f 1 >align2
paste <(awk -F' ' '{print $1}' align1) <(awk -F' ' '{print $1}' align2) >alignmentRates
rm align1 align2
less alignmentRates

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

launcher_creator.py -j s2b -n s2b -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB s2b.slurm

ls *bam | wc -l

zipper.py -a -9 -f sam --launcher -e studivanms@gmail.com
sbatch zip.slurm

mkdir ../mappedReads
mv *.sam.gz ../mappedReads
mv *.al ../mappedReads


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
# Next, run 2bRAD clones.R script


#------------------------------
## Hard Genotyping (GATK; 2bRAD and WGS together)

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

# Ok, let's genotype!
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

# Create a sample/vcf lookup tab-delimited file 'GenomicsDBImport' on your local machine and scp it to KoKo
scp vcfs mstudiva@koko-login.hpc.fau.edu:~/resist/GATK/

# Create a lookup table of genome scaffolds
cd ~/db/ofavgenome/
grep -e ">" Orbicella_faveolata_gen_17.scaffolds.fa | awk 'sub(/^>/, "")' > intervals
mv intervals ~/resist/GATK/intervals
cd ~/resist/GATK/

# Combining all vcf files into a genomics workspace
echo "gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ofav_database \
       --batch-size 50 \
       -L intervals \
       --sample-name-map vcfs \
       --tmp-dir=~/scratch/ \
       --reader-threads 5" > combine.sh
chmod +x combine.sh
sbatch --partition=longq7 -o combine.o%j -e combine.e%j combine.sh --cpus-per-task=5
