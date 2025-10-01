## 2bRAD/WGS reads processing pipeline, version November 29, 2023
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

# Trim Galore (a wrapper of cutadapt and fastQC)
git clone https://github.com/FelixKrueger/TrimGalore.git
mv TrimGalore/* .
rm -rf TrimGalore/

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py


#------------------------------
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load anaconda3-2021.05-gcc-9.4.0-llhdqho
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n cutadapt cutadapt
conda activate cutadapt # this is specifically for cutadapt, which doesn't play well with other KoKo modules

conda create -n GATKenv gatk4 picard openjdk
conda activate GATKenv # this is specifically for GATK, which requires a different version of java than KoKo


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

# Count raw reads
echo '#!/bin/bash' >rawReads.sh
echo readCounts.sh -e gz -o Raw >>rawReads.sh
sbatch -o rawReads.o%j -e rawReads.e%j --mail-type=ALL --mail-user=studivanms@gmail.com rawReads.sh
# scp RawReadCounts to local machine


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
echo '#!/bin/bash' >rawReads.sh
echo readCounts.sh -e gz -o Raw >>rawReads.sh
sbatch -o rawReads.o%j -e rawReads.e%j --mail-type=ALL --mail-user=studivanms@gmail.com rawReads.sh
# scp RawReadCounts to local machine


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
echo 'conda activate cutadapt' >> trimse.sh
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file/.tr0/}.trimlog.txt" >> trimse.sh;
done

# Non-parallel job submission
sbatch -o trimse.o%j -e trimse.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse.sh
sbatch -o trimse2.o%j -e trimse2.e%j --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com trimse2.sh

conda deactivate

# Do we have the correct number of files?
ls -l *.trim | wc -l

# Counting the trimmed reads
echo '#!/bin/bash' >cleanReads
echo readCounts.sh -e trim -o Filt >>cleanReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com cleanReads
# scp FiltReadCounts to local machine


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

conda activate cutadapt
module load fastqc/v0.11.9

# does not work with launcher_creator, consider breaking up script and running multiple jobs
chmod +x *.sh
sbatch -o trim.o%j -e trim.e%j --mail-type=ALL --mail-user=studivanms@gmail.com trim.sh # run sbatch command with all the other versions of your script

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
echo '#!/bin/bash' >cleanReads
echo readCounts.sh -e gz -o Filt >>cleanReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com cleanReads
cat FiltReadCounts
# scp FiltReadCounts to local machine


#------------------------------
## Genome formatting

cd ~/db/

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load samtools-1.10-gcc-8.3.0-khgksad

echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build Orbicella_faveolata_gen_17.scaffolds.fa OfaveolataGenome >>genomeBuild.sh
echo samtools faidx Orbicella_faveolata_gen_17.scaffolds.fa >>genomeBuild.sh
sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mail-type=ALL --mail-user=studivanms@gmail.com genomeBuild.sh

# For GATK (Hard call genotyping) only
conda activate GATKenv
cd ~/bin/
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar

cd ~/db/ofavgenome/
java -jar ~/bin/picard/build/libs/picard.jar CreateSequenceDictionary R=Orbicella_faveolata_gen_17.scaffolds.fa  O=Orbicella_faveolata_gen_17.scaffolds.dict

# Concatenated symbiont genomes
mkdir ~/db/symGenomes
# Using concatenated Symbiodiniaceae reference from Ryan Eckert GitHub
scp symbConcatGenome.fasta mstudiva@koko-login.hpc.fau.edu:~/db/symGenome/

echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build symbConcatGenome.fasta symbConcatGenome >>genomeBuild.sh
echo samtools faidx symbConcatGenome.fasta >>genomeBuild.sh
sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mail-type=ALL --mail-user=studivanms@gmail.com genomeBuild.sh


#------------------------------
## Host alignment (2bRAD)

cd ~/project/directory/2bRAD/filteredReads
mkdir symbionts

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
2bRAD_bowtie2_launcher.py -g ~/db/ofavgenome/OfaveolataGenome -f .trim --split -u un -a al --undir symbionts --launcher -e studivanms@gmail.com
sbatch --mem=200GB maps.slurm

ls *.sam | wc -l

echo '#!/bin/bash' >mappedReads
echo readCounts.sh -e al -o Host >>mappedReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedReads
# scp HostReadCounts to local machine

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
## Host alignment (WGS)

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

echo '#!/bin/bash' >mappedReads
echo readCounts.sh -e al -o Host >>mappedReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedReads
# scp HostReadCounts to local machine

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
## Symbiont Alignment (2bRAD)

cd ~/project/directory/filteredReads/symbionts
# if your samples are gzipped:
zipper.py -a -9 -f gz --gunzip --launcher -e studivanms@gmail.com
sbatch zip.slurm

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5

# aligning reads to concatenated symbiont reference; separating out unmapped reads to alignment rates calculations
mkdir trash
SYMGENOME=~/db/symGenomes/symbConcatGenome
2bRAD_bowtie2_launcher.py -g $SYMGENOME -f .un -n zooxMaps --split -u junk -a zoox --undir trash --launcher -e studivanms@gmail.com
sbatch zooxMaps.slurm

# some housekeeping
mkdir ../../mappedReads/symbionts
mv *.sam ../../mappedReads/symbionts
mv *.zoox ../../mappedReads/symbionts
cd ../../mappedReads/symbionts

# Counting the mapped zoox reads
# calculate mapping efficiency from these values compared to trimmed reads in Excel
echo '#!/bin/bash' >mappedZooxReads
echo readCounts.sh -e zoox -o Zoox >>mappedZooxReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedZooxReads
# scp ZooxReadCounts to local machine

module load samtools-1.10-gcc-8.3.0-khgksad

# making script to generate indexed sam files
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done
launcher_creator.py -j s2b -n s2b -t 6:00:00 -N 5 -e studivanms@gmail.com -q shortq7
sbatch s2b.slurm

# counting the symbiont reads by genera
>ZooxReads
for i in *.bam; do
echo $i >>ZooxReads;
samtools idxstats $i | cut -f 1,3 >>ZooxReads;
done


#------------------------------
## Symbiont Alignment (WGS)

cd ~/project/directory/filteredReads/symbionts
# if your samples are gzipped:
zipper.py -a -9 -f gz --gunzip --launcher -e studivanms@gmail.com
sbatch zip.slurm

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5

# aligning reads to concatenated symbiont reference; separating out unmapped reads to alignment rates calculations
# on a local machine, create a bowtie job script like this:
# bowtie2 --score-min L,16,1 --local -L 16 -p 2 -x /mnt/beegfs/home/mstudiva/db/symGenomes/symbConcatGenome -1 s001.1.un -2 s001.2.un -S s001.un.bt2.sam --no-unal --al-conc ./s001.zoox --un-conc trash/s001.junk
# scp to KoKo
chmod +x bowtieZoox.sh

mkdir trash
# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
launcher_creator.py -j bowtieZoox.sh -n zooxMaps -q mediumq7 -t 24:00:00 -e studivanms@gmail.com -N 16
sbatch zooxMaps.slurm

# some housekeeping
mkdir ../../mappedReads/symbionts
mv *.sam ../../mappedReads/symbionts
mv *.zoox ../../mappedReads/symbionts
cd ../../mappedReads/symbionts

# Counting the mapped zoox reads
# calculate mapping efficiency from these values compared to trimmed reads in Excel
echo '#!/bin/bash' >mappedZooxReads
echo readCounts.sh -e 1.zoox -o Zoox >>mappedZooxReads
sbatch --mem=200GB --mail-type=ALL --mail-user=studivanms@gmail.com mappedZooxReads
# scp ZooxReadCounts to local machine

module load samtools-1.10-gcc-8.3.0-khgksad

# making script to generate indexed sam files
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done
launcher_creator.py -j s2b -n s2b -t 6:00:00 -N 10 -e studivanms@gmail.com -q shortq7
sbatch s2b.slurm

# counting the symbiont reads by genera
>ZooxReads
for i in *.bam; do
echo $i >>ZooxReads;
samtools idxstats $i | cut -f 1,3 >>ZooxReads;
done


#------------------------------
# Now proceed with GATK_processing_README or ANGSD_processing_README depending on your data type
