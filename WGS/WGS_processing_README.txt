## WGS reads processing pipeline, version September 20, 2023
# Created by Michael Studivan (studivanms@gmail.com)

#------------------------------
## Scripts

# Trim Galore (a wrapper of cutadapt and fastQC)
git clone https://github.com/FelixKrueger/TrimGalore.git
mv TrimGalore/* .
rm -rf TrimGalore/

# Makes all bin scripts executable
chmod +x *.sh *.pl *.py


#------------------------------
## Getting reads

# scp sequences to project directory on KoKo
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
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n cutadaptenv cutadapt
conda activate cutadaptenv # this is specifically for cutadapt, which doesn't play well with other KoKo modules


#------------------------------
## Removing adaptors and low quality reads

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
## Aligning to reference genome

cd ~/db/

module load bowtie2-2.3.5.1-gcc-8.3.0-63cvhw5
module load samtools-1.10-gcc-8.3.0-khgksad

echo '#!/bin/bash' >genomeBuild.sh
echo bowtie2-build Orbicella_faveolata_gen_17.scaffolds.fa OfaveolataGenome >>genomeBuild.sh
echo samtools faidx Orbicella_faveolata_gen_17.scaffolds.fa >>genomeBuild.sh

sbatch -o genomeBuild.o%j -e genomeBuild.e%j --mem=200GB genomeBuild.sh

cd ~/project/directory/filteredReads
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
## Genotyping

mkdir ../ANGSD
cd ../ANGSD
mv ../mappedReads/*.bam* .
# If working with WGS and 2bRAD samples together, copy all *.bam* files to a master directory 'ANGSD'

ls *bam >bamsClones

module load angsd-0.933-gcc-9.2.0-65d64pp

# change the -minInd flag to the number of samples you have
export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 2260 -minInd 365"
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

# change the -minInd flag to the number of samples you have
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -dosnpstat 1 -doHWE 1 -hwe_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 185 -snp_pval 1e-6 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

echo '#!/bin/bash' > ofavClones.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out ofavClones >>ofavClones.sh

sbatch --mem=200GB -o ofavClones.o%j -e ofavClones.e%j -p shortq7 --mail-type=ALL --mail-user=studivanms@gmail.com ofavClones.sh

# scp to local machine
cd local/directory/
scp mstudiva@koko-login.hpc.fau.edu:~/resist/ANGSD/ofavClones.ibsMat .
scp mstudiva@koko-login.hpc.fau.edu:~/resist/ANGSD/bamsClones .
