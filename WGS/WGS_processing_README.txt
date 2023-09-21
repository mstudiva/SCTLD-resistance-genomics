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
## Unzipping reads

> unzip
for F in *.gz; do
echo "gunzip $F" >> unzip;
done

module load launcher
launcher_creator.py -j unzip -n unzip -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch --mem=200GB unzip.slurm

# Count raw reads
echo '#!/bin/bash' >rawReads.sh
echo readCounts.sh -e .fastq -o resistRaw >>rawReads.sh
sbatch -o rawReads.o%j -e rawReads.e%j rawReads.sh --mail-type=ALL --mail-user=studivanms@gmail.com


#------------------------------
## Creating conda environments for specialized modules

# uncomment and run below if you don't have conda set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n trimgalore parallel cutadapt fastqc
conda activate trimgalore # this is specifically for cutadapt, which doesn't play well with other KoKo modules


#------------------------------
## Removing adaptors and low quality reads

conda activate trimgalore

# does not work with launcher_creator, consider breaking up script and running multiple jobs
chmod +x trim.sh
sbatch -o trim.o%j -e trim.e%j --mem=200GB trim.sh

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


#------------------------------
## Genotyping

mkdir ../ANGSD
cd ../ANGSD
mv ../mappedReads/*.bam* .

ls *bam >bamsClones

module load angsd-0.933-gcc-9.2.0-65d64pp

export FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 2260 -minInd 113"
export TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

echo '#!/bin/bash' >ofavDD.sh
echo angsd -b bamsClones -GL 1 $FILTERS $TODO -P 1 -out dd >>ofavDD.sh

sbatch --mem=200GB -o ofavDD.o%j -e ofavDD.e%j --mail-user=studivanms@gmail.com --mail-type=ALL ofavDD.sh

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
