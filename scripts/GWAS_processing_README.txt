## Genome-Wide Association Study (GWAS) pipeline, version April 1, 2024
# Created by Michael Studivan (studivanms@gmail.com) based on Andries Marees' GitHub and Vollmer et al. 2023 (Science)
https://github.com/MareesAT/GWA_tutorial/
https://doi.org/10.1126/science.adi3601

# Getting some R scripts for QAQC and visualization
git clone https://github.com/MareesAT/GWA_tutorial.git
unzip GWA_tutorial/1_QC_GWAS.zip
mv GWA_tutorial/1_QC_GWAS/*.R .


#------------------------------
## Removing clones

mkdir project/directory/GWAS
cp project/directory/GATK/ofav_wgs_snp_passing.vcf.gz project/directory/GWAS

module load bcftools-1.9-gcc-8.3.0-gjzr3wl
module load vcftools-0.1.14-gcc-8.3.0-safy5vc
module load plink-1.07-gcc-8.3.0-azf4a6i

# Indexes the vcf file for faster processing
echo '#!/bin/sh' > index.sh
echo 'module load bcftools-1.9-gcc-8.3.0-gjzr3wl' >> index.sh
echo 'bcftools index ofav_wgs_snp_passing.vcf.gz' >> index.sh
sbatch --partition=shortq7 -o index.o%j -e index.e%j index.sh

# Convert vcf to PLINK data format
echo '#!/bin/sh' > plink.sh
echo 'module load vcftools-0.1.14-gcc-8.3.0-safy5vc' >> plink.sh
echo 'vcftools --gzvcf ofav_wgs_snp_passing.vcf.gz --plink --chrom-map scaffolds_rename.txt --out ofav_wgs_snp_passing' >> plink.sh
sbatch --partition=shortq7 --mem=200GB -o plink.o%j -e plink.e%j plink.sh

# Missingness per individual and per SNP, and make histograms using PLINK
echo '#!/bin/sh' > missing.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> missing.sh
echo 'plink --file ofav_wgs_snp_passing --missing' >> missing.sh
sbatch --partition=shortq7 --mem=200GB -o missing.o%j -e missing.e%j missing.sh

# scp plink.imiss and plink.lmiss to your local machine for individuals missing SNP loci and SNP loci missing individuals, respectively
# Use plink.imiss specifically to determine which clone-mates should be removed (remove the clonal genotypes with the most missing data)
# Create a single-column text file containing all the clonal samples to remove and scp it to your working directory

# Removing all but one of each of the clonal genotypes
echo '#!/bin/sh' > noclones.sh
echo 'module load vcftools-0.1.14-gcc-8.3.0-safy5vc' >> noclones.sh
echo 'vcftools --gzvcf ofav_wgs_snp_passing.vcf.gz --remove clones_remove.txt --plink --chrom-map scaffolds_rename.txt --out ofav_wgs_snp_passing_noclones' >> noclones.sh
sbatch --partition=shortq7 --mem=200GB -o noclones.o%j -e noclones.e%j noclones.sh

# Exporting the clones-removed file as .vcf.gz for later
echo '#!/bin/sh' > exportvcf.sh
echo 'module load vcftools-0.1.14-gcc-8.3.0-safy5vc' >> exportvcf.sh
echo 'vcftools --gzvcf ofav_wgs_snp_passing.vcf.gz --remove clones_remove.txt --recode --out ofav_wgs_snp_passing_noclones' >> exportvcf.sh
sbatch --partition=shortq7 --mem=200GB -o exportvcf.o%j -e exportvcf.e%j exportvcf.sh
# scp ofav_wgs_snp_passing_noclones.vcf to your local machine

bgzip -c ofav_wgs_snp_passing_noclones.recode.vcf > ofav_wgs_snp_passing_noclones.vcf.gz
tabix -p vcf ofav_wgs_snp_passing_noclones.vcf.gz


#------------------------------
## Missingness

# Generate plots to visualize the missingness distributions
conda activate R
echo '#!/bin/bash' >Rmiss.sh
echo 'conda activate R' >>Rmiss.sh
echo Rscript --no-save hist_miss.R >>Rmiss.sh
sbatch --partition=shortq7 --mem=200GB -o Rmiss.o%j -e Rmiss.e%j Rmiss.sh
# scp *miss.pdf to your local machine

# Remove SNPs with missingness >0.2
echo '#!/bin/sh' > filter1.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> filter1.sh
echo 'plink --file ofav_wgs_snp_passing_noclones --geno 0.2 --make-bed --out ofav_wgs_snp_passing_noclones_filtered1' >> filter1.sh
sbatch --partition=shortq7 --mem=200GB -o filter1.o%j -e filter1.e%j filter1.sh

# Remove individuals with missingness >0.2
echo '#!/bin/sh' > filter2.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> filter2.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered1 --mind 0.2 --make-bed --out ofav_wgs_snp_passing_noclones_filtered2' >> filter2.sh
sbatch --partition=shortq7 --mem=200GB -o filter2.o%j -e filter2.e%j filter2.sh


#------------------------------
## Minor allele frequency (MAF)

echo '#!/bin/sh' > lowmaf1.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> lowmaf1.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered2 --freq --out MAF_check' >> lowmaf1.sh
sbatch --partition=shortq7 --mem=200GB -o lowmaf1.o%j -e lowmaf1.e%j lowmaf1.sh

# Generate plots to visualize the MAF distributions
conda activate R
echo '#!/bin/bash' >Rmaf.sh
echo 'conda activate R' >>Rmaf.sh
echo Rscript --no-save MAF_check.R >>Rmaf.sh
sbatch --partition=shortq7 --mem=200GB -o Rmaf.o%j -e Rmaf.e%j Rmaf.sh
# scp MAF_distribution.pdf to your local machine

# Remove SNPs with a low minor allele frequency (MAF)
# A conventional MAF threshold for GWAS is between 0.01 - 0.05
echo '#!/bin/sh' > lowmaf2.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> lowmaf2.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered2 --maf 0.05 --make-bed --out ofav_wgs_snp_passing_noclones_filtered_lowmaf' >>lowmaf2.sh
sbatch --partition=shortq7 --mem=200GB -o lowmaf2.o%j -e lowmaf2.e%j lowmaf2.sh


#------------------------------
## Hardy-Weinberg equilibrium (HWE)

echo '#!/bin/sh' > hwe1.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> hwe1.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf --hardy' >> hwe1.sh
sbatch --partition=shortq7 --mem=200GB -o hwe1.o%j -e hwe1.e%j hwe1.sh

# Selecting SNPs with HWE p-value below 0.00001, required for one of the two plots generated by the next Rscript, allows to zoom in on strongly deviating SNPs.
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe

# Generate plots to visualize the HWE distributions
conda activate R
echo '#!/bin/bash' >Rhwe.sh
echo 'conda activate R' >>Rhwe.sh
echo Rscript --no-save hwe.R >>Rhwe.sh
sbatch --partition=shortq7 --mem=200GB -o Rhwe.o%j -e Rhwe.e%j Rhwe.sh
# scp histhwe_below_theshold.pdf to your local machine

# Remove SNPs which are not in HWE
echo '#!/bin/sh' > hwe2.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> hwe2.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf --hwe 1e-6 --hwe-all --make-bed --out ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe' >> hwe2.sh
sbatch --partition=shortq7 --mem=200GB -o hwe2.o%j -e hwe2.e%j hwe2.sh


#------------------------------
## Heterozygosity

# The parameters '50 5 0.2' stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously
echo '#!/bin/sh' > het1.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> het1.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe --indep-pairwise 50 5 0.2 --out indepSNP' >> het1.sh
sbatch --partition=shortq7 --mem=200GB -o het1.o%j -e het1.e%j het1.sh

echo '#!/bin/sh' > het2.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> het2.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe --extract indepSNP.prune.in --het --out R_check' >> het2.sh
sbatch --partition=shortq7 --mem=200GB -o het2.o%j -e het2.e%j het2.sh

# Plot of the heterozygosity rate distribution
conda activate R
echo '#!/bin/bash' >Rhet1.sh
echo 'conda activate R' >>Rhet1.sh
echo Rscript --no-save check_heterozygosity_rate.R >>Rhet1.sh
sbatch --partition=shortq7 --mem=200GB -o Rhet1.o%j -e Rhet1.e%j Rhet1.sh
# scp heterozygosity.pdf to your local machine

# Generate a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean
conda activate R
echo '#!/bin/bash' >Rhet2.sh
echo 'conda activate R' >>Rhet2.sh
echo Rscript --no-save heterozygosity_outliers_list.R >>Rhet2.sh
sbatch --partition=shortq7 --mem=200GB -o Rhet2.o%j -e Rhet2.e%j Rhet2.sh

# Generates fail-het-qc.txt
# Modify this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers
echo '#!/bin/sh' > het3.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> het3.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe --remove het_fail_ind.txt --make-bed --out ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe_het' >>het3.sh
sbatch --partition=shortq7 --mem=200GB -o het3.o%j -e het3.e%j het3.sh


#------------------------------
## Export for GWAS

echo '#!/bin/sh' > export.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> export.sh
echo 'plink --bfile ofav_wgs_snp_passing_noclones_filtered_lowmaf_hwe_het --recode --out ofav_wgs_gwas' >>export.sh
sbatch --partition=shortq7 --mem=200GB -o export.o%j -e export.e%j export.sh
# scp ofav_wgs_gwas.ped and ofav_wgs_gwas.map to your local machine

# Continue with the R script GWAS_analysis.R


#------------------------------
## GWAS_analysis.R on KoKo

# Generating new dissimilarity matrix for dendrogram plotting with susceptibility data
# scp GWAS_analysis_koko.R to your working directory

# Start with our handy dandy R conda environment
conda activate R
R # to activate R shell

# Now create job script to run R scipt
echo '#!/bin/bash' >Rgwas.sh
echo 'conda activate R' >>Rgwas.sh
echo Rscript GWAS_analysis_koko.R >>Rgwas.sh
sbatch --partition=shortq7 -o Rgwas.o%j -e Rgwas.e%j --constraint="epyc7702" --mem=0 --mail-type=ALL --mail-user=studivanms@gmail.com Rgwas.sh
# --constraint="epyc7702" --mem=0 specifies a node with 1Tb memory, and allows use of all the memory
