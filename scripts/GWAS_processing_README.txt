## Genome-Wide Association Study (GWAS) pipeline, version March 13, 2024
# Created by Michael Studivan (studivanms@gmail.com) based on Vollmer et al. 2023 (Science)
https://doi.org/10.1126/science.adi3601


#------------------------------
## Summary statistics

mkdir project/directory/GWAS
cp project/directory/GATK/ofav_wgs_snp_passing.vcf.gz project/directory/GWAS

module load bcftools-1.9-gcc-8.3.0-gjzr3wl
module load vcftools-0.1.14-gcc-8.3.0-safy5vc
module load plink-1.07-gcc-8.3.0-azf4a6i

# indexes the vcf file for faster processing
echo '#!/bin/sh' > index.sh
echo 'module load bcftools-1.9-gcc-8.3.0-gjzr3wl' >> index.sh
echo 'bcftools index ofav_wgs_snp_passing.vcf.gz' >> index.sh
chmod +x *.sh
sbatch --partition=shortq7 -o index.o%j -e index.e%j index.sh

# convert vcf to PLINK data format
echo '#!/bin/sh' > plink.sh
echo 'module load vcftools-0.1.14-gcc-8.3.0-safy5vc' >> plink.sh
echo 'vcftools --gzvcf ofav_wgs_snp_passing.vcf.gz --plink --chrom-map scaffolds_rename.txt --out ofav_wgs_snp_passing' >> plink.sh
chmod +x *.sh
sbatch --partition=shortq7 --mem=200GB -o plink.o%j -e plink.e%j plink.sh

# missingness per individual and per SNP and make histograms using PLINK
echo '#!/bin/sh' > missing.sh
echo 'module load plink-1.07-gcc-8.3.0-azf4a6i' >> missing.sh
echo 'plink --file ofav_wgs_snp_passing --missing' >> missing.sh
chmod +x *.sh
sbatch --partition=shortq7 --mem=200GB -o missing.o%j -e missing.e%j missing.sh
