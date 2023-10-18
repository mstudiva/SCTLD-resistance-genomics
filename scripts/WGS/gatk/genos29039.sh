#!/bin/bash
conda activate GATKenv
gatk --java-options -Xmx4g HaplotypeCaller    -R /mnt/beegfs/home/mstudiva/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa    -ERC GVCF    -I s153.bt2.bam.rg   -O s153.bt2.bam.rg.vcf
