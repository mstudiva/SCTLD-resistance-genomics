#!/bin/bash
conda activate GATKenv
gatk --java-options -Xmx4g HaplotypeCaller    -R /mnt/beegfs/home/mstudiva/db/ofavgenome/Orbicella_faveolata_gen_17.scaffolds.fa    -ERC GVCF    -I s046.bt2.bam.rg   -O s046.bt2.bam.rg.vcf
