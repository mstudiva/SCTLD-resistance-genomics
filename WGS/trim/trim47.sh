#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s095_R1_001.fastq.gz s095_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s096_R1_001.fastq.gz s096_R2_001.fastq.gz --fastqc -o trimmedReads/
