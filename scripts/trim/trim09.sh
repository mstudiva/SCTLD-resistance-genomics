#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s019_R1_001.fastq.gz s019_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s020_R1_001.fastq.gz s020_R2_001.fastq.gz --fastqc -o trimmedReads/
