#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s087_R1_001.fastq.gz s087_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s088_R1_001.fastq.gz s088_R2_001.fastq.gz --fastqc -o trimmedReads/
