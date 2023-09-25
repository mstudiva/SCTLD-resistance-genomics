#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s023_R1_001.fastq.gz s023_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s024_R1_001.fastq.gz s024_R2_001.fastq.gz --fastqc -o trimmedReads/
