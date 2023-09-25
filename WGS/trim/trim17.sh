#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s035_R1_001.fastq.gz s035_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s036_R1_001.fastq.gz s036_R2_001.fastq.gz --fastqc -o trimmedReads/
