#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s057_R1_001.fastq.gz s057_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s058_R1_001.fastq.gz s058_R2_001.fastq.gz --fastqc -o trimmedReads/
