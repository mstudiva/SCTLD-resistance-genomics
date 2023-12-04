#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s027_R1_001.fastq.gz s027_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s028_R1_001.fastq.gz s028_R2_001.fastq.gz --fastqc -o trimmedReads/
