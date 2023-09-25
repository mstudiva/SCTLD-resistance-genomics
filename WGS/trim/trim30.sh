#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s061_R1_001.fastq.gz s061_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s062_R1_001.fastq.gz s062_R2_001.fastq.gz --fastqc -o trimmedReads/
