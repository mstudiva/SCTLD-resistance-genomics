#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s157_R1_001.fastq.gz s157_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s158_R1_001.fastq.gz s158_R2_001.fastq.gz --fastqc -o trimmedReads/
