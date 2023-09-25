#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s065_R1_001.fastq.gz s065_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s066_R1_001.fastq.gz s066_R2_001.fastq.gz --fastqc -o trimmedReads/
