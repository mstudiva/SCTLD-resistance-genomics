#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s067_R1_001.fastq.gz s067_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s068_R1_001.fastq.gz s068_R2_001.fastq.gz --fastqc -o trimmedReads/
