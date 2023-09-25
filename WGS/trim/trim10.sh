#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s021_R1_001.fastq.gz s021_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s022_R1_001.fastq.gz s022_R2_001.fastq.gz --fastqc -o trimmedReads/
