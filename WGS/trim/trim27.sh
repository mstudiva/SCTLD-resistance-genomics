#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s055_R1_001.fastq.gz s055_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s056_R1_001.fastq.gz s056_R2_001.fastq.gz --fastqc -o trimmedReads/
