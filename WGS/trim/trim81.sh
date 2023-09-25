#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s163_R1_001.fastq.gz s163_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s164_R1_001.fastq.gz s164_R2_001.fastq.gz --fastqc -o trimmedReads/
