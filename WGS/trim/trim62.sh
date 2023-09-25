#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s125_R1_001.fastq.gz s125_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s126_R1_001.fastq.gz s126_R2_001.fastq.gz --fastqc -o trimmedReads/
