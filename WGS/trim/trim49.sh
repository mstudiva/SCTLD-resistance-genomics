#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s099_R1_001.fastq.gz s099_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s100_R1_001.fastq.gz s100_R2_001.fastq.gz --fastqc -o trimmedReads/
