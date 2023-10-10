#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s013_R1_001.fastq.gz s013_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s014_R1_001.fastq.gz s014_R2_001.fastq.gz --fastqc -o trimmedReads/
