#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s097_R1_001.fastq.gz s097_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s098_R1_001.fastq.gz s098_R2_001.fastq.gz --fastqc -o trimmedReads/
