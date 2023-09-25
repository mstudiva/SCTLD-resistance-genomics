#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s069_R1_001.fastq.gz s069_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s070_R1_001.fastq.gz s070_R2_001.fastq.gz --fastqc -o trimmedReads/
