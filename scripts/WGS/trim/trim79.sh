#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s159_R1_001.fastq.gz s159_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s160_R1_001.fastq.gz s160_R2_001.fastq.gz --fastqc -o trimmedReads/
