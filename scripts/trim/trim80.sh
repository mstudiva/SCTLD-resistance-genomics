#!/bin/bash
conda activate cutadaptenv
trim_galore --paired s161_R1_001.fastq.gz s161_R2_001.fastq.gz --fastqc -o trimmedReads/
trim_galore --paired s162_R1_001.fastq.gz s162_R2_001.fastq.gz --fastqc -o trimmedReads/
